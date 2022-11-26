mod simulation;
mod snapshots;

use std::{f64::consts::PI, str::FromStr};

use chrono::{DateTime, Datelike, Duration, NaiveDate, Utc};
use clap::{App, Arg};
use lazy_static::lazy_static;
use nalgebra::Vector3;
use numeric_algs::symplectic::integration::{Integrator, StepSize, SuzukiIntegrator};
use simulation::SimState;
use snapshots::Snapshots;

const STEP: f64 = 60.0;
const OMEGA: f64 = 24.06570982441908 / 12.0 * PI / 86400.0;
const ANGLE_2022_01_01: f64 = 3.3271216795200544;
lazy_static! {
    static ref YEAR: Duration = Duration::seconds(365 * 24 * 3600 + 6 * 3600);
    static ref EPOCH: DateTime<Utc> = DateTime::<Utc>::from_str("2022-01-01T00:00:00Z").unwrap();
}

struct Himawari {
    axis: Vector3<f64>,
    v1: Vector3<f64>,
    v2: Vector3<f64>,
}

impl Himawari {
    fn new() -> Self {
        let obliquity = 23.45_f64.to_radians();

        let axis = Vector3::new(0.0, obliquity.sin(), obliquity.cos());
        let v1 = Vector3::new(0.0, -obliquity.cos(), obliquity.sin());
        let v2 = Vector3::new(1.0, 0.0, 0.0);

        Self { axis, v1, v2 }
    }

    fn looking_dir(&self, t: f64) -> Vector3<f64> {
        let beta = ANGLE_2022_01_01 + OMEGA * t + 140.7_f64.to_radians();
        -self.v1 * beta.cos() - self.v2 * beta.sin()
    }

    fn pos(&self, t: f64) -> Vector3<f64> {
        let r = 42171.0;
        -r * self.looking_dir(t)
    }

    fn within_frame(&self, sim: &SimState) -> bool {
        let t = (sim.time() - *EPOCH).num_nanoseconds().unwrap() as f64 / 1e9;

        let dir = self.looking_dir(t);
        let up = self.axis;
        let right = dir.cross(&up);

        let moon = sim.body_by_name("Moon").unwrap();
        let earth = sim.body_by_name("Earth").unwrap();
        let pos = self.pos(t);
        let himawari_pos = earth.pos + pos;
        let dir_to_moon = (moon.pos - himawari_pos).normalize();

        let x = dir_to_moon.dot(&right);
        let y = dir_to_moon.dot(&up);
        let z = dir_to_moon.dot(&dir);

        let tan_fov2 = 8.7_f64.to_radians().tan();

        (x / z).abs() < tan_fov2 && (y / z).abs() < tan_fov2 && z > 0.0
    }

    fn ang_to_moon(&self, sim: &SimState) -> f64 {
        let t = (sim.time() - *EPOCH).num_nanoseconds().unwrap() as f64 / 1e9;

        let moon = sim.body_by_name("Moon").unwrap();
        let earth = sim.body_by_name("Earth").unwrap();
        let pos = self.pos(t);
        let dir = self.looking_dir(t);
        let himawari_pos = earth.pos + pos;
        let dir_to_moon = (moon.pos - himawari_pos).normalize();
        dir_to_moon.dot(&dir).acos()
    }
}

fn nearest_month_start(datetime: DateTime<Utc>) -> DateTime<Utc> {
    let month = datetime.month();
    let prev_month = DateTime::<Utc>::from_utc(
        NaiveDate::from_ymd(datetime.year(), month, 1).and_hms(0, 0, 0),
        Utc,
    );
    let next_month = DateTime::<Utc>::from_utc(
        if month == 12 {
            NaiveDate::from_ymd(datetime.year() + 1, 1, 1)
        } else {
            NaiveDate::from_ymd(datetime.year(), month + 1, 1)
        }
        .and_hms(0, 0, 0),
        Utc,
    );
    if datetime - prev_month < next_month - datetime {
        prev_month
    } else {
        next_month
    }
}

fn close_to_month(datetime: DateTime<Utc>) -> bool {
    let nearest_month = nearest_month_start(datetime);
    ((datetime - nearest_month).num_nanoseconds().unwrap() as f64 / 1e9).abs() < STEP
}

fn propagate_to<I: Integrator<SimState>>(
    integrator: &mut I,
    sim: &mut SimState,
    snapshots: &mut Snapshots,
    datetime: DateTime<Utc>,
) {
    while datetime > sim.time() {
        let step_size = STEP.min((datetime - sim.time()).num_nanoseconds().unwrap() as f64 / 1e9);
        integrator.propagate_in_place(
            sim,
            SimState::position_derivative,
            SimState::momentum_derivative,
            StepSize::Step(step_size),
        );
        maybe_save_snapshot(integrator, sim, snapshots);
    }
    while datetime < sim.time() {
        let step_size = STEP.min((sim.time() - datetime).num_nanoseconds().unwrap() as f64 / 1e9);
        integrator.propagate_in_place(
            sim,
            SimState::position_derivative,
            SimState::momentum_derivative,
            StepSize::Step(-step_size),
        );
        maybe_save_snapshot(integrator, sim, snapshots);
    }
}

fn maybe_save_snapshot<I: Integrator<SimState>>(
    integrator: &mut I,
    sim: &SimState,
    snapshots: &mut Snapshots,
) {
    if close_to_month(sim.time()) {
        let mut sim_clone = sim.clone();
        let nearest_month = nearest_month_start(sim.time());
        if nearest_month != sim_clone.time() {
            propagate_to(integrator, &mut sim_clone, snapshots, nearest_month);
        }
        snapshots.insert(sim_clone);
    }
}

fn generate(start_date: DateTime<Utc>, end_date: DateTime<Utc>) {
    let mut snapshots = Snapshots::new();

    let mut integrator = SuzukiIntegrator::new(STEP);

    let mut sim = snapshots.get_closest(start_date);
    propagate_to(&mut integrator, &mut sim, &mut snapshots, start_date);

    let mut currently_visible = false;
    let himawari = Himawari::new();

    while sim.time() < end_date {
        sim.step_forwards(&mut integrator);

        maybe_save_snapshot(&mut integrator, &sim, &mut snapshots);

        let ang_to_moon = himawari.ang_to_moon(&sim);
        let obscured = ang_to_moon < 8.45_f64.to_radians();
        let within_frame = himawari.within_frame(&sim);

        let date = sim.time();

        match (within_frame, obscured, currently_visible) {
            (true, false, false) => {
                println!("Becoming visible: {}", date);
                currently_visible = true;
            }
            (true, false, true) => (),
            (true, true, false) => (),
            (true, true, true) => {
                println!("Becoming obscured: {}", date);
                currently_visible = false;
            }
            (false, false, false) => (),
            (false, false, true) => {
                println!("Leaving frame: {}\n", date);
                currently_visible = false;
            }
            (false, true, false) => (),
            (false, true, true) => {
                println!("Becoming obscured outside of the frame? {}", date);
                currently_visible = false;
            }
        }
    }
}

fn main() {
    let matches = App::new("Himawari Moon Predictor")
        .about("App predicting when the Moon should be visible in Himawari photos")
        .arg(
            Arg::with_name("start_date")
                .long("start-date")
                .short("s")
                .value_name("DATETIME")
                .help("Starting date in the format YYYY-MM-DDTHH:MM:SSZ (default 2022-01-01T00:00:00Z)")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("end_date")
                .long("end-date")
                .short("e")
                .value_name("DATETIME")
                .help("End date in the format YYYY-MM-DDTHH:MM:SSZ (default 2024-01-01T00:00:00Z)")
                .takes_value(true),
        )
        .get_matches();

    let start_date = matches
        .value_of("start_date")
        .unwrap_or("2022-01-01T00:00:00Z");
    let end_date = matches
        .value_of("end_date")
        .unwrap_or("2024-01-01T00:00:00Z");
    generate(
        DateTime::<Utc>::from_str(start_date).unwrap(),
        DateTime::<Utc>::from_str(end_date).unwrap(),
    );
}
