mod simulation;
mod snapshots;

use std::{f64::consts::PI, str::FromStr};

use chrono::{DateTime, Duration, Utc};
use nalgebra::Vector3;
use numeric_algs::symplectic::integration::{Integrator, StepSize, SuzukiIntegrator};
use simulation::SimState;
use snapshots::Snapshots;

const STEP: f64 = 60.0;
const YEAR: f64 = 365.25 * 24.0 * 3600.0;
const OMEGA: f64 = 24.06570982441908 / 12.0 * PI / 86400.0;
const ANGLE_2022_01_01: f64 = 3.3271216795200544;

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

    fn within_frame(&self, sim: &SimState, t: f64) -> bool {
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

    fn ang_to_moon(&self, sim: &SimState, t: f64) -> f64 {
        let moon = sim.body_by_name("Moon").unwrap();
        let earth = sim.body_by_name("Earth").unwrap();
        let pos = self.pos(t);
        let dir = self.looking_dir(t);
        let himawari_pos = earth.pos + pos;
        let dir_to_moon = (moon.pos - himawari_pos).normalize();
        dir_to_moon.dot(&dir).acos()
    }
}

struct SimWithTimestamp {
    timestamp: f64,
    state: SimState,
    integrator: SuzukiIntegrator,
}

impl SimWithTimestamp {
    fn new(timestamp: i64, state: SimState) -> Self {
        Self {
            timestamp: timestamp as f64,
            state,
            integrator: SuzukiIntegrator::new(STEP),
        }
    }

    fn propagate_to(&mut self, timestamp: f64) {
        while timestamp > self.timestamp {
            let step_size = STEP.min(timestamp - self.timestamp);
            self.integrator.propagate_in_place(
                &mut self.state,
                SimState::position_derivative,
                SimState::momentum_derivative,
                StepSize::Step(step_size),
            );
            self.timestamp += step_size;
        }
        while timestamp < self.timestamp {
            let step_size = STEP.min(self.timestamp - timestamp);
            self.integrator.propagate_in_place(
                &mut self.state,
                SimState::position_derivative,
                SimState::momentum_derivative,
                StepSize::Step(-step_size),
            );
            self.timestamp -= step_size;
        }
    }

    fn step_forwards(&mut self) {
        self.integrator.propagate_in_place(
            &mut self.state,
            SimState::position_derivative,
            SimState::momentum_derivative,
            StepSize::UseDefault,
        );
        self.timestamp += STEP;
    }

    fn time_since(&self, epoch: DateTime<Utc>) -> f64 {
        self.timestamp - (epoch.timestamp() as f64)
    }
}

fn main() {
    let epoch = DateTime::<Utc>::from_str("2022-01-01T00:00:00Z").unwrap();

    let snapshots = Snapshots::new();

    let (timestamp, sim) = snapshots.get_closest(epoch.timestamp());
    let mut sim = SimWithTimestamp::new(timestamp, sim);
    sim.propagate_to(epoch.timestamp() as f64);

    let mut currently_visible = false;
    let himawari = Himawari::new();

    while sim.time_since(epoch) < 2.0 * YEAR {
        sim.step_forwards();
        let time = sim.time_since(epoch);

        let ang_to_moon = himawari.ang_to_moon(&sim.state, time);
        let obscured = ang_to_moon < 8.45_f64.to_radians();
        let within_frame = himawari.within_frame(&sim.state, time);

        let date = epoch + Duration::seconds(time as i64);

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
