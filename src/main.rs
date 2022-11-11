mod simulation;

use std::{f64::consts::PI, str::FromStr};

use chrono::{DateTime, Duration, Utc};
use nalgebra::Vector3;
use numeric_algs::symplectic::integration::{Integrator, StepSize, SuzukiIntegrator};
use simulation::{Body, SimState};

const STEP: f64 = 60.0;
const YEAR: f64 = 365.25 * 24.0 * 3600.0;
const OMEGA: f64 = 24.06570982441908 / 12.0 * PI / 86400.0;

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
        let beta = 0.697374558 * PI / 12.0 + OMEGA * (t + 43200.0) + 140.0_f64.to_radians();
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

fn main() {
    let epoch = DateTime::<Utc>::from_str("2000-01-01T00:00:00Z").unwrap();

    let mut sim = SimState::new()
        .with_body(Body {
            name: "Sun".to_owned(),
            gm: 132712440041.93938,
            pos: Vector3::new(
                -1.068108951496322E+06,
                -4.177210908491462E+05,
                3.086887010002915E+04,
            ),
            vel: Vector3::new(
                9.305302656256911E-03,
                -1.283177282717393E-02,
                -1.631700118015769E-04,
            ),
            radius: 696000.0,
        })
        .with_body(Body {
            name: "Mercury".to_owned(),
            gm: 22031.86855,
            pos: Vector3::new(
                -2.212073002393702E+07,
                -6.682435921338345E+07,
                -3.461577076477692E+06,
            ),
            vel: Vector3::new(
                3.666229234452722E+01,
                -1.230266984222893E+01,
                -4.368336206255391E+00,
            ),
            radius: 2440.0,
        })
        .with_body(Body {
            name: "Venus".to_owned(),
            gm: 324858.592,
            pos: Vector3::new(
                -1.085736592234813E+08,
                -3.784241757371509E+06,
                6.190088659339075E+06,
            ),
            vel: Vector3::new(
                8.984650886248794E-01,
                -3.517203951420625E+01,
                -5.320225928762774E-01,
            ),
            radius: 6052.0,
        })
        .with_body(Body {
            name: "Earth".to_owned(),
            gm: 398600.435436,
            pos: Vector3::new(
                -2.627903751048988E+07,
                1.445101984929515E+08,
                3.025245352813601E+04,
            ),
            vel: Vector3::new(
                -2.983052803412253E+01,
                -5.220465675237847E+00,
                -1.014855999592612E-04,
            ),
            radius: 6371.0,
        })
        .with_body(Body {
            name: "Moon".to_owned(),
            gm: 4902.800066,
            pos: Vector3::new(
                -2.659668775178492E+07,
                1.442683153167126E+08,
                6.680827660505474E+04,
            ),
            vel: Vector3::new(
                -2.926974096801152E+01,
                -6.020397935372383E+00,
                -1.740818643718001E-03,
            ),
            radius: 1737.0,
        })
        .with_body(Body {
            name: "Mars".to_owned(),
            gm: 42828.375214,
            pos: Vector3::new(
                2.069270543147017E+08,
                -3.560689745239088E+06,
                -5.147936537447235E+06,
            ),
            vel: Vector3::new(
                1.304308833322233E+00,
                2.628158890420931E+01,
                5.188465740839767E-01,
            ),
            radius: 3390.0,
        })
        .with_body(Body {
            name: "Jupiter".to_owned(),
            gm: 126686531.900,
            pos: Vector3::new(
                5.978410555886381E+08,
                4.387048655696349E+08,
                -1.520164176015472E+07,
            ),
            vel: Vector3::new(
                -7.892632213479861E+00,
                1.115034525890079E+01,
                1.305100448596264E-01,
            ),
            radius: 69911.0,
        })
        .with_body(Body {
            name: "Saturn".to_owned(),
            gm: 37931206.159,
            pos: Vector3::new(
                9.576383364792708E+08,
                9.821475307689621E+08,
                -5.518981181311160E+07,
            ),
            vel: Vector3::new(
                -7.419580382572883E+00,
                6.725982471305630E+00,
                1.775012039800541E-01,
            ),
            radius: 58232.0,
        })
        .with_body(Body {
            name: "Uranus".to_owned(),
            gm: 5793951.322,
            pos: Vector3::new(
                2.157706590772995E+09,
                -2.055242872276605E+09,
                -3.559274281048691E+07,
            ),
            vel: Vector3::new(
                4.646953838324629E+00,
                4.614361336011624E+00,
                -4.301369677250144E-02,
            ),
            radius: 25362.0,
        })
        .with_body(Body {
            name: "Neptune".to_owned(),
            gm: 6835099.97,
            pos: Vector3::new(
                2.513785451779509E+09,
                -3.739265135509532E+09,
                1.907027540535474E+07,
            ),
            vel: Vector3::new(
                4.475107938022004E+00,
                3.062850546988970E+00,
                -1.667293921151841E-01,
            ),
            radius: 24624.0,
        });

    let mut integrator = SuzukiIntegrator::new(STEP);
    let mut time = 0.0;

    let mut currently_visible = false;
    let himawari = Himawari::new();

    while time < 23.0 * YEAR {
        integrator.propagate_in_place(
            &mut sim,
            SimState::position_derivative,
            SimState::momentum_derivative,
            StepSize::UseDefault,
        );
        time += STEP;

        let ang_to_moon = himawari.ang_to_moon(&sim, time);
        let obscured = ang_to_moon < 8.45_f64.to_radians();
        let within_frame = himawari.within_frame(&sim, time);

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
