mod simulation;

use std::{f64::consts::PI, str::FromStr};

use chrono::{DateTime, Duration, Utc};
use nalgebra::Vector3;
use numeric_algs::symplectic::integration::{Integrator, StepSize, SuzukiIntegrator};
use simulation::{Body, SimState};

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

fn main() {
    let epoch = DateTime::<Utc>::from_str("2022-01-01T00:00:00Z").unwrap();

    let mut sim = SimState::new()
        .with_body(Body {
            name: "Sun".to_owned(),
            gm: 132712440041.93938,
            pos: Vector3::new(
                -1.283674643550172E+06,
                5.007104996950605E+05,
                2.589397504295033E+04,
            ),
            vel: Vector3::new(
                -5.809369653802155E-03,
                -1.461959576560110E-02,
                2.513455442031695E-04,
            ),
            radius: 696000.0,
        })
        .with_body(Body {
            name: "Mercury".to_owned(),
            gm: 22031.86855,
            pos: Vector3::new(
                5.242617205495467E+07,
                -5.596063357617276E+06,
                -5.398976570474024E+06,
            ),
            vel: Vector3::new(
                -3.931719860392732E+00,
                5.056613955108243E+01,
                4.493726800433638E+00,
            ),
            radius: 2440.0,
        })
        .with_body(Body {
            name: "Venus".to_owned(),
            gm: 324858.592,
            pos: Vector3::new(
                -1.143612889654620E+07,
                1.076180391552140E+08,
                2.081921801192194E+06,
            ),
            vel: Vector3::new(
                -3.498958532524220E+01,
                -3.509011592387367E+00,
                1.971012081662609E+00,
            ),
            radius: 6052.0,
        })
        .with_body(Body {
            name: "Earth".to_owned(),
            gm: 398600.435436,
            pos: Vector3::new(
                -2.741147560901964E+07,
                1.452697499646169E+08,
                1.907499306293577E+04,
            ),
            vel: Vector3::new(
                -2.981801522121922E+01,
                -5.415519940416356E+00,
                1.781036907294364E-03,
            ),
            radius: 6371.0,
        })
        .with_body(Body {
            name: "Moon".to_owned(),
            gm: 4902.800066,
            pos: Vector3::new(
                -2.750334415265606E+07,
                1.449229071220244E+08,
                1.107672420730442E+04,
            ),
            vel: Vector3::new(
                -2.875697243341937E+01,
                -5.673128746911559E+00,
                -9.368517608037941E-02,
            ),
            radius: 1737.0,
        })
        .with_body(Body {
            name: "Mars".to_owned(),
            gm: 42828.375214,
            pos: Vector3::new(
                -1.309510737126251E+08,
                -1.893127398896606E+08,
                -7.714450109843910E+05,
            ),
            vel: Vector3::new(
                2.090994471204196E+01,
                -1.160503586188451E+01,
                -7.557181497936503E-01,
            ),
            radius: 3390.0,
        })
        .with_body(Body {
            name: "Jupiter".to_owned(),
            gm: 126686531.900,
            pos: Vector3::new(
                6.955554713494443E+08,
                -2.679620040967891E+08,
                -1.444959769995748E+07,
            ),
            vel: Vector3::new(
                4.539612624165795E+00,
                1.280513202430234E+01,
                -1.547160200183022E-01,
            ),
            radius: 69911.0,
        })
        .with_body(Body {
            name: "Saturn".to_owned(),
            gm: 37931206.159,
            pos: Vector3::new(
                1.039929082221698E+09,
                -1.056650148100382E+09,
                -2.303098768547428E+07,
            ),
            vel: Vector3::new(
                6.345150014839902E+00,
                6.756117343710409E+00,
                -3.704447196861729E-01,
            ),
            radius: 58232.0,
        })
        .with_body(Body {
            name: "Uranus".to_owned(),
            gm: 5793951.322,
            pos: Vector3::new(
                2.152570437700128E+09,
                2.016888245555490E+09,
                -2.039611192913723E+07,
            ),
            vel: Vector3::new(
                -4.705853565766252E+00,
                4.652144641704226E+00,
                7.821724397220797E-02,
            ),
            radius: 25362.0,
        })
        .with_body(Body {
            name: "Neptune".to_owned(),
            gm: 6835099.97,
            pos: Vector3::new(
                4.431790029647159E+09,
                -6.114486878261740E+08,
                -8.954348455674592E+07,
            ),
            vel: Vector3::new(
                7.066248666117450E-01,
                5.417073900824804E+00,
                -1.271380125837618E-01,
            ),
            radius: 24624.0,
        });

    let mut integrator = SuzukiIntegrator::new(STEP);
    let mut time = 0.0;

    let mut currently_visible = false;
    let himawari = Himawari::new();

    while time < 2.0 * YEAR {
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
