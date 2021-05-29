mod simulation;

use std::str::FromStr;

use chrono::{DateTime, Datelike, Duration, TimeZone, Utc};
use nalgebra::Vector3;
use numeric_algs::symplectic::{
    integration::{Integrator, StepSize, SuzukiIntegrator},
    State, StateDerivative,
};
use simulation::{Body, SimState};

const STEP: f64 = 300.0;
const YEAR: f64 = 365.25 * 24.0 * 3600.0;

#[derive(Debug, Clone, Copy, PartialEq)]
enum Eclipse {
    PenumbralLunar,
    PartialLunar,
    TotalLunar,
    PartialSolar,
    TotalSolar,
    AnnularSolar,
}

struct EclipseDetector {
    light_dirs: Vec<(f64, Vector3<f64>)>,
}

impl EclipseDetector {
    fn new() -> Self {
        Self {
            light_dirs: Default::default(),
        }
    }

    fn light_dir_for(&self, time: f64) -> Option<Vector3<f64>> {
        match self
            .light_dirs
            .binary_search_by(|entry| entry.0.partial_cmp(&time).unwrap())
        {
            Err(0) => None,
            Err(i) if i == self.light_dirs.len() => None,
            Err(i) => {
                let (t1, vec1) = self.light_dirs[i - 1];
                let (t2, vec2) = self.light_dirs[i];
                Some(vec1 + (vec2 - vec1) / (t2 - t1) * (time - t1))
            }
            Ok(i) => Some(self.light_dirs[i].1),
        }
    }

    fn save_light_dir(&mut self, time: f64, dir: Vector3<f64>) {
        self.light_dirs.push((time, dir));
        if time - self.light_dirs[0].0 > 600.0 + STEP {
            let _ = self.light_dirs.remove(0);
        }
    }

    fn detect_eclipse(&self, sim: &SimState, time: f64) -> Option<Eclipse> {
        let sun = sim.body_by_name("Sun").unwrap();
        let earth = sim.body_by_name("Earth").unwrap();
        let moon = sim.body_by_name("Moon").unwrap();

        // correction for shadow enlargement: https://eclipse.gsfc.nasa.gov/LEcat5/shadow.html
        let re = earth.radius * 1.011;

        let dist = earth.distance_from(sun);
        let delay = dist / 299_792.458;

        let light_dir = self.light_dir_for(time - delay)?.normalize();

        // lunar eclipses

        let shadow_cone_height = re * dist / (sun.radius - re);
        let moon_rel = moon.pos - earth.pos;

        let a = (re / shadow_cone_height).asin();

        let h = moon_rel.dot(&light_dir);
        let r_vec = moon_rel - light_dir * h;
        let r = r_vec.dot(&r_vec).sqrt() / a.cos();
        let h2 = shadow_cone_height - h + r * a.sin();

        if h > 0.0 && (r + moon.radius) / h2 < re / shadow_cone_height {
            return Some(Eclipse::TotalLunar);
        }

        if h > 0.0 && (r - moon.radius) / h2 < re / shadow_cone_height {
            return Some(Eclipse::PartialLunar);
        }

        None
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
    let mut current_eclipse = None;

    let mut eclipse_detector = EclipseDetector::new();

    while time < 23.0 * YEAR {
        integrator.propagate_in_place(
            &mut sim,
            SimState::position_derivative,
            SimState::momentum_derivative,
            StepSize::UseDefault,
        );
        time += STEP;

        let sun = sim.body_by_name("Sun").unwrap();
        let earth = sim.body_by_name("Earth").unwrap();

        eclipse_detector.save_light_dir(time, earth.pos - sun.pos);

        let new_eclipse = eclipse_detector.detect_eclipse(&sim, time);

        if new_eclipse != current_eclipse {
            // find the actual transition moment with 1 sec precision
            let mut time2 = time;
            let mut sim2 = sim.clone();
            let step2 = 1.0;
            loop {
                integrator.propagate_in_place(
                    &mut sim2,
                    SimState::position_derivative,
                    SimState::momentum_derivative,
                    StepSize::Step(-step2),
                );
                time2 -= step2;
                if time - time2 > STEP {
                    panic!("wtf");
                }
                let eclipse = eclipse_detector.detect_eclipse(&sim2, time2);
                if eclipse == current_eclipse {
                    break;
                }
            }

            let date = epoch + Duration::seconds(time2 as i64);
            // correction TT -> UT
            let t = date.year() as f64 + (date.month() as f64 - 0.5) / 12.0 - 2000.0;
            let delta_t = if date.year() < 2005 {
                63.86 + 0.3345 * t - 0.060374 * t * t + 0.0017275 * t * t * t
            } else {
                62.92 + 0.32217 * t + 0.005589 * t * t
            };
            let date = date - Duration::seconds(delta_t as i64);
            //println!("{:#?}", sim.body_by_name("Earth").unwrap());
            //println!("{:#?}", sim.body_by_name("Moon").unwrap());
            //println!("{:#?}", sim.body_by_name("Sun").unwrap());
            if let Some(eclipse) = new_eclipse {
                println!("{:?}: date = {}", eclipse, date);
            } else {
                println!("Eclipse ends: date = {}\n", date);
            }
        }
        current_eclipse = new_eclipse;
    }
}
