#![allow(unused)]
mod body;

use std::{
    fmt,
    fs::File,
    io::{Read, Write},
    ops::{Add, Div, Mul, Neg, Sub},
    path::Path,
};

use chrono::{DateTime, Duration, NaiveDateTime, Utc};
use nalgebra::{DVector, Vector3};
use num::Zero;
use numeric_algs::symplectic::{
    integration::{Integrator, StepSize},
    State, StateDerivative,
};
use toml::{map::Map, Value};

use super::STEP;

pub use body::Body;

type Position = Vector3<f64>;
type Velocity = Vector3<f64>;

const DIM: usize = 3;

#[derive(Clone)]
pub struct SimState {
    time: DateTime<Utc>,
    bodies: Vec<Body>,
}

impl Default for SimState {
    fn default() -> Self {
        Self::new()
    }
}

impl SimState {
    pub fn new() -> Self {
        Self {
            time: DateTime::<Utc>::from_utc(NaiveDateTime::from_timestamp(0, 0), Utc),
            bodies: Vec::new(),
        }
    }

    pub fn time(&self) -> DateTime<Utc> {
        self.time
    }

    pub fn with_timestamp(mut self, timestamp: f64) -> Self {
        let timestamp_i64 = (timestamp * 1e9) as i64;
        let timestamp_s = timestamp_i64 / 1_000_000_000;
        let timestamp_ns = (timestamp_i64 % 1_000_000_000) as u32;
        self.time = DateTime::<Utc>::from_utc(
            NaiveDateTime::from_timestamp(timestamp_s, timestamp_ns),
            Utc,
        );
        self
    }

    pub fn with_datestr(mut self, datestr: &str) -> Self {
        let naive_datetime = NaiveDateTime::parse_from_str(datestr, "%Y-%m-%dT%H:%M:%SZ").unwrap();
        self.time = DateTime::<Utc>::from_utc(naive_datetime, Utc);
        self
    }

    pub fn with_body(mut self, body: Body) -> Self {
        self.bodies.push(body);
        self
    }

    pub fn body_index_by_name(&self, name: &str) -> Option<usize> {
        self.bodies
            .iter()
            .enumerate()
            .find(|(_, body)| body.name == name)
            .map(|(idx, _)| idx)
    }

    pub fn body_by_name(&self, name: &str) -> Option<&Body> {
        self.bodies.iter().find(|body| body.name == name)
    }

    pub fn position_derivative(&self) -> SimDerivative {
        let mut derivative = Vec::with_capacity(self.bodies.len() * DIM);
        for body in &self.bodies {
            for i in 0..DIM {
                derivative.push(body.vel[i]);
            }
        }
        SimDerivative(DVector::from_vec(derivative))
    }

    pub fn momentum_derivative(&self) -> SimDerivative {
        let mut derivative = vec![0.0; DIM * self.bodies.len()];
        for (i, body) in self.bodies.iter().enumerate() {
            let mut accel: Vector3<f64> = Zero::zero();
            for (i2, body2) in self.bodies.iter().enumerate() {
                if i2 == i {
                    continue;
                }
                let diff = body2.pos - body.pos;
                let dist = body.distance_from(body2);
                let part_accel = body2.gm / (dist * dist);
                accel += part_accel * diff / dist;
            }
            for j in 0..DIM {
                derivative[i * DIM + j] = accel[j];
            }
        }
        SimDerivative(DVector::from_vec(derivative))
    }

    pub fn bodies(&self) -> impl Iterator<Item = &Body> {
        self.bodies.iter()
    }

    pub fn get_body(&self, idx: usize) -> &Body {
        &self.bodies[idx]
    }

    pub fn save<P: AsRef<Path>>(self, path: P) {
        let mut file = File::create(path).expect("should create file");
        let value = Value::from(self);
        file.write_all(value.to_string().as_bytes());
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Self {
        let mut file = File::open(path).expect("should open file");
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer);
        toml::from_slice::<Value>(&buffer)
            .expect("should parse TOML")
            .into()
    }

    pub fn propagate_to<I: Integrator<Self>>(
        &mut self,
        integrator: &mut I,
        datetime: DateTime<Utc>,
    ) {
        while datetime > self.time {
            let step_size =
                STEP.min((datetime - self.time).num_nanoseconds().unwrap() as f64 / 1e9);
            integrator.propagate_in_place(
                self,
                SimState::position_derivative,
                SimState::momentum_derivative,
                StepSize::Step(step_size),
            );
        }
        while datetime < self.time {
            let step_size =
                STEP.min((self.time - datetime).num_nanoseconds().unwrap() as f64 / 1e9);
            integrator.propagate_in_place(
                self,
                SimState::position_derivative,
                SimState::momentum_derivative,
                StepSize::Step(-step_size),
            );
        }
    }

    pub fn step_forwards<I: Integrator<Self>>(&mut self, integrator: &mut I) {
        integrator.propagate_in_place(
            self,
            SimState::position_derivative,
            SimState::momentum_derivative,
            StepSize::UseDefault,
        );
    }

    pub fn time_since(&self, epoch: DateTime<Utc>) -> Duration {
        self.time - epoch
    }
}

impl State for SimState {
    type PositionDerivative = SimDerivative;
    type MomentumDerivative = SimDerivative;

    fn shift_position_in_place(&mut self, dir: &SimDerivative, amount: f64) {
        for (i, body) in self.bodies.iter_mut().enumerate() {
            for j in 0..DIM {
                body.pos[j] += dir.0[i * DIM + j] * amount;
            }
        }
        self.time = self.time + Duration::nanoseconds((amount * 1e9) as i64);

        // rounding hack
        let time_ns = self.time.timestamp_subsec_nanos() as i64;
        if time_ns > 999_999_000 {
            self.time = self.time + Duration::nanoseconds(1_000_000_000 - time_ns);
        } else if time_ns < 1_000 {
            self.time = self.time - Duration::nanoseconds(time_ns);
        }
    }

    fn shift_momentum_in_place(&mut self, dir: &SimDerivative, amount: f64) {
        for (i, body) in self.bodies.iter_mut().enumerate() {
            for j in 0..DIM {
                body.vel[j] += dir.0[i * DIM + j] * amount;
            }
        }
    }
}

impl fmt::Debug for SimState {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        for (i, body) in self.bodies.iter().enumerate() {
            writeln!(formatter, "{}. {:?}", i + 1, body)?;
        }
        Ok(())
    }
}

impl From<SimState> for Value {
    fn from(sim: SimState) -> Value {
        let mut map = Map::new();
        map.insert(
            "bodies".to_string(),
            Value::Array(sim.bodies.into_iter().map(Value::from).collect()),
        );
        Value::Table(map)
    }
}

impl From<Value> for SimState {
    fn from(val: Value) -> SimState {
        match val {
            Value::Table(mut map) => match map.remove("bodies") {
                Some(Value::Array(bodies)) => SimState {
                    time: DateTime::<Utc>::from_utc(NaiveDateTime::from_timestamp(0, 0), Utc),
                    bodies: bodies.into_iter().map(Body::from).collect(),
                },
                _ => panic!("should have an array"),
            },
            _ => panic!("should be a table"),
        }
    }
}

#[derive(Clone)]
pub struct SimDerivative(DVector<f64>);

impl Add<SimDerivative> for SimDerivative {
    type Output = SimDerivative;

    fn add(self, other: SimDerivative) -> SimDerivative {
        SimDerivative(self.0 + other.0)
    }
}

impl Sub<SimDerivative> for SimDerivative {
    type Output = SimDerivative;

    fn sub(self, other: SimDerivative) -> SimDerivative {
        SimDerivative(self.0 - other.0)
    }
}

impl Mul<f64> for SimDerivative {
    type Output = SimDerivative;

    fn mul(self, other: f64) -> SimDerivative {
        SimDerivative(self.0 * other)
    }
}

impl Div<f64> for SimDerivative {
    type Output = SimDerivative;

    fn div(self, other: f64) -> SimDerivative {
        SimDerivative(self.0 / other)
    }
}

impl Neg for SimDerivative {
    type Output = SimDerivative;

    fn neg(self) -> SimDerivative {
        SimDerivative(-self.0)
    }
}

impl StateDerivative for SimDerivative {}
