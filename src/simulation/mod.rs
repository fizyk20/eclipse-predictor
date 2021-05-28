mod body;

pub use body::Body;
use nalgebra::{DVector, Vector3};
use num::Zero;
use numeric_algs::symplectic::{State, StateDerivative};

use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

type Position = Vector3<f64>;
type Velocity = Vector3<f64>;

const DIM: usize = 3;

#[derive(Clone)]
pub struct SimState {
    bodies: Vec<Body>,
}

impl SimState {
    pub fn new() -> Self {
        Self { bodies: Vec::new() }
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
        let mut derivative = Vec::with_capacity(self.bodies.len() * DIM);
        for _ in 0..DIM * self.bodies.len() {
            derivative.push(0.0);
        }
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
