use std::fmt;

use nalgebra::Vector3;
use toml::{map::Map, Value};

use super::{Position, Velocity};

#[derive(Clone)]
pub struct Body {
    pub name: String,

    pub gm: f64,
    pub pos: Position,
    pub vel: Velocity,

    pub radius: f64,
}

impl Body {
    pub fn distance_from(&self, other: &Body) -> f64 {
        let diff = self.pos - other.pos;
        diff.dot(&diff).sqrt()
    }
}

impl fmt::Debug for Body {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        write!(
            formatter,
            "{} {{ GM = {}; pos = [{:10.3} ; {:10.3} ; {:10.3}]; \
            vel = [{:10.3} ; {:10.3} ; {:10.3} ]}}",
            self.name,
            self.gm,
            self.pos[0],
            self.pos[1],
            self.pos[2],
            self.vel[0],
            self.vel[1],
            self.vel[2]
        )
    }
}

impl From<Body> for Value {
    fn from(body: Body) -> Value {
        let mut map = Map::new();
        map.insert("name".to_string(), Value::String(body.name));
        map.insert("gm".to_string(), Value::Float(body.gm));
        map.insert("radius".to_string(), Value::Float(body.radius));
        map.insert(
            "position".to_string(),
            Value::Array(vec![
                Value::Float(body.pos.x),
                Value::Float(body.pos.y),
                Value::Float(body.pos.z),
            ]),
        );
        map.insert(
            "velocity".to_string(),
            Value::Array(vec![
                Value::Float(body.vel.x),
                Value::Float(body.vel.y),
                Value::Float(body.vel.z),
            ]),
        );
        Value::Table(map)
    }
}

impl From<Value> for Body {
    fn from(val: Value) -> Body {
        let mut map = match val {
            Value::Table(map) => map,
            val => panic!("should be table, got {:?}", val),
        };
        let name = match map.remove("name").expect("should have a name") {
            Value::String(name) => name,
            _ => panic!("should be a string"),
        };
        let gm = match map.remove("gm").expect("should have a gm") {
            Value::Float(gm) => gm,
            _ => panic!("should be a float"),
        };
        let radius = match map.remove("radius").expect("should have a radius") {
            Value::Float(radius) => radius,
            _ => panic!("should be a float"),
        };
        let pos = match map.remove("position").expect("should have a position") {
            Value::Array(pos) => {
                assert!(pos.len() == 3);
                Vector3::new(
                    pos[0].as_float().expect("should be a float"),
                    pos[1].as_float().expect("should be a float"),
                    pos[2].as_float().expect("should be a float"),
                )
            }
            _ => panic!("should be an array"),
        };
        let vel = match map.remove("velocity").expect("should have a velocity") {
            Value::Array(vel) => {
                assert!(vel.len() == 3);
                Vector3::new(
                    vel[0].as_float().expect("should be a float"),
                    vel[1].as_float().expect("should be a float"),
                    vel[2].as_float().expect("should be a float"),
                )
            }
            _ => panic!("should be an array"),
        };
        Body {
            name,
            gm,
            radius,
            pos,
            vel,
        }
    }
}
