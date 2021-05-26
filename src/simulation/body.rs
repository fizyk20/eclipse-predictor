use super::{Position, Velocity};
use std::fmt;

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
