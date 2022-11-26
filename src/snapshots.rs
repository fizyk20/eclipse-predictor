use std::fs;

use chrono::{DateTime, Utc};

use super::simulation::SimState;

#[derive(Clone, Debug)]
pub struct Snapshots {
    snapshots: Vec<SimState>,
}

impl Snapshots {
    pub fn new() -> Self {
        let mut snapshots = vec![];
        for entry in fs::read_dir("snapshots").unwrap() {
            let entry = entry.expect("should be a DirEntry");
            let name = entry.file_name().into_string().expect("should be a string");
            if name.ends_with("state") {
                let datestr: &str = name.split('.').next().expect("should have timestamp");
                let state = SimState::load(&format!("snapshots/{}", name)).with_datestr(datestr);
                snapshots.push(state);
            }
        }
        snapshots.sort_by_key(|x| x.time());
        Snapshots { snapshots }
    }

    pub fn get_closest(&self, datetime: DateTime<Utc>) -> SimState {
        match self.snapshots.binary_search_by_key(&datetime, |x| x.time()) {
            Ok(i) => self.snapshots[i].clone(),
            Err(i) => {
                if i == 0 {
                    self.snapshots[0].clone()
                } else if i == self.snapshots.len() {
                    self.snapshots[self.snapshots.len() - 1].clone()
                } else {
                    let t1 = self.snapshots[i - 1].time();
                    let t2 = self.snapshots[i].time();
                    if t2 - datetime > datetime - t1 {
                        self.snapshots[i - 1].clone()
                    } else {
                        self.snapshots[i].clone()
                    }
                }
            }
        }
    }

    pub fn insert(&mut self, sim: SimState) {
        match self
            .snapshots
            .binary_search_by_key(&sim.time(), |x| x.time())
        {
            Ok(_) => (), // exists already - nothing to do
            Err(i) => {
                let name = format!(
                    "snapshots/{}.state",
                    sim.time().format("%Y-%m-%dT%H:%M:%SZ")
                );
                sim.clone().save(&name);
                self.snapshots.insert(i, sim);
            }
        }
    }
}
