use std::fs;

use super::simulation::SimState;

#[derive(Clone, Debug)]
pub struct Snapshots {
    snapshots: Vec<(i64, SimState)>,
}

impl Snapshots {
    pub fn new() -> Self {
        let mut snapshots = vec![];
        for entry in fs::read_dir("snapshots").unwrap() {
            let entry = entry.expect("should be a DirEntry");
            let name = entry.file_name().into_string().expect("should be a string");
            if name.ends_with("state") {
                let timestamp: i64 = name
                    .split(".")
                    .next()
                    .expect("should have timestamp")
                    .parse()
                    .expect("should be a valid i64");
                let state = SimState::load(&format!("snapshots/{}", name));
                snapshots.push((timestamp, state));
            }
        }
        snapshots.sort_by_key(|x| x.0);
        Snapshots { snapshots }
    }

    pub fn get_closest(&self, timestamp: i64) -> (i64, SimState) {
        match self.snapshots.binary_search_by_key(&timestamp, |x| x.0) {
            Ok(i) => self.snapshots[i].clone(),
            Err(i) => {
                if i == 0 {
                    self.snapshots[0].clone()
                } else if i == self.snapshots.len() {
                    self.snapshots[self.snapshots.len() - 1].clone()
                } else {
                    let t1 = self.snapshots[i - 1].0;
                    let t2 = self.snapshots[i].0;
                    if t2 - timestamp > timestamp - t1 {
                        self.snapshots[i - 1].clone()
                    } else {
                        self.snapshots[i].clone()
                    }
                }
            }
        }
    }
}
