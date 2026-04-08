use std::sync::Arc;
use std::thread;
use std::time::Duration;

use crate::output::Output;

use super::output::TuiOutput;

pub fn spawn_demo(out: Arc<TuiOutput>) {
    thread::spawn(move || {
        for i in 0..3 {
            out.status(&format!("demo status {i}"));
            thread::sleep(Duration::from_millis(200));
        }
        let p = out.progress(100, "demo");
        for _ in 0..20 {
            p.inc(5);
            thread::sleep(Duration::from_millis(60));
        }
        p.finish("demo done");
        out.success("demo finished");
    });
}
