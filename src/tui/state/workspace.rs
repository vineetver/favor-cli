use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, Sender, TryRecvError};
use std::thread;
use std::time::{Duration, Instant};

use walkdir::WalkDir;

use super::artifacts::{classify, group_order, keep_entry, Artifact};

const BATCH_SIZE: usize = 32;
const FLUSH_INTERVAL: Duration = Duration::from_millis(80);

enum ScanMsg {
    Batch(Vec<Artifact>),
    Done,
}

pub struct WorkspaceState {
    pub cwd: PathBuf,
    pub extra_roots: Vec<PathBuf>,
    pub artifacts: Vec<Artifact>,
    pub focus: usize,
    pub scanning: bool,
    pub scan_started: Instant,
    rx: Option<Receiver<ScanMsg>>,
}

impl WorkspaceState {
    pub fn new(cwd: PathBuf, extra_roots: Vec<PathBuf>) -> Self {
        let mut s = Self {
            cwd,
            extra_roots,
            artifacts: Vec::new(),
            focus: 0,
            scanning: false,
            scan_started: Instant::now(),
            rx: None,
        };
        s.start_scan();
        s
    }

    pub fn rescan(&mut self) {
        self.artifacts.clear();
        self.focus = 0;
        self.start_scan();
    }

    fn start_scan(&mut self) {
        let (tx, rx) = mpsc::channel();
        self.rx = Some(rx);
        self.scanning = true;
        self.scan_started = Instant::now();
        let cwd = self.cwd.clone();
        let extra = self.extra_roots.clone();
        thread::spawn(move || background_scan(cwd, extra, tx));
    }

    /// Pull any batches that have arrived. Returns true if `artifacts` changed.
    /// Sorts the list once when the scan completes; while in flight the list
    /// stays in arrival order so streaming feels live.
    pub fn drain_scan(&mut self) -> bool {
        let Some(rx) = &self.rx else {
            return false;
        };
        let mut changed = false;
        let mut finished = false;
        loop {
            match rx.try_recv() {
                Ok(ScanMsg::Batch(items)) => {
                    self.artifacts.extend(items);
                    changed = true;
                }
                Ok(ScanMsg::Done) => {
                    finished = true;
                    changed = true;
                    break;
                }
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => {
                    finished = true;
                    break;
                }
            }
        }
        if finished {
            self.rx = None;
            self.scanning = false;
            self.artifacts.sort_by(|a, b| {
                group_order(&a.kind)
                    .cmp(&group_order(&b.kind))
                    .then_with(|| a.display_name.cmp(&b.display_name))
            });
        }
        changed
    }

    pub fn focused(&self) -> Option<&Artifact> {
        self.artifacts.get(self.focus)
    }

    pub fn move_focus(&mut self, delta: isize) {
        if self.artifacts.is_empty() {
            self.focus = 0;
            return;
        }
        let len = self.artifacts.len() as isize;
        let next = (self.focus as isize + delta).rem_euclid(len);
        self.focus = next as usize;
    }
}

fn background_scan(cwd: PathBuf, extra: Vec<PathBuf>, tx: Sender<ScanMsg>) {
    let mut seen: HashSet<PathBuf> = HashSet::new();
    let mut batch: Vec<Artifact> = Vec::with_capacity(BATCH_SIZE);
    let mut last_flush = Instant::now();

    walk(&cwd, 3, &mut seen, &mut batch, &mut last_flush, &tx);
    for root in &extra {
        walk(root, 1, &mut seen, &mut batch, &mut last_flush, &tx);
    }

    if !batch.is_empty() {
        let _ = tx.send(ScanMsg::Batch(batch));
    }
    let _ = tx.send(ScanMsg::Done);
}

fn walk(
    root: &Path,
    depth: usize,
    seen: &mut HashSet<PathBuf>,
    batch: &mut Vec<Artifact>,
    last_flush: &mut Instant,
    tx: &Sender<ScanMsg>,
) {
    if !root.exists() {
        return;
    }
    for entry in WalkDir::new(root)
        .max_depth(depth)
        .follow_links(false)
        .into_iter()
        .filter_entry(|e| keep_entry(e.path(), e.depth()))
    {
        let Ok(entry) = entry else { continue };
        let path = entry.path().to_path_buf();
        if !seen.insert(path.clone()) {
            continue;
        }
        if let Some(art) = classify(&path) {
            batch.push(art);
        }
        if batch.len() >= BATCH_SIZE || last_flush.elapsed() >= FLUSH_INTERVAL {
            let chunk = std::mem::take(batch);
            if !chunk.is_empty() && tx.send(ScanMsg::Batch(chunk)).is_err() {
                return;
            }
            *last_flush = Instant::now();
        }
    }
}
