//! Process-wide dependency container. One Engine per command — built
//! once at the top of `main.rs::run` (or inside the TUI worker thread)
//! and threaded into every command. Owns Store, DfEngine, Resources,
//! Config, and a memoized AnnotationRegistry.

use std::path::PathBuf;
use std::sync::OnceLock;

use crate::config::Config;
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::resource::Resources;
use crate::store::annotation::AnnotationRegistry;
use crate::store::cohort::{CohortHandle, CohortId};
use crate::store::config::StoreConfig;
use crate::store::Store;

pub struct Engine {
    store: Store,
    df: DfEngine,
    resources: Resources,
    config: Option<Config>,
    annotation_registry: OnceLock<AnnotationRegistry>,
}

impl Engine {
    /// Configured mode — every command except ingest, init, setup, data,
    /// and uninstall. Errors with `DataMissing` if `cohort setup` has not
    /// been run.
    pub fn open(store_path: Option<PathBuf>) -> Result<Self, CohortError> {
        let config = Config::load_configured()?;
        Self::build(store_path, Some(config))
    }

    /// Unconfigured mode — ingest only. Allowed to run before `cohort
    /// setup` because ingest just builds variant sets from raw VCFs/TSVs.
    /// If a config file exists and is configured, the engine still picks
    /// it up so build-detection in `apply_build_override` keeps working.
    pub fn open_unconfigured(store_path: Option<PathBuf>) -> Result<Self, CohortError> {
        let config = Config::load().ok().filter(|c| c.is_configured());
        Self::build(store_path, config)
    }

    fn build(
        store_path: Option<PathBuf>,
        config: Option<Config>,
    ) -> Result<Self, CohortError> {
        let resources = match &config {
            Some(c) => Resources::detect_with_config(&c.resources),
            None => Resources::detect(),
        };
        let store = Store::open(StoreConfig::resolve(store_path)?)?;
        let df = DfEngine::new(&resources)?;

        // Process-global rayon pool. `build_global` returns Err on the
        // second call within one process — that's idempotent by design,
        // since the pool is a one-shot install. Multiple Engine::open
        // calls in one process (TUI dispatching ingest then annotate)
        // safely no-op on calls 2..n.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(resources.threads)
            .build_global();

        Ok(Self {
            store,
            df,
            resources,
            config,
            annotation_registry: OnceLock::new(),
        })
    }

    pub fn store(&self) -> &Store {
        &self.store
    }

    pub fn df(&self) -> &DfEngine {
        &self.df
    }

    pub fn resources(&self) -> &Resources {
        &self.resources
    }

    /// Configured-mode accessor. Panics if called on a bare engine — call
    /// sites that hold a bare engine must use `config_opt` instead. The
    /// type-level split between `open` and `open_unconfigured` makes this
    /// invariant easy to satisfy at the call site.
    pub fn config(&self) -> &Config {
        self.config
            .as_ref()
            .expect("Engine::config called on unconfigured engine")
    }

    pub fn config_opt(&self) -> Option<&Config> {
        self.config.as_ref()
    }

    /// Memoized annotation registry. The snapshot is bound to the engine
    /// lifetime — `cohort data pull` mutates `refs.toml` but does not
    /// hold an engine, so within one command's lifetime the registry is
    /// stable. Errors with `DataMissing` on a bare engine.
    pub fn annotation_registry(&self) -> Result<&AnnotationRegistry, CohortError> {
        if let Some(reg) = self.annotation_registry.get() {
            return Ok(reg);
        }
        let config = self.config_opt().ok_or_else(|| {
            CohortError::DataMissing("Not configured. Run `cohort setup` first.".into())
        })?;
        let reg = self.store.annotations(config)?;
        Ok(self.annotation_registry.get_or_init(|| reg))
    }

    pub fn cohort(&self, id: &CohortId) -> CohortHandle<'_> {
        self.store.cohort(id)
    }
}
