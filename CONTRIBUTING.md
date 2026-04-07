# Contributing

## Build

```bash
git clone https://github.com/vineetver/cohort-cli.git
cd cohort-cli
cargo check           # type-check (~10s)
cargo test            # run test suite
cargo build --release # full build (~5 min, pulls DataFusion + Arrow)
```

> **HPC users:** never build on login nodes. Use `srun -p <partition> -c 8 --mem=16G cargo build --release`.

## Workflow

1. Fork the repo and create a branch from `master`
2. Make your changes
3. `cargo check` with zero warnings (`RUSTFLAGS="-D warnings"`)
4. `cargo test` passes
5. Open a pull request

## Code standards

**Every line ships.** No commented-out code, placeholder TODOs, or dead code.

| Principle | In practice |
|-----------|-------------|
| Strong types | Enums over strings. Parse at the boundary, convert immediately. |
| Actionable errors | Every error tells the user what command fixes it. |
| Output contract | All output through the `Output` trait. Never bare `println!`. |
| Memory-aware | Derive from `Resources::detect()`. Batch when it doesn't fit. Never crash. |
| Flat control flow | Early returns, guard clauses, match arms that return directly. |
| Thin commands | Commands validate and dispatch. Logic lives in modules. |
| Statistical invariance | Refactoring must never change p-values. |

## Architecture

```
src/
  main.rs              CLI entry point
  cli.rs               clap definitions
  column.rs            Col enum, single source of truth for all column names
  engine.rs            DataFusion query engine
  types.rs             Domain types: Chromosome, Consequence, AnnotatedVariant
  error.rs             CohortError with exit codes and recovery hints
  output.rs            Output trait: human (colored) + machine (JSON)

  commands/            One file per CLI command (thin dispatch)
  ingest/              Format detection, VCF parsing, SQL generation
  data/                Parquet data layer, MinIO transfer, annotation DB
  staar/               STAAR association pipeline, scoring, sparse genotype store
  setup/               Interactive configuration wizard
```

## Adding a command

1. Create `src/commands/<name>.rs` with a typed config struct
2. Add a variant to `Command` in `src/cli.rs`
3. Route in `src/main.rs`
4. Support `--dry-run` (validate + emit JSON plan) and `--format json`
5. Use `Resources::detect_with_config()` for memory-aware sizing

## Adding a data pack

One entry in `PACKS` in `src/data/mod.rs`. Everything else picks it up.

## Bugs

Use the [bug report template](https://github.com/vineetver/cohort-cli/issues/new?template=bug_report.yml). Include `cohort --version`, platform, and memory allocation.

## License

By contributing you agree your contributions are licensed under [GPL-3.0](LICENSE).
