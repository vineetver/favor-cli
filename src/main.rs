use clap::Parser;

mod cli;
pub mod column;
mod commands;
mod config;
mod data;
mod engine;
mod error;
mod ingest;
mod output;
mod resource;
mod runtime;
mod staar;
mod store;
#[cfg(test)]
mod test_fixtures;
mod setup;
mod types;

use cli::{Cli, Command};
use error::CohortError;
use output::OutputMode;

fn main() {
    let cli = Cli::parse();
    // FAVOR_MACHINE forces machine mode regardless of --format / tty.
    // Resolved here at the composition root so the output adapter stays
    // env-free.
    let mode = if std::env::var("FAVOR_MACHINE").is_ok_and(|v| !v.is_empty()) {
        OutputMode::Machine
    } else {
        OutputMode::detect(&cli.format)
    };
    let out = output::create(&mode);

    if let Some(t) = cli.threads {
        std::env::set_var("FAVOR_THREADS", t.to_string());
    }

    let dry_run = cli.dry_run;
    let store_path = cli.store_path.clone();
    let result = run(cli.command, store_path, &*out, &mode, dry_run);

    if let Err(e) = result {
        out.error(&e);
        std::process::exit(e.exit_code());
    }
}

fn run(
    command: Option<Command>,
    store_path: Option<std::path::PathBuf>,
    out: &dyn output::Output,
    mode: &OutputMode,
    dry_run: bool,
) -> Result<(), CohortError> {
    let command = command.ok_or_else(|| {
        CohortError::Input(
            "No command given. Run `favor --help` to see available commands.".into(),
        )
    })?;
    match command {
        Command::Init { path, force } => setup::init(path, force, out, mode),
        Command::Setup {
            root,
            tier,
            packs,
            environment,
            memory_budget,
        } => setup::setup(out, root, tier, packs, environment, memory_budget),
        Command::Data { action } => data::transfer::run(action, out),
        Command::Annotation { action } => match action {
            cli::AnnotationAction::Attach { name, kind, path } => {
                let engine = runtime::Engine::open_unconfigured(store_path)?;
                commands::annotation::attach(&engine, &name, kind, path, out)
            }
        },
        Command::Store { action } => match action {
            cli::StoreAction::Gc => {
                let engine = runtime::Engine::open(store_path)?;
                commands::store::gc(&engine, out)
            }
        },
        Command::Uninstall => setup::uninstall(out),

        Command::Ingest {
            inputs,
            output,
            emit_sql,
            build,
            annotations,
            cohort_id,
            rebuild,
            chromosome,
        } => {
            let engine = runtime::Engine::open_unconfigured(store_path)?;
            let config = commands::ingest::build_config(
                inputs,
                output,
                emit_sql,
                build,
                annotations,
                cohort_id,
                rebuild,
                chromosome,
            )?;
            commands::ingest::run_ingest(&engine, &config, out, dry_run)
        }
        Command::Annotate {
            input,
            output: output_path,
            full,
        } => {
            let engine = runtime::Engine::open(store_path)?;
            let config = commands::annotate::build_config(&engine, input, output_path, full)?;
            commands::annotate::run_annotate(&engine, &config, out, dry_run)
        }
        Command::Enrich {
            input,
            tissue,
            output: output_path,
        } => {
            let engine = runtime::Engine::open(store_path)?;
            let config = commands::enrich::build_config(input, tissue, output_path)?;
            commands::enrich::run_enrich(&engine, &config, out, dry_run)
        }
        Command::Interpret {
            input,
            tissue,
            disease,
            output: output_path,
        } => {
            let engine = runtime::Engine::open(store_path)?;
            commands::interpret::run(&engine, input, tissue, disease, output_path, out)
        }
        Command::Staar {
            genotypes,
            cohort,
            phenotype,
            trait_name,
            covariates,
            annotations,
            masks,
            maf_cutoff,
            window_size,
            individual,
            spa,
            ancestry_col,
            ai_base_tests,
            ai_seed,
            scang_lmin,
            scang_lmax,
            scang_step,
            kinship,
            kinship_groups,
            known_loci,
            null_model,
            emit_sumstats,
            rebuild_store,
            cohort_id,
            no_store: _no_store,
            adaptive: _adaptive,
            column_map,
            output: output_path,
        } => {
            let engine = runtime::Engine::open(store_path)?;
            commands::staar::run(
                &engine,
                commands::staar::StaarArgs {
                    genotypes,
                    cohort,
                    phenotype,
                    trait_names: trait_name,
                    covariates,
                    annotations,
                    masks,
                    maf_cutoff,
                    window_size,
                    individual,
                    spa,
                    ancestry_col,
                    ai_base_tests,
                    ai_seed,
                    scang_lmin,
                    scang_lmax,
                    scang_step,
                    kinship,
                    kinship_groups,
                    known_loci,
                    null_model,
                    emit_sumstats,
                    rebuild_store,
                    column_map,
                    output_path,
                    cohort_id,
                },
                out,
                dry_run,
            )
        }
        Command::Grm {
            cohort,
            king_seg,
            degree,
            n_pcs,
            block_size,
            output: output_path,
        } => {
            let engine = runtime::Engine::open(store_path)?;
            commands::grm::run(
                &engine,
                commands::grm::GrmArgs {
                    cohort,
                    king_seg,
                    degree,
                    n_pcs,
                    block_size,
                    output: output_path,
                },
                out,
                dry_run,
            )
        }
        Command::LdPrune {
            cohort,
            phenotype,
            trait_name,
            covariates,
            variants,
            maf_cutoff,
            cond_p_thresh,
            column_map,
            output: output_path,
        } => {
            let engine = runtime::Engine::open(store_path)?;
            commands::ld_prune::run(
                &engine,
                commands::ld_prune::LdPruneArgs {
                    cohort,
                    phenotype,
                    trait_name,
                    covariates,
                    variants,
                    maf_cutoff,
                    cond_p_thresh,
                    column_map,
                    output: output_path,
                },
                out,
                dry_run,
            )
        }
        Command::MetaStaar {
            studies,
            masks,
            maf_cutoff,
            window_size,
            known_loci,
            conditional_model,
            output: output_path,
        } => {
            let engine = runtime::Engine::open(store_path)?;
            let config = commands::meta_staar::build_config(
                studies,
                masks,
                maf_cutoff,
                window_size,
                known_loci,
                conditional_model,
                output_path,
            )?;
            commands::meta_staar::run_meta_staar(&engine, &config, out, dry_run)
        }
        Command::Schema { table } => {
            let engine = runtime::Engine::open(store_path)?;
            commands::inspect::schema(&engine, table, out)
        }
        Command::Manifest => {
            let engine = runtime::Engine::open(store_path)?;
            commands::inspect::manifest(&engine, out)
        }
    }
}
