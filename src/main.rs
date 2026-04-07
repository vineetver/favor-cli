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
mod staar;
#[cfg(test)]
mod test_fixtures;
mod setup;
mod types;

use cli::{Cli, Command};
use error::FavorError;
use output::OutputMode;

fn main() {
    let cli = Cli::parse();
    let mode = OutputMode::detect(&cli.format);
    let out = output::create(&mode);

    if let Some(t) = cli.threads {
        std::env::set_var("FAVOR_THREADS", t.to_string());
    }

    let dry_run = cli.dry_run;
    let result = run(cli.command, &*out, &mode, dry_run);

    if let Err(e) = result {
        out.error(&e);
        std::process::exit(e.exit_code());
    }
}

fn run(
    command: Command,
    out: &dyn output::Output,
    mode: &OutputMode,
    dry_run: bool,
) -> Result<(), FavorError> {
    match command {
        Command::Init { path, force } => setup::init(path, force, out, mode),
        Command::Setup { environment, memory_budget } => {
            setup::setup(out, mode, environment, memory_budget)
        }
        Command::Data { action } => data::transfer::run(action, out),
        Command::Ingest { inputs, output, emit_sql, build } => {
            commands::ingest::run(inputs, output, emit_sql, build, out, dry_run)
        }
        Command::Annotate {
            input,
            output: output_path,
            full,
        } => commands::annotate::run(input, output_path, full, out, dry_run),
        Command::Enrich {
            input,
            tissue,
            output: output_path,
        } => commands::enrich::run(input, tissue, output_path, out, dry_run),
        Command::Interpret {
            input,
            tissue,
            disease,
            output: output_path,
        } => commands::interpret::run(input, tissue, disease, output_path, out),
        Command::Staar {
            genotypes,
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
            emit_sumstats,
            rebuild_store,
            no_store: _no_store,
            adaptive: _adaptive,
            column_map,
            output: output_path,
        } => commands::staar::run(
            commands::staar::StaarArgs {
                genotypes,
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
                emit_sumstats,
                rebuild_store,
                column_map,
                output_path,
            },
            out,
            dry_run,
        ),
        Command::MetaStaar {
            studies,
            masks,
            maf_cutoff,
            window_size,
            output: output_path,
        } => commands::meta_staar::run(
            studies, masks, maf_cutoff, window_size, output_path, out, dry_run,
        ),
        Command::Schema { table } => commands::inspect::schema(table, out),
        Command::Manifest => commands::inspect::manifest(out),
        Command::Uninstall => setup::uninstall(out),
    }
}
