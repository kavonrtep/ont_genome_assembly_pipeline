# Repository Guidelines

## Project Structure & Module Organization

This repository contains a Snakemake-based Oxford Nanopore genome assembly and analysis pipeline. The main workflow is `Snakefile`; read mapping against an existing assembly is `Snakefile_mapping`. Entry points are `run_pipeline.py` and `run_mapping.py`, designed to execute inside the Singularity container. Custom Python, R, and shell helpers live in `scripts/`. Conda environments are in `envs/`. Example inputs and small references are in `data/`; generated files belong under the configured `output_dir`, `logs/`, or temporary directories.

## Build, Test, and Development Commands

- `singularity build images/genome_assembly_pipeline.sif Singularity`: build the container image locally.
- `singularity run -B /path/to/data -B $PWD images/genome_assembly_pipeline.sif -c config.yaml -t 20`: run the full assembly workflow with 20 threads.
- `snakemake --snakefile Snakefile --configfile config_testing.yaml --cores 2 --use-conda --dry-run`: validate the main workflow DAG without running jobs.
- `snakemake --snakefile Snakefile_mapping --configfile config_mapping_only.yaml --cores 2 --use-conda --dry-run`: validate the mapping-only workflow.

Use `config_template.yaml` as the starting point for new runs.

## Coding Style & Naming Conventions

Python scripts use 4-space indentation, `snake_case` names, and explicit argument parsing with `argparse`. Keep command-line scripts executable and include a `#!/usr/bin/env python3` shebang when run directly. Snakemake rule names should be descriptive action names, such as `detect_contamination_mito` or `map_hifi_reads_to_assembly`. Prefer paths built with `os.path.join` in workflow Python blocks. R scripts should keep input/output arguments explicit.

## Testing Guidelines

There is no formal unit test suite or coverage threshold. Before submitting changes, run a Snakemake dry-run for any affected workflow and, when practical, execute the smallest sample-data run using `config_testing.yaml`. For script changes, run the script with `--help` or a minimal input to confirm argument parsing and expected outputs. Check `logs/` when rules fail.

## Commit & Pull Request Guidelines

Recent commits use short, plain-language summaries such as `read mapping added` and `bugfixes`. Keep commit messages concise and focused on one change. Pull requests should describe the workflow or script affected, list the validation commands run, and note any required config changes, external data paths, or container rebuilds. Include screenshots only for visual outputs such as dotplots or reports.

## Security & Configuration Tips

Do not commit large generated assemblies, BAM files, private datasets, or machine-specific absolute paths. Verify Singularity bind mounts before running workflows so all configured input and output paths are visible inside the container.
