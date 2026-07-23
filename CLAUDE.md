# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Snakemake-based Oxford Nanopore / HiFi genome assembly and analysis pipeline, shipped as a Singularity container. It exists in three parallel workflow variants (see Architecture) selected at runtime via `run_pipeline.py -w`.

## Commands

Validate a workflow's DAG without running anything (fastest way to check a Snakefile edit is syntactically and logically sound):

```bash
snakemake --snakefile Snakefile --configfile config_testing.yaml --cores 2 --use-conda --dry-run
snakemake --snakefile Snakefile_generic --configfile tests/fixtures/config_generic_ci.yaml --cores 2 --dry-run
snakemake --snakefile Snakefile_mapping --configfile config_mapping_only.yaml --cores 2 --use-conda --dry-run
```

Check Python syntax and version consistency (what CI runs first):

```bash
python -m py_compile run_pipeline.py run_mapping.py version.py
python version.py
```

Run the full pipeline through the entry point (mirrors how it runs inside the container):

```bash
python run_pipeline.py -w generic -c config_generic_template.yaml -t 4 -S "--dry-run"
python run_pipeline.py --print-config-template generic   # inspect a template
python run_pipeline.py --version
```

Build and smoke-test the Singularity image locally:

```bash
singularity build images/genome_assembly_pipeline.sif Singularity
singularity run -B /path/to/data -B $PWD images/genome_assembly_pipeline.sif -c config.yaml -t 20
```

Fetch the CI fixture and run the generic workflow against it:

```bash
tests/fixtures/fetch_bioinformatics_hifi_fixture.sh
singularity run -B $PWD images/genome_assembly_pipeline_v12.sif \
  -w generic -c tests/fixtures/config_bioinformatics_hifi.yaml -t 2
```

There is no unit test suite. Validation is: a Snakemake `--dry-run` for any affected workflow, plus (when practical) a real run against `config_testing.yaml` or the fixtures in `tests/fixtures/`. There is a `pipeline` GitHub Actions workflow that dry-runs `Snakefile_generic` against the CI fixtures on every PR/push to `main`, and a `release` workflow that builds and publishes the `.sif` on tag push (tag must equal `version.py`'s `__version__` exactly).

## Architecture

Three independent Snakemake workflows share one `scripts/` directory and one Singularity container, but are otherwise separate pipelines with separate config schemas:

- **`Snakefile`** (workflow name `pisum`/`legacy`) — the original, Pisum-specific full pipeline: contamination filtering → hifiasm assembly → Quast → painting-probe BLAST + dotplots → ONT read-back mapping → clipping/zero-coverage analysis → contig-pair extraction → TideCluster → IGV session. Flat config keys (`output_dir`, `oxford_nanopore_reads`, `painting_probes`, `assembly_fai_ref`, ...).
- **`Snakefile_generic`** — a species-agnostic assembly + per-haplotype annotation workflow. Contamination filtering (optional) → hifiasm → GFA→FASTA per haplotype → then every downstream step **fans out over the haplotype set**. Two orthogonal axes drive it: (1) `assembly.mode` (`primary | phased | hic | trio | triploid_hic`) determines the hifiasm invocation *and* the candidate haplotype set (`MODE_HAP_GFA` maps each label to its hifiasm GFA suffix). A `discover_haplotypes` **checkpoint** after assembly prunes that candidate set to the haplotypes hifiasm actually emitted (a homozygous genome or a purge-off run may yield only the primary — `bp.p_ctg` + empty `bp.a_ctg`, no `bp.hap1/hap2`), and every downstream step fans out over that **realized** set via `realized_haplotypes()` rather than the static candidates. (2) an optional `parental_assignment` block labels the assembled contigs by parental origin/subgenome via `yak triobin` (Hi-C phases but cannot assign parentage — parental Illumina is required for that). Per-haplotype steps, each toggled by its own `enabled`: `quast`, `map_reads_back`, `tidecluster` (TR + optional custom library, TideCluster 1.17), `telomere` (tidk, configurable motif), `reference_assignment` (RagTag scaffold vs a single external `reference.fasta`), `dotplot_vs_reference` (minimap2 → PAF → `scripts/paf_dotplot.py`). Nested config schema (`reads.*`, `assembly.*`, `parental_assignment.*`, `reference.fasta`, per-step blocks, `outputs.output_dir`); legacy `hifiasm.*` / `steps.*` keys are still read as fallbacks (old configs default to `mode: primary`). Trio uses `-1/-2` parental yak (needs both parents); Hi-C uses `--h1/--h2` (+`--enzyme`, `--n-hap 3` for `triploid_hic`) — the pinned hifiasm fork accepts both flag sets but does not jointly phase (trio takes precedence).
- **`Snakefile_mapping`** — mapping-only workflow: takes an existing assembly plus ONT and HiFi reads, maps both back, and runs the same clipping/zero-coverage/contig-pair/TideCluster/IGV analysis as the full pipeline. Flat config keys (`genome_assembly`, `oxford_nanopore_reads`, `hifi_reads`, `tandem_repeats`, ...).

`run_pipeline.py` is the single entry point (`%runscript` in `Singularity`). Its `WORKFLOWS` dict is the authoritative map of workflow name → snakefile, config template, output-path config keys, and required-input config keys. It validates required inputs exist, creates the output dir, prints Singularity bind-mount hints on failure, and shells out to `snakemake --use-conda`. `run_mapping.py` is a thin wrapper specific to the `mapping` workflow. **When adding a new workflow variant or renaming config keys, update `WORKFLOWS` in `run_pipeline.py` in the same change** — the required-file and output-path lookups there are the only place that logic lives.

Each rule declares its own `conda:` env (see `envs/*.yaml`: `pysam`, `basic_tools`, `quast`, `seqkit`, `tidecluster`, plus `tidk`, `ragtag`, `yak`, `dotplot` for the generic annotation steps) rather than one shared environment — keep new rules consistent with this per-rule-env pattern. Rules also do `scripts_dir=$(realpath scripts); export PATH=$scripts_dir:$PATH` to call helpers in `scripts/` by name; hifiasm comes from the container's system PATH, not conda. All ONT read-to-assembly mapping uses `minimap2 -x lr:hq` (from `envs/basic_tools.yaml`) — there is no `dorado` dependency.

`version.py` is the single source of truth for the pipeline version (`__version__`, PEP-440-lite). The Singularity build stamps it into `/.singularity.d/labels.json`; the release workflow requires the git tag to equal it exactly.

`config_*_template.yaml` files are the canonical, checked-in defaults for each workflow (`config_template.yaml`/`config_generic_template.yaml`/`config_mapping_template.yaml`); user-specific configs like `config.yaml`, `config_testing.yaml`, `config_mapping_only.yaml` are working copies and may contain machine-specific paths — don't assume they're portable.

`scripts/` holds standalone Python/R/shell helpers invoked from Snakemake `shell:` blocks (contamination detection, BAM filtering/clipping analysis, dotplot generation, contig-pair extraction). They're independent CLI tools with `argparse`, not a shared library — check a script's `--help` before changing its interface, since Snakefile rules call it by fixed flag names.

## Conventions

- Snakemake rule names are descriptive verb phrases (`detect_contamination_mito`, `map_nanopore_reads_to_assembly`), not generic names.
- Python scripts: 4-space indent, `snake_case`, explicit `argparse`, `#!/usr/bin/env python3` shebang, executable.
- Build workflow paths with `os.path.join` in Snakefile Python blocks.
- Don't commit large generated assemblies, BAM files, private datasets, machine-specific absolute paths, or `.sif` images.
