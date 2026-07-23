#!/usr/bin/env python3
"""
This script reads a user-provided config file and launches the Snakemake pipeline.
It is designed to run within your Singularity container.
"""

import argparse
import os
import shlex
import subprocess
import sys
import yaml

from version import __version__


PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))

WORKFLOWS = {
    "pisum": {
        "snakefile": "Snakefile",
        "template": "config_template.yaml",
        "fallback_template": "config.yaml",
        "output_path": ["output_dir"],
        "required_files": [
            "oxford_nanopore_reads", "painting_probes", "assembly_fai_ref",
            "blast_probes_ref", "tandem_repeats"
        ],
    },
    "legacy": {
        "snakefile": "Snakefile",
        "template": "config_template.yaml",
        "fallback_template": "config.yaml",
        "output_path": ["output_dir"],
        "required_files": [
            "oxford_nanopore_reads", "painting_probes", "assembly_fai_ref",
            "blast_probes_ref", "tandem_repeats"
        ],
    },
    "generic": {
        "snakefile": "Snakefile_generic",
        "template": "config_generic_template.yaml",
        "fallback_template": None,
        "output_path": ["outputs", "output_dir"],
        "required_files": [
            "reads.ont_fastq", "contamination_filter.plastid_db",
            "contamination_filter.mitochondrial_db"
        ],
    },
    "mapping": {
        "snakefile": "Snakefile_mapping",
        "template": "config_mapping_template.yaml",
        "fallback_template": None,
        "output_path": ["output_dir"],
        "required_files": [
            "genome_assembly", "oxford_nanopore_reads", "hifi_reads",
            "tandem_repeats"
        ],
    },
}


def pipeline_path(filename):
    return os.path.join(PIPELINE_DIR, filename)


def config_get(config_object, dotted_key, default=None):
    value = config_object
    for key in dotted_key.split("."):
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def config_get_path(config_object, path, default=None):
    value = config_object
    for key in path:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def resolve_required_files(workflow_name, config_object):
    """Return the list of config keys that must point to existing files.

    For the generic workflow this depends on the chosen assembly.mode and which
    per-haplotype steps are enabled, mirroring the validation in
    Snakefile_generic. Other workflows use their static required_files list.
    """
    workflow = WORKFLOWS[workflow_name]
    if workflow_name != "generic":
        return list(workflow["required_files"])

    def gget(key, default=None):
        return config_get(config_object, key, default)

    # Long-read input. reads.hifi_fastq is preferred for HiFi, but older configs
    # (and the CI fixture) put HiFi reads under reads.ont_fastq with
    # read_type: hifi, so fall back to ont_fastq when hifi_fastq is unset.
    read_type = gget("assembly.read_type", gget("hifiasm.read_type", "ont"))
    if read_type == "hifi" and gget("reads.hifi_fastq"):
        reads_key = "reads.hifi_fastq"
    else:
        reads_key = "reads.ont_fastq"
    required = [reads_key]

    if gget("contamination_filter.enabled", True) is not False:
        required += ["contamination_filter.plastid_db",
                     "contamination_filter.mitochondrial_db"]

    mode = gget("assembly.mode", "primary")
    if mode in ("hic", "triploid_hic"):
        required += ["reads.hic_r1", "reads.hic_r2"]
    if mode == "trio":
        required += ["assembly.trio.paternal_illumina",
                     "assembly.trio.maternal_illumina"]

    if gget("reference_assignment.enabled", False) or \
            gget("dotplot_vs_reference.enabled", False):
        required += ["reference.fasta"]

    if gget("parental_assignment.enabled", False):
        required += ["parental_assignment.paternal_illumina",
                     "parental_assignment.maternal_illumina"]

    return required


def read_template(workflow_name):
    workflow = WORKFLOWS[workflow_name]
    template_candidates = [workflow["template"], workflow.get("fallback_template")]
    for template_name in template_candidates:
        if not template_name:
            continue
        template_path = pipeline_path(template_name)
        if os.path.exists(template_path):
            with open(template_path, 'r') as f:
                return f.read()
    return f"Could not read {workflow['template']}."


def show_singularity_settings(config_object, workflow_name):
    """
    Print the Singularity bind options based on directories from required files.
    """
    workflow = WORKFLOWS[workflow_name]
    dirs = set()
    # Always include the output directory.
    output_dir = config_get_path(config_object, workflow["output_path"], "output")
    dirs.add(os.path.abspath(output_dir))
    for key in resolve_required_files(workflow_name, config_object):
        value = config_get(config_object, key)
        for path in str(value).split(",") if value else []:
            path = path.strip()
            if path:
                dirs.add(os.path.dirname(os.path.abspath(path)))
    bind_string = " ".join([f"-B {d}" for d in dirs])
    print("Run Singularity with the following bind options:")
    print(f"singularity run {bind_string} ...")


def main():
    parser = argparse.ArgumentParser(
        description="Pipeline for assembly and analysis using Snakemake within a "
                    "Singularity container.",
        epilog="Use --print-config-template to inspect a workflow-specific template.",
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument(
        '-w', '--workflow', choices=sorted(WORKFLOWS), default='pisum',
        help='Workflow to run. Default keeps the original Pisum/full pipeline.'
        )
    parser.add_argument('-c', '--config', required=False, help='Path to config file')
    parser.add_argument(
        '-t', '--threads', required=False, default=2, type=int,
        help='Number of threads to use'
        )
    parser.add_argument(
        '-S', '--snakemake_args', nargs=argparse.REMAINDER, default=[],
        help='Additional snakemake arguments (e.g., "--dry-run --reason"). Do not '
             'include options that are set by the script. Put this option last.'
        )
    parser.add_argument(
        '--conda-frontend', choices=['mamba', 'conda'], default='mamba',
        help='Conda frontend passed to Snakemake. The container default is mamba.'
        )
    parser.add_argument(
        '--print-config-template', choices=sorted(WORKFLOWS), metavar='WORKFLOW',
        help='Print a config template for the selected workflow and exit.'
        )
    parser.add_argument(
        '--version', action='store_true',
        help='Print the pipeline version and exit.'
        )
    args = parser.parse_args()

    if args.version:
        print(f"genome_assembly_pipeline {__version__}")
        return

    if args.print_config_template:
        print(read_template(args.print_config_template))
        return

    if not args.config:
        parser.error("-c/--config is required unless --print-config-template is used.")

    workflow_name = args.workflow
    workflow = WORKFLOWS[workflow_name]
    print(f"Genome assembly pipeline version: {__version__}")
    print(f"Workflow: {workflow_name}")

    # Load configuration.
    try:
        with open(args.config, 'r') as f:
            config_object = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Cannot open config file {args.config}")
        sys.exit(1)

    # Create output directory if it does not exist.
    output_dir = os.path.abspath(
        config_get_path(config_object, workflow["output_path"], "output")
        )
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except PermissionError:
            print(f"Cannot create output directory {output_dir}")
            show_singularity_settings(config_object, workflow_name)
            sys.exit(1)

    # Check that required input files exist. The required set is resolved from
    # the workflow and (for generic) the chosen mode / enabled steps. Values may
    # be comma-separated lists (e.g. paired parental Illumina files).
    missing_files = []
    for key in resolve_required_files(workflow_name, config_object):
        configured_value = config_get(config_object, key)
        if not configured_value:
            missing_files.append(key)
            continue
        for configured_path in str(configured_value).split(","):
            configured_path = configured_path.strip()
            if not configured_path:
                continue
            path = os.path.abspath(configured_path)
            if not os.path.exists(path):
                missing_files.append(f"{key} ({path})")

    if missing_files:
        print("Missing or inaccessible required input(s):")
        for missing in missing_files:
            print(f"  - {missing}")
        show_singularity_settings(config_object, workflow_name)
        sys.exit(1)

    # Set the snakefile path (the file is provided via Singularity %files).
    snakefile = pipeline_path(workflow["snakefile"])
    if not os.path.exists(snakefile):
        print(f"Cannot find workflow Snakefile: {snakefile}")
        sys.exit(1)

    # Retrieve the conda environments path from the environment.
    conda_envs_path = os.environ.get('CONDA_ENVS_PATH')
    if not conda_envs_path:
        print("CONDA_ENVS_PATH is not set in the environment.")
        sys.exit(1)

    # Set cache directory for Snakemake to use.
    cache_dir = os.path.join(output_dir, ".cache")

    # Construct the Snakemake command. do not delete file is --keep-going is used
    cmd = [
        "snakemake", "--snakefile", snakefile, "--configfile", args.config,
        "--cores", str(args.threads), "--use-conda", "--conda-prefix",
        conda_envs_path, "--conda-frontend", args.conda_frontend, "--show-failed-logs",
        "--keep-going", "--verbose", "--keep-incomplete"
        ]
    if args.snakemake_args:
        extra_args = []
        for value in args.snakemake_args:
            extra_args.extend(shlex.split(value))
        cmd.extend(extra_args)

    # Update environment with cache directory.
    env = os.environ.copy()
    env['XDG_CACHE_HOME'] = cache_dir

    print("Running Snakemake with the following command:")
    print(" ".join(shlex.quote(part) for part in cmd))

    subprocess.check_call(cmd, env=env)


if __name__ == "__main__":
    main()
