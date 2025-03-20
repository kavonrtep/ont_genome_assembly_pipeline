#!/usr/bin/env python3
"""
This script reads a user-provided config file and launches the Snakemake pipeline.
It is designed to run within your Singularity container.
"""

import argparse
import os
import subprocess
import sys
import yaml


def show_singularity_settings(config_object):
    """
    Print the Singularity bind options based on directories from required files.
    """
    required_keys = ["oxford_nanopore_reads", "painting_probes", "assembly_fai_ref",
                     "blast_probes_ref", "tandem_repeats"]
    dirs = set()
    # Always include the output directory.
    dirs.add(os.path.abspath(config_object['output_dir']))
    for key in required_keys:
        if key in config_object:
            dirs.add(os.path.dirname(os.path.abspath(config_object[key])))
    bind_string = " ".join([f"-B {d}" for d in dirs])
    print("Run Singularity with the following bind options:")
    print(f"singularity run {bind_string} ...")


def main():
    # Default config template file inside container.
    config_template = "/opt/pipeline/config.yaml"
    try:
        with open(config_template, 'r') as f:
            config_string = f.read()
    except Exception:
        config_string = "Could not read default config template."

    parser = argparse.ArgumentParser(
        description="Pipeline for assembly and analysis using Snakemake within a "
                    "Singularity container.",
        epilog=f"Example config.yaml file:\n\n{config_string}",
        formatter_class=argparse.RawDescriptionHelpFormatter
        )
    parser.add_argument('-c', '--config', required=True, help='Path to config file')
    parser.add_argument(
        '-t', '--threads', required=False, default=2, type=int,
        help='Number of threads to use'
        )
    parser.add_argument(
        '-S', '--snakemake_args', type=str, nargs='?', default="",
        help='Additional snakemake arguments (e.g., "--dry-run --reason"). Do not '
             'include options that are set by the script.'
        )
    args = parser.parse_args()

    # Load configuration.
    try:
        with open(args.config, 'r') as f:
            config_object = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Cannot open config file {args.config}")
        sys.exit(1)

    # Create output directory if it does not exist.
    output_dir = os.path.abspath(config_object.get('output_dir', 'output'))
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except PermissionError:
            print(f"Cannot create output directory {output_dir}")
            show_singularity_settings(config_object)
            sys.exit(1)

    # Check that required input files exist.
    required_files = ["oxford_nanopore_reads", "painting_probes", "assembly_fai_ref",
                      "blast_probes_ref", "tandem_repeats"]
    for key in required_files:
        if key in config_object:
            path = os.path.abspath(config_object[key])
            if not os.path.exists(path):
                print(
                    f"Required file for '{key}' at {path} does not exist or is not "
                    f"accessible."
                    )
                show_singularity_settings(config_object)
                sys.exit(1)

    # Set the snakefile path (the file is provided via Singularity %files).
    snakefile = "/opt/pipeline/Snakefile"

    # Retrieve the conda environments path from the environment.
    conda_envs_path = os.environ.get('CONDA_ENVS_PATH')
    if not conda_envs_path:
        print("CONDA_ENVS_PATH is not set in the environment.")
        sys.exit(1)

    # Set cache directory for Snakemake to use.
    cache_dir = os.path.join(output_dir, ".cache")

    # Construct the Snakemake command.
    cmd = (f"snakemake --snakefile {snakefile} --configfile {args.config} "
           f"--cores {args.threads} --use-conda --conda-prefix {conda_envs_path} "
           f"--conda-frontend mamba --show-failed-logs {args.snakemake_args}")

    # Update environment with cache directory.
    env = os.environ.copy()
    env['XDG_CACHE_HOME'] = cache_dir

    print("Running Snakemake with the following command:")
    print(cmd)

    subprocess.check_call(cmd, shell=True, env=env)


if __name__ == "__main__":
    main()
