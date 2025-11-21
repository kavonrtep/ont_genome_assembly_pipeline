#!/usr/bin/env python3
"""
This script reads a user-provided config file and launches the Snakemake pipeline
for read mapping and analysis. It is designed to run within your Singularity container.
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
    required_keys = ["genome_assembly", "oxford_nanopore_reads", "hifi_reads",
                     "tandem_repeats"]
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
        config_string = """# Example config.yaml for read mapping pipeline
# Required input files
genome_assembly: "/path/to/your/genome_assembly.fasta"
oxford_nanopore_reads: "/path/to/ont_reads.fastq"
hifi_reads: "/path/to/hifi_reads.fastq.gz"


# Output directory
output_dir: "results"

# Optional parameters (uncomment and modify if needed)
# min_mapq: 1  # Minimum mapping quality for filtering
"""

    parser = argparse.ArgumentParser(
        description="Pipeline for read mapping and analysis using Snakemake within a "
                    "Singularity container. Maps ONT and HiFi reads to a pre-assembled genome.",
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
    required_files = ["genome_assembly", "oxford_nanopore_reads", "hifi_reads",
                      "tandem_repeats"]
    missing_files = []

    for key in required_files:
        if key not in config_object:
            print(f"Required key '{key}' is missing from config file.")
            missing_files.append(key)
        else:
            path = os.path.abspath(config_object[key])
            if not os.path.exists(path):
                print(
                    f"Required file for '{key}' at {path} does not exist or is not "
                    f"accessible."
                    )
                missing_files.append(key)

    if missing_files:
        print(f"\nMissing or inaccessible files: {', '.join(missing_files)}")
        show_singularity_settings(config_object)
        sys.exit(1)

    # Validate file formats (basic checks)
    genome_path = os.path.abspath(config_object["genome_assembly"])
    if not (genome_path.endswith('.fasta') or genome_path.endswith(
            '.fa'
            ) or genome_path.endswith(
        '.fna'
        )):
        print(
            f"Warning: Genome assembly file '{genome_path}' does not have a "
            f"standard FASTA extension (.fasta, .fa, .fna)"
            )

    ont_reads_path = os.path.abspath(config_object["oxford_nanopore_reads"])
    if not (ont_reads_path.endswith('.fastq') or ont_reads_path.endswith('.fq')):
        print(
            f"Warning: ONT reads file '{ont_reads_path}' does not have a "
            f"standard FASTQ extension (.fastq, .fq)"
            )

    hifi_reads_path = os.path.abspath(config_object["hifi_reads"])
    if not (hifi_reads_path.endswith('.fastq') or hifi_reads_path.endswith(
            '.fq'
            ) or hifi_reads_path.endswith(
        '.fastq.gz'
        ) or hifi_reads_path.endswith('.fq.gz')):
        print(
            f"Warning: HiFi reads file '{hifi_reads_path}' does not have a "
            f"standard FASTQ extension (.fastq, .fq, .fastq.gz, .fq.gz)"
            )

    # Set the snakefile path (the file is provided via Singularity %files).
    snakefile = "/opt/pipeline/Snakefile_mapping"

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
           f"--conda-frontend mamba --show-failed-logs {args.snakemake_args}"
           " --keep-going --verbose --keep-incomplete")

    # Update environment with cache directory.
    env = os.environ.copy()
    env['XDG_CACHE_HOME'] = cache_dir

    print("=" * 60)
    print("READ MAPPING AND ANALYSIS PIPELINE")
    print("=" * 60)
    print(f"Config file: {args.config}")
    print(f"Threads: {args.threads}")
    print(f"Output directory: {output_dir}")
    print(f"Genome assembly: {config_object['genome_assembly']}")
    print(f"ONT reads: {config_object['oxford_nanopore_reads']}")
    print(f"HiFi reads: {config_object['hifi_reads']}")
    print(f"Tandem repeats library: {config_object['tandem_repeats']}")
    print("=" * 60)
    print("\nRunning Snakemake with the following command:")
    print(cmd)
    print("=" * 60)

    try:
        subprocess.check_call(cmd, shell=True, env=env)
        print("\n" + "=" * 60)
        print("PIPELINE COMPLETED SUCCESSFULLY!")
        print("=" * 60)
        print(f"Results can be found in: {output_dir}")
        print("\nKey output files:")
        print(
            f"  - ONT mapping: {output_dir}/analysis/mapping_ONT_reads/ont_mapped.sorted.bam"
            )
        print(
            f"  - HiFi mapping: {output_dir}/analysis/mapping_HiFi_reads/hifi_mapped.sorted.bam"
            )
        print(f"  - Clipping analysis: {output_dir}/analysis/clipping_info/")
        print(f"  - Contig pairs: {output_dir}/analysis/contig_pairs/")
        print(f"  - TideCluster results: {output_dir}/analysis/tidecluster/")
        print(f"  - IGV session: {output_dir}/analysis/igv/igv.xml")
    except subprocess.CalledProcessError as e:
        print(f"\nPipeline failed with error code {e.returncode}")
        print("Check the log files for detailed error information.")
        sys.exit(1)


if __name__ == "__main__":
    main()