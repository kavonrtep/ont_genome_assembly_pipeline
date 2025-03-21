# Genome Assembly and Analysis Pipeline

This pipeline performs genome assembly from Oxford Nanopore reads, contamination detection, quality assessment, and several downstream analyses. The workflow is orchestrated by Snakemake and is executed within a Singularity container.

## Requirements

- **Singularity:** The container runtime is required to run the pipeline. You can install Singularity using conda as follows:
  
 ```bash
  conda create -n singularity3 -c conda-forge "singularity>=3.6"
  conda activate singularity3
```

- **Input Data:**  Ensure that your input files (e.g., FASTQ, FASTA, BLAST databases) and configuration file are accessible to the container via proper bind mounts.


## Container and Files 

The container is built from a Singularity definition file that installs:
- Snakemake (workflow management)
- Hifiasm (for genome assembly)
- Dorado (for read mapping)
- Additional tools including seqkit, BLAST, Quast, samtools, and custom Python/R scripts in the `scripts` folder


The following key files are included inside the container:
 
- **Singularity definition file**  (used to build the image)
- **Snakefile**  – defines the pipeline rules
- **config_template.yaml**  – a template for user configuration
- **run_pipeline.py**  – the script that launches the pipeline


## Configuration File 

An example `config.yaml` file is shown below. Adjust file paths as necessary:

```yaml
# Input files
oxford_nanopore_reads: "data/read_all_sample.fastq"      
plastid_db: "data/P.fulvum_MW160430.1_plast.fasta"                     
mitochondrial_db: "data/Pisum_fulvum_NC_059792.1_min.fasta"         
tandem_repeats: "data/FabTR_all_sequences_231010.db.RM_format"
painting_probes: "data/Cameor_v2_release_2_chromosomes_only.fasta_all.NGSfilter_CamIllumina.selected.CLEAN.fasta"
assembly_fai_ref: "data/JI1006_fulvum_v2.1.fasta.fai"
blast_probes_ref: "data/JI1006_v2.1_x_oligos_CAMv2r2.blast_out"

# Directories
quast_dir: "quast"
analysis_dir: "analysis"
output_dir: "output"
```


## Quick Start 

A Singularity image (e.g. `genome_assembly_pipeline.sif`) can be built from the provided Singularity definition file. Once built (or downloaded), run the pipeline with the following command:


```bash
singularity run -B /path/to/data -B $PWD genome_assembly_pipeline.sif -c config.yaml -t 20
```

 
- The `-B` flag binds host directories so that input/output files are accessible inside the container.
- The `-c` option specifies the path to your configuration file.
- The `-t` option sets the number of threads to use during execution.


## Running on a Cluster 
TODO

If you plan to run the pipeline on a cluster (e.g., Metacentrum), you can use the provided script in `scripts/run_pipeline_metacentrum.sh` (adjust file paths and parameters as needed) to submit jobs.

## Output Structure 

After successful execution, the output directory (as specified in `config.yaml`) will contain:
 
- **Assembly Outputs:** 
  - Final assembly file (e.g. `hifiasm_assembly.bp.p_ctg.gfa.fasta`)
- **Quality Assessment:** 
  - A `quast` folder with Quast results and a marker file (e.g. `quast_done.txt`)
 
- **Downstream Analyses:** 
  - Subdirectories within the `analysis` folder for painting probes mapping, Oxford Nanopore read mapping, clipping information, contig pairs, and TideCluster results.


## Building the Singularity Container 

To build the container locally, execute the following command (adjust the image name and Singularity path as needed):

```bash
SINGULARITY=$(which singularity)
sudo ionice -c3 $SINGULARITY build images/genome_assembly_pipeline_v1.sif Singularity
```



## Troubleshooting 
- **Input Files Not Found:** 

Verify that the paths in your `config.yaml` are correct and that the directories are correctly bound (using `-B`) when running the container.
- **Permission Issues:** 
Ensure that the output directory is writable. If you receive permission errors, adjust your bind mounts or output directory paths.
- **Pipeline Errors:** 
Check the log files in the `.snakemake/log/` directory for detailed error messages and troubleshooting hints.


