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
- Hifiasm (for genome assembly; primary, Hi-C, trio and `--n-hap` phasing)
- Additional tools including minimap2, seqkit, BLAST, Quast, samtools, TideCluster, tidk (telomeres), RagTag (reference scaffolding), yak (parental k-mers), fastp (Illumina QC), and custom Python/R scripts in the `scripts` folder


The following key files are included inside the container:
 
- **Singularity definition file**  (used to build the image)
- **Snakefile**  – defines the original Pisum-oriented/full pipeline rules
- **Snakefile_generic**  – defines the generic lightweight assembly workflow
- **Snakefile_mapping**  – defines the mapping-only workflow
- **config_template.yaml**  – a template for user configuration
- **config_generic_template.yaml**  – a template for generic assembly runs
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
- The default workflow is the original Pisum/full workflow.

To run the generic lightweight assembly workflow, select it explicitly:

```bash
singularity run -B /path/to/data -B $PWD genome_assembly_pipeline.sif \
  -w generic -c config_generic.yaml -t 20
```

To inspect a workflow-specific config template from the container:

```bash
singularity run genome_assembly_pipeline.sif --print-config-template generic
singularity run genome_assembly_pipeline.sif --print-config-template pisum
singularity run genome_assembly_pipeline.sif --version
```

See **Generic workflow: assembly modes and inputs** below for the phasing modes, required input data, and per-haplotype annotation steps the generic workflow supports.

## Workflows

The image ships three independent workflows, selected with `-w`:

| `-w` | Snakefile | Purpose | Config template |
|------|-----------|---------|-----------------|
| `pisum` / `legacy` (default) | `Snakefile` | Original *Pisum*-specific full pipeline (contamination → hifiasm → painting probes → ONT mapping → TideCluster → IGV) | `config_template.yaml` |
| `generic` | `Snakefile_generic` | Species-agnostic assembly + per-haplotype annotation (see below) | `config_generic_template.yaml` |
| `mapping` | `Snakefile_mapping` | Map ONT + HiFi reads back to an existing assembly and run clipping / zero-coverage / TideCluster analysis | `config_mapping_template.yaml` |

## Generic workflow: assembly modes and inputs

The generic workflow assembles with hifiasm and then runs a set of optional
per-haplotype annotation steps. Two things are chosen in the config: the
**assembly mode** (how the genome is phased, and therefore which haplotypes
exist) and which **annotation steps** to run.

### Assembly modes (`assembly.mode`)

| Mode | Phasing | Extra input required | hifiasm flags | Haplotypes produced |
|------|---------|----------------------|---------------|---------------------|
| `primary` | none (collapsed) | – | *(default)* | `primary` |
| `phased` | hifiasm graph phasing | – | *(default)* | `primary`, `hap1`, `hap2` |
| `hic` | Hi-C | Hi-C R1/R2 | `--h1 --h2` (`--enzyme`) | `primary`, `hap1`, `hap2` |
| `trio` | parental k-mers (**both** parents) | paternal + maternal Illumina | `-1 -2` (yak) | `hap1`, `hap2` (parent-labelled) |
| `triploid_hic` | Hi-C, 3 haplotypes | Hi-C R1/R2 | `--n-hap 3 --h1 --h2` | `hap1`, `hap2`, `hap3` |

After assembly a `discover_haplotypes` checkpoint keeps only the haplotypes
hifiasm actually emitted, so the annotation steps fan out over that **realized**
set. For example a homozygous genome run in `phased` mode gracefully falls back
to just `primary` instead of failing. Hi-C phasing separates haplotypes but
**cannot** assign them to a parent — use `parental_assignment` (below) with
parental Illumina for that.

### Input data (config keys)

| Config key | Feeds | Required when |
|------------|-------|---------------|
| `reads.ont_fastq` / `reads.hifi_fastq` | long reads for assembly | always (per `assembly.read_type`) |
| `reads.hic_r1` / `reads.hic_r2` | Hi-C phasing | `mode: hic` or `triploid_hic` |
| `assembly.trio.paternal_illumina` / `maternal_illumina` | trio binning (yak) | `mode: trio` |
| `parental_assignment.paternal_illumina` / `maternal_illumina` | label contigs by parent | `parental_assignment.enabled` |
| `reference.fasta` | RagTag assignment + dotplot | `reference_assignment` or `dotplot_vs_reference` enabled |
| `contamination_filter.plastid_db` / `mitochondrial_db` | plastid/mito read filtering | `contamination_filter.enabled` |
| `tidecluster.library` | custom tandem-repeat annotation library | optional |

Paired-end Illumina is given as a comma-separated pair, e.g.
`paternal_illumina: "pat_R1.fastq.gz,pat_R2.fastq.gz"`. All parental Illumina is
fastp quality/adapter-trimmed before yak counting (`illumina_qc`).

### Per-haplotype annotation steps

Each step is toggled by its own `enabled` flag and runs once per realized haplotype.

| Step | Toggle | Tool | Key output |
|------|--------|------|-----------|
| Contamination filter | `contamination_filter.enabled` | BLAST | filtered reads |
| Illumina QC | `illumina_qc.enabled` | fastp | `data/qc/*.qc.fastq.gz` (+ report) |
| Assembly QC | `quast.enabled` | Quast | `quast/{hap}/` |
| Read-back mapping | `map_reads_back.enabled` | minimap2 | `mapping/reads_to_{hap}.sorted.bam` |
| Tandem repeats + rDNA | `tidecluster.enabled` | TideCluster 1.17 | `analysis/tidecluster/{hap}/` (incl. `tc_rdna.tsv`) |
| Telomeres | `telomere.enabled` | tidk | `analysis/telomere/{hap}/` |
| Reference assignment | `reference_assignment.enabled` | RagTag | `analysis/reference_assignment/{hap}/ragtag.scaffold.*` |
| Dotplot vs reference | `dotplot_vs_reference.enabled` | minimap2 + `paf_dotplot.py` | `analysis/dotplots/{hap}_vs_reference.png` |
| Parental assignment | `parental_assignment.enabled` | yak triobin | `analysis/parental_assignment/parental_assignment.tsv` |

Print the full annotated template with `--print-config-template generic`.

### Examples

```bash
# Primary ONT assembly, telomere + reference dotplot annotation
singularity run -B /path/to/data -B $PWD genome_assembly_pipeline.sif \
  -w generic -c config_generic.yaml -t 20

# Trio-phased HiFi assembly (parental Illumina in the config)
#   assembly: { mode: trio, read_type: hifi,
#               trio: { paternal_illumina: pat_R1,pat_R2, maternal_illumina: mat_R1,mat_R2 } }
singularity run -B /path/to/data -B $PWD genome_assembly_pipeline.sif \
  -w generic -c config_trio.yaml -t 20
```


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

Release images are built automatically by GitHub Actions when a version tag is pushed. The version is defined in `version.py`, and tags must match it exactly:

```bash
python version.py
git tag 0.1.0
git push origin 0.1.0
```

The release workflow builds the SIF, runs the Bioinformatics HiFi fixture inside the built container, pushes the image to GHCR, and creates a GitHub Release. Users can pull a released image with:

```bash
apptainer pull oras://ghcr.io/kavonrtep/ont_genome_assembly_pipeline/sif:0.1.0
apptainer pull oras://ghcr.io/kavonrtep/ont_genome_assembly_pipeline/sif:latest
```

Inspect the embedded version labels with:

```bash
singularity inspect --labels genome_assembly_pipeline_0.1.0.sif
```

To build the container locally, execute the following command (adjust the image name and Singularity path as needed):

```bash
SINGULARITY=$(which singularity)
sudo ionice -c3 $SINGULARITY build images/genome_assembly_pipeline_v9.sif Singularity
sudo ionice -c3 $SINGULARITY build images/genome_assembly_pipeline_v10.sif Singularity
# hifiasm replaced with my fork
sudo ionice -c3 $SINGULARITY build images/genome_assembly_pipeline_v11.sif Singularity
sudo ionice -c3 $SINGULARITY build images/genome_assembly_pipeline_v12.sif Singularity
```

## Testing the Pipeline
```bash
singularity run -B /mnt -B $PWD images/genome_assembly_pipeline_v12.sif -c config.yaml -t 6

singularity run -B /mnt -B $PWD images/genome_assembly_pipeline_v12.sif \
  -w generic -c config_generic.yaml -t 6 -S "--dry-run"
```

The release fixture uses the PacBio HiFi training data from `kavonrtep/bioinformatics`:

```bash
tests/fixtures/fetch_bioinformatics_hifi_fixture.sh
singularity run -B $PWD images/genome_assembly_pipeline_v12.sif \
  -w generic -c tests/fixtures/config_bioinformatics_hifi.yaml -t 2
```

To exercise the whole generic workflow (assembly + every annotation step) on the
*S. aureus* NCTC 8325 HiFi dataset against its real reference, prepare the
fixture data and run `config_generic_fulltest.yaml`:

```bash
tests/fixtures/prepare_fulltest.sh   # fetches reads + reference, makes parental subsets
singularity run -B $PWD images/genome_assembly_pipeline_v12.sif \
  -w generic -c tests/fixtures/config_generic_fulltest.yaml -t 8
```


## Troubleshooting 
- **Input Files Not Found:** 

Verify that the paths in your `config.yaml` are correct and that the directories are correctly bound (using `-B`) when running the container.
- **Permission Issues:** 
Ensure that the output directory is writable. If you receive permission errors, adjust your bind mounts or output directory paths.
- 06
