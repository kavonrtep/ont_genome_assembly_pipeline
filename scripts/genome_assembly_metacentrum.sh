#!/bin/bash
#PBS -N ont_genome_assembly
#PBS -l select=1:ncpus=32:mem=256gb:scratch_local=600gb
#PBS -q elixircz@pbs-m1.metacentrum.cz
#PBS -l walltime=120:00:00
#PBS -j oe
#PBS -m bae

## CONFIGURATION - use absolute paths for input files and Singularity image
SINGULARITY_IMAGE="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/genome_assembly_pipeline_v1.sif"

# Input files (absolute paths)
OXFORD_READS="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/IPIP200579_ULK250311_dorado_0.9.1_dupl_sup_5mC_filt_Q12_m30k.fastq"
PLASTID_DB="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/P.sativum_NC_014057.1_plast.fasta"
MITO_DB="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/P.sativum_elatius_OQ078753.1_mit.fasta"
TANDEM_REPEATS="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/FabTR_all_sequences_231010.db.RM_format"
PAINTING_PROBES="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/Cameor_v2_release_2_chromosomes_only.fasta_all.NGSfilter_CamIllumina.selected.CLEAN.fasta"

# Optional for generating dotplots
ASSEMBLY_FAI_REF="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/cameor_v2.fasta.fai"    # leave empty if not used
BLAST_PROBES_REF="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/CAMv2r3_x_oligos_CAMv2r2.blast_out"  # leave empty if not used

# Directories (output will be created inside SCRATCHDIR and then copied back)
OUTPUT_DIR="/auto/brno2/home/kavonrtep/IPIP200579_ULK250311/hifiasm"    # permanent storage location
QUAST_DIR="quast"         # relative directory name for Quast results
ANALYSIS_DIR="analysis"   # relative directory name for additional analysis

# Pipeline parameters
MIN_LIMIT=80000
MAX_LIMIT=120000
TRIM_START=100
TRIM_END=10
## END OF CONFIGURATION

# Copy Singularity image and input files to SCRATCHDIR
SINGULARITY_IMAGE_BASENAME=$(basename "$SINGULARITY_IMAGE")

rsync -avt "$SINGULARITY_IMAGE" "$SCRATCHDIR/"
rsync -avt "$OXFORD_READS" "$SCRATCHDIR/$(basename "$OXFORD_READS")"
rsync -avt "$PLASTID_DB" "$SCRATCHDIR/$(basename "$PLASTID_DB")"
rsync -avt "$MITO_DB" "$SCRATCHDIR/$(basename "$MITO_DB")"
rsync -avt "$TANDEM_REPEATS" "$SCRATCHDIR/$(basename "$TANDEM_REPEATS")"
rsync -avt "$PAINTING_PROBES" "$SCRATCHDIR/$(basename "$PAINTING_PROBES")"

if [ -n "$ASSEMBLY_FAI_REF" ]; then
    rsync -avt "$ASSEMBLY_FAI_REF" "$SCRATCHDIR/$(basename "$ASSEMBLY_FAI_REF")"
fi

if [ -n "$BLAST_PROBES_REF" ]; then
    rsync -avt "$BLAST_PROBES_REF" "$SCRATCHDIR/$(basename "$BLAST_PROBES_REF")"
fi

cd "$SCRATCHDIR"

# Generate config.yaml for the pipeline
cat > config.yaml <<EOF
oxford_nanopore_reads: $(basename "$OXFORD_READS")
plastid_db: $(basename "$PLASTID_DB")
mitochondrial_db: $(basename "$MITO_DB")
tandem_repeats: $(basename "$TANDEM_REPEATS")
painting_probes: $(basename "$PAINTING_PROBES")
EOF

# Append optional keys if provided
if [ -n "$ASSEMBLY_FAI_REF" ]; then
    echo "assembly_fai_ref: $(basename "$ASSEMBLY_FAI_REF")" >> config.yaml
fi

if [ -n "$BLAST_PROBES_REF" ]; then
    echo "blast_probes_ref: $(basename "$BLAST_PROBES_REF")" >> config.yaml
fi

cat >> config.yaml <<EOF
quast_dir: $QUAST_DIR
analysis_dir: $ANALYSIS_DIR
output_dir: output
min_limit: $MIN_LIMIT
max_limit: $MAX_LIMIT
trim_start: $TRIM_START
trim_end: $TRIM_END
EOF

# Prepare temporary directory for Singularity to use
tmpdir="$SCRATCHDIR/tmp"
mkdir -p "$tmpdir"
export TMPDIR="$tmpdir"

# Run the pipeline inside the Singularity container
singularity run -B "$SCRATCHDIR" --env TMPDIR="$tmpdir" "$SINGULARITY_IMAGE_BASENAME" -c config.yaml -t "$PBS_NCPUS"

# Copy results from SCRATCHDIR back to permanent storage
mkdir -p "$OUTPUT_DIR"
rsync -avt "$SCRATCHDIR/output" "$OUTPUT_DIR/"

# Optionally, copy the config and log files for record keeping
rsync -avt "$SCRATCHDIR/config.yaml" "$OUTPUT_DIR/"
rsync -avt "$SCRATCHDIR/.snakemake" "$OUTPUT_DIR/"

# Archive the output directory
zip -y -fz -r "$SCRATCHDIR/output.zip" output
rsync -avt "$SCRATCHDIR/output.zip" "$OUTPUT_DIR/"
