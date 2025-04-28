#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh
conda activate pysam
exec python /opt/pipeline/scripts/extract_reads_from_bam3.py "$@"
