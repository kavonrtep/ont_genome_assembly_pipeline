#!/usr/bin/env bash
# Generate the two SYNTHETIC "parental" read sets used by
# config_generic_fulltest.yaml to exercise the parental_assignment (yak triobin)
# branch end to end. These are simply two disjoint 1500-read subsets of the
# S. aureus HiFi fixture reads -- they verify the yak/triobin/merge mechanics
# and output format, NOT real trio biology.
set -euo pipefail

data_dir="${1:-tests/fixtures/data}"
reads="$data_dir/pacbio_hifi_reads25.fastq"

if [[ ! -s "$reads" ]]; then
    echo "Missing $reads -- run tests/fixtures/fetch_bioinformatics_hifi_fixture.sh first." >&2
    exit 1
fi

head -n 6000 "$reads" > "$data_dir/parent_pat.fastq"          # reads 1..1500
sed -n '6001,12000p' "$reads" > "$data_dir/parent_mat.fastq"  # reads 1501..3000
echo "Wrote $data_dir/parent_pat.fastq and $data_dir/parent_mat.fastq (1500 reads each)."
