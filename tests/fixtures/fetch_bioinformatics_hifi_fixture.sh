#!/usr/bin/env bash
set -euo pipefail

out_dir="${1:-tests/fixtures/data}"
mkdir -p "$out_dir"

url="https://raw.githubusercontent.com/kavonrtep/bioinformatics/master/data/genome_assembly/pacbio_hifi_reads25.fastq.gz"
gz="$out_dir/pacbio_hifi_reads25.fastq.gz"
out="$out_dir/pacbio_hifi_reads25.fastq"

if [[ -s "$out" ]]; then
    echo "Fixture already exists: $out"
    exit 0
fi

curl -L --fail --retry 3 --retry-delay 5 -o "$gz" "$url"
gzip -dc "$gz" > "$out"
ls -lh "$out"
