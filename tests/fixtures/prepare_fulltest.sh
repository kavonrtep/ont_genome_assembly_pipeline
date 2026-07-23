#!/usr/bin/env bash
# Fetch/derive everything config_generic_fulltest.yaml needs (data is gitignored,
# so this makes the full end-to-end fixture reproducible from scratch):
#   * S. aureus NCTC 8325 HiFi reads (~25x)
#   * the matching GCF_000013425.1 reference genome
#   * two synthetic "parental" read subsets for the yak/triobin branch
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="$here/data"
mkdir -p "$data_dir"

"$here/fetch_bioinformatics_hifi_fixture.sh" "$data_dir"

ref="$data_dir/GCF_000013425.1_ASM1342v1_genomic.fna"
if [[ ! -s "$ref" ]]; then
    curl -L --fail --retry 3 --retry-delay 5 -o "$ref" \
      "https://raw.githubusercontent.com/kavonrtep/bioinformatics/master/data/genome_assembly/GCF_000013425.1_ASM1342v1_genomic.fna"
fi
echo "reference: $ref"

"$here/make_parents.sh" "$data_dir"

echo "Fixture ready. Run with run_pipeline.py -w generic -c tests/fixtures/config_generic_fulltest.yaml"
