#!/usr/bin/env bash
# Fast, container-free validation of the pipeline sources: Python syntax +
# Snakemake dry-runs across workflows, fixture configs, and every assembly.mode.
# This is the single source of truth for "is a change sound" -- CI and local
# runs both call it, so pre-tag checks match CI exactly.
#
# Usage:  scripts/validate.sh
# Needs `snakemake` on PATH (dry-run only -- no conda envs are built). Override
# the binary with $SNAKEMAKE (e.g. a scratch snakemake-minimal env locally).
set -euo pipefail

cd "$(dirname "$0")/.."
SMK="${SNAKEMAKE:-snakemake}"
command -v "$SMK" >/dev/null 2>&1 || { echo "ERROR: snakemake not found (set \$SNAKEMAKE)"; exit 1; }

fail=0
ok()   { echo "  ok:   $*"; }
bad()  { echo "  FAIL: $*"; fail=1; }
step() { echo; echo "== $* =="; }

dry() {  # dry(label, snakefile, args...)
  local label="$1"; shift
  if "$SMK" --cores 2 --dry-run "$@" >/dev/null 2>/tmp/validate.$$.err; then
    ok "$label"
  else
    bad "$label"; sed 's/^/      /' /tmp/validate.$$.err | tail -8
  fi
  rm -f /tmp/validate.$$.err
}

step "Python syntax + version"
python -m py_compile run_pipeline.py run_mapping.py version.py && ok "py_compile"
python -c "from version import __version__, parse_version; parse_version(__version__); print('  version', __version__)"

step "Build-envs Snakefile enumerates envs"
dry "Snakefile_build_envs parses" -s Snakefile_build_envs

step "Generic workflow: committed fixture configs"
dry "config_generic_ci" -s Snakefile_generic --configfile tests/fixtures/config_generic_ci.yaml
# HiFi fixture: override reads to the committed minimal.fastq so no 142 MB fetch
dry "config_bioinformatics_hifi" -s Snakefile_generic \
  --configfile tests/fixtures/config_bioinformatics_hifi.yaml \
  --config reads='{"ont_fastq":"tests/fixtures/minimal.fastq"}'

step "Generic workflow: assembly.mode matrix + full annotation (dummy inputs)"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/in"
for f in reads.fq hic_R1.fq hic_R2.fq pat.fq mat.fq ref.fa; do : > "$tmp/in/$f"; done

emit_cfg() {  # emit_cfg(mode) -> path
  local mode="$1" cfg="$tmp/$1.yaml"
  {
    echo "reads: { ont_fastq: \"$tmp/in/reads.fq\", hic_r1: \"$tmp/in/hic_R1.fq\", hic_r2: \"$tmp/in/hic_R2.fq\" }"
    echo "contamination_filter: { enabled: false }"
    echo "assembly:"
    echo "  mode: \"$mode\""
    echo "  read_type: \"ont\""
    echo "  trio: { paternal_illumina: \"$tmp/in/pat.fq\", maternal_illumina: \"$tmp/in/mat.fq\" }"
    echo "switch_analysis: { enabled: true }"
    echo "outputs: { output_dir: \"$tmp/out_$mode\" }"
  } > "$cfg"
  echo "$cfg"
}
for m in primary phased hic trio triploid_hic; do
  rm -rf "$tmp/out_$m"
  dry "mode: $m" -s Snakefile_generic --configfile "$(emit_cfg "$m")"
done

# full annotation DAG (every step enabled) on a dummy reference + parents
cat > "$tmp/full.yaml" <<EOF
reads: { ont_fastq: "$tmp/in/reads.fq" }
contamination_filter: { enabled: false }
assembly: { mode: "primary", read_type: "ont" }
reference: { fasta: "$tmp/in/ref.fa" }
quast: { enabled: true }
map_reads_back: { enabled: true }
tidecluster: { enabled: true }
telomere: { enabled: true }
reference_assignment: { enabled: true }
dotplot_vs_reference: { enabled: true }
parental_assignment:
  enabled: true
  paternal_illumina: "$tmp/in/pat.fq"
  maternal_illumina: "$tmp/in/mat.fq"
outputs: { output_dir: "$tmp/out_full" }
EOF
dry "all annotation steps enabled" -s Snakefile_generic --configfile "$tmp/full.yaml"

echo
if [ "$fail" -eq 0 ]; then echo "VALIDATE: OK"; else echo "VALIDATE: FAILED"; exit 1; fi
