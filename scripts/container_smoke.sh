#!/usr/bin/env bash
# Container smoke test: drive a built .sif through the runscript and assert the
# image is complete and produces the expected generic-workflow outputs. Shared
# by release.yml and local pre-release checks so the two cannot drift.
#
# Usage:  scripts/container_smoke.sh <image.sif>
# Env:    THREADS (default 2)
set -euo pipefail

cd "$(dirname "$0")/.."
SIF="${1:?usage: container_smoke.sh <image.sif>}"
[ -s "$SIF" ] || { echo "ERROR: image not found: $SIF"; exit 1; }
RT="$(command -v apptainer || command -v singularity || true)"
[ -n "$RT" ] || { echo "ERROR: neither apptainer nor singularity found on PATH"; exit 1; }
echo "runtime: $RT"
echo "image:   $SIF"

echo; echo "== runscript smoke =="
"$RT" run -B "$PWD" "$SIF" --version
"$RT" run -B "$PWD" "$SIF" --print-config-template generic >/dev/null

echo; echo "== all conda envs baked (no runtime creation) =="
# HOME/XDG_CACHE_HOME must be writable -- the runscript sets these, but a bare
# `exec snakemake` would otherwise die writing its source cache to /root/.cache.
mkdir -p .envcheck/cache
out="$("$RT" exec -B "$PWD" \
  --env HOME="$PWD/.envcheck" --env XDG_CACHE_HOME="$PWD/.envcheck/cache" \
  "$SIF" snakemake -s /opt/pipeline/Snakefile_build_envs --use-conda \
    --conda-prefix /opt/conda/envs --conda-create-envs-only --cores 2 2>&1 || true)"
echo "$out" | tail -40
if echo "$out" | grep -qi "Creating conda environment"; then
  echo "::error::one or more conda envs are missing from the image (would be built at runtime)"
  exit 1
fi
echo "  all conda environments present"

echo; echo "== generic HiFi fixture (primary mode) =="
tests/fixtures/fetch_bioinformatics_hifi_fixture.sh
"$RT" run -B "$PWD" "$SIF" \
  -w generic -c tests/fixtures/config_bioinformatics_hifi.yaml -t "${THREADS:-2}"

echo; echo "== assert expected fixture outputs =="
O=tests/fixtures/output_bioinformatics_hifi
test -s "$O/assembly/assembly.primary.fa"        && echo "  ok: assembly.primary.fa"
test -s "$O/assembly/assembly.outputs.tsv"       && echo "  ok: assembly.outputs.tsv"
test -f "$O/quast/primary/quast_done.txt"        && echo "  ok: quast/primary/quast_done.txt"

echo; echo "CONTAINER-SMOKE: OK"
