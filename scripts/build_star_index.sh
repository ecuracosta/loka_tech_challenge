#!/usr/bin/env bash
set -euo pipefail

# Build STAR genome index.
#
# Usage:
#   build_star_index.sh <GENOME_FASTA> <GTF> <OUT_INDEX_DIR> <THREADS> [SJDB_OVERHANG]
#
# Notes:
# - sjdbOverhang should be (readLength - 1). If read length is unknown, 100 is a reasonable default for 101bp reads.

GENOME_FASTA="${1:?Missing genome FASTA}"
GTF="${2:?Missing GTF}"
OUT_DIR="${3:?Missing output index dir}"
THREADS="${4:-8}"
SJDB_OVERHANG="${5:-100}"

mkdir -p "${OUT_DIR}"

# Skip only if key index files exist
for f in Genome SA SAindex genomeParameters.txt; do
  if [[ ! -s "${OUT_DIR}/${f}" ]]; then
    missing=1
  fi
done

if [[ "${missing:-0}" -eq 0 ]]; then
  echo "STAR index looks complete at ${OUT_DIR}. Skipping build."
  exit 0
fi

# If a previous attempt produced partial files, keep them for inspection.
# STAR does not reliably resume, but this avoids silently redoing work unless you explicitly clean.
if [[ -f "${OUT_DIR}/.build_failed" ]]; then
  echo "WARNING: previous STAR index build failed (marker found: ${OUT_DIR}/.build_failed)."
  echo "Inspect logs/files in ${OUT_DIR}, or delete the directory to rebuild from scratch."
  exit 2
fi

echo "Building STAR index..."
echo "  FASTA: ${GENOME_FASTA}"
echo "  GTF:   ${GTF}"
echo "  OUT:   ${OUT_DIR}"

STAR \
  --runThreadN "${THREADS}" \
  --runMode genomeGenerate \
  --genomeDir "${OUT_DIR}" \
  --genomeFastaFiles "${GENOME_FASTA}" \
  --sjdbGTFfile "${GTF}" \
  --sjdbOverhang "${SJDB_OVERHANG}" \
  --genomeSAsparseD 4

echo "Done. Index created at ${OUT_DIR}"
