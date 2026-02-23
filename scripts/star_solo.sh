#!/usr/bin/env bash
set -euo pipefail

R1="${1:?R1 fastq required}"
R2="${2:?R2 fastq required}"
STAR_INDEX="${3:?STAR index dir required}"
OUTDIR="${4:?output dir required}"
SAMPLE="${5:?sample name required}"
THREADS="${6:-8}"

SOLO_CB_WHITELIST="${SOLO_CB_WHITELIST:-None}"

echo "[star_solo.sh] R1=${R1}"
echo "[star_solo.sh] R2=${R2}"
echo "[star_solo.sh] STAR_INDEX=${STAR_INDEX}"
echo "[star_solo.sh] OUTDIR=${OUTDIR}"
echo "[star_solo.sh] SAMPLE=${SAMPLE}"
echo "[star_solo.sh] THREADS=${THREADS}"
echo "[star_solo.sh] SOLO_CB_WHITELIST=${SOLO_CB_WHITELIST}"

mkdir -p "${OUTDIR}"

STAR \
  --genomeDir "${STAR_INDEX}" \
  --readFilesIn "${R1}" "${R2}" \
  --readFilesCommand zcat \
  --runThreadN "${THREADS}" \
  --outFileNamePrefix "${OUTDIR}/${SAMPLE}." \
  --outSAMtype BAM SortedByCoordinate \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "${SOLO_CB_WHITELIST}" \
  --soloCBstart 1 --soloCBlen 16 \
  --soloUMIstart 17 --soloUMIlen 12 \
  --soloBarcodeReadLength 0 \
  --soloFeatures Gene
