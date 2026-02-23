set -euo pipefail

# Usage:
#   star_align.sh <R1.fastq(.gz)> <R2.fastq(.gz)> <STAR_INDEX_DIR> <OUT_DIR> <SAMPLE_ID> <THREADS>
#
# Outputs (in OUT_DIR):
#   <SAMPLE_ID>.Aligned.sortedByCoord.out.bam
#   <SAMPLE_ID>.Log.final.out
#   <SAMPLE_ID>.Log.out
#   <SAMPLE_ID>.Log.progress.out
#   <SAMPLE_ID>.SJ.out.tab

R1="${1:?Missing R1 FASTQ}"
R2="${2:?Missing R2 FASTQ}"
STAR_INDEX_DIR="${3:?Missing STAR index dir}"
OUT_DIR="${4:?Missing output dir}"
SAMPLE_ID="${5:?Missing sample id}"
THREADS="${6:-8}"

mkdir -p "${OUT_DIR}"

# Decide whether input is gzipped
READ_CMD="cat"
if [[ "${R1}" == *.gz ]]; then
  READ_CMD="zcat"
fi

STAR \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX_DIR}" \
  --readFilesIn "${R1}" "${R2}" \
  --readFilesCommand "${READ_CMD}" \
  --outFileNamePrefix "${OUT_DIR}/${SAMPLE_ID}." \
  --outSAMtype BAM SortedByCoordinate

# STAR writes: <prefix>Aligned.sortedByCoord.out.bam and logs with same prefix
EOF