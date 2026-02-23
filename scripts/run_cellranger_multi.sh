#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------
# STAR-based "primary" for 10x GEX (Perturb-seq)
# Replaces CellRanger multi with STARsolo.
# ---------------------------------------

ROOT="${ROOT:-/work}"

RUN_ID="${1:-sc3_crispr_10k}"
GEX_FASTQ_DIR="${2:-$ROOT/data/fastq/SC3_v3_NextGem_DI_CRISPR_10K_fastqs/SC3_v3_NextGem_DI_CRISPR_10K_gex_fastqs}"
REF_DIR="${3:-$ROOT/data/ref/refdata-gex-GRCh38-2020-A}"
OUT_DIR="${4:-$ROOT/outputs/primary/$RUN_ID}"
THREADS="${THREADS:-8}"

# 10x chemistry for this dataset: NextGem 3' v3.1 / v3 (often treated similarly in STARsolo presets)
CHEMISTRY="${CHEMISTRY:-SC3Pv3}"   # STARsolo preset name; we set explicit lengths below to be safe

# Where to find whitelists (10x barcodes list)
# CellRanger ref has barcodes whitelist usually under: lib/python/cellranger/barcodes/
# We don't have that; so we ship a whitelist file ourselves or download from 10x.
# For 10x 3' v3: 3M-february-2018.txt.gz is used commonly.
WHITELIST="${WHITELIST:-$ROOT/data/ref/10x_whitelists/3M-february-2018.txt.gz}"

# Read structure: For 10x GEX, usually:
#   R1 = cell barcode + UMI (e.g., 16bp CB + 12bp UMI = 28bp)
#   R2 = cDNA sequence
CB_LEN="${CB_LEN:-16}"
UMI_LEN="${UMI_LEN:-12}"

mkdir -p "$OUT_DIR"
mkdir -p "$(dirname "$WHITELIST")"

log() { echo "==> $*"; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "❌ Missing required command: $1"
    exit 1
  }
}

# ---------------------------------------
# Checks
# ---------------------------------------
need_cmd STAR
need_cmd pigz

[[ -d "$GEX_FASTQ_DIR" ]] || { echo "❌ GEX_FASTQ_DIR not found: $GEX_FASTQ_DIR"; exit 1; }
[[ -d "$REF_DIR/star" ]] || { echo "❌ STAR index not found at: $REF_DIR/star"; exit 1; }
[[ -f "$REF_DIR/fasta/genome.fa" ]] || { echo "❌ genome.fa missing under: $REF_DIR/fasta"; exit 1; }
[[ -f "$REF_DIR/genes/genes.gtf" ]] || { echo "❌ genes.gtf missing under: $REF_DIR/genes"; exit 1; }

# ---------------------------------------
# Whitelist handling
# ---------------------------------------
if [[ ! -f "$WHITELIST" ]]; then
  log "Whitelist not found. Downloading 10x v3 whitelist (3M-february-2018) ..."
  # This is the standard 10x whitelist widely used for v2/v3 chemistries.
  # If you prefer a different one, override WHITELIST env var.
  wget -q -O "$WHITELIST" "https://cf.10xgenomics.com/supp/cell-exp/3M-february-2018.txt.gz"
fi

# ---------------------------------------
# Discover FASTQs (Illumina naming)
# ---------------------------------------
# Expecting files like: *_S1_L00X_R1_001.fastq.gz and *_R2_001.fastq.gz
mapfile -t R1S < <(ls -1 "$GEX_FASTQ_DIR"/*_R1_001.fastq.gz 2>/dev/null | sort || true)
mapfile -t R2S < <(ls -1 "$GEX_FASTQ_DIR"/*_R2_001.fastq.gz 2>/dev/null | sort || true)

if [[ ${#R1S[@]} -eq 0 || ${#R2S[@]} -eq 0 ]]; then
  echo "❌ Could not find R1/R2 FASTQs in: $GEX_FASTQ_DIR"
  echo "   Found:"
  ls -lh "$GEX_FASTQ_DIR" | tail -n 30 || true
  exit 1
fi

if [[ ${#R1S[@]} -ne ${#R2S[@]} ]]; then
  echo "❌ Mismatch number of R1 and R2 files:"
  echo "   R1=${#R1S[@]} R2=${#R2S[@]}"
  exit 1
fi

# STAR accepts comma-separated lists for multiple lanes
R1_CSV="$(IFS=,; echo "${R1S[*]}")"
R2_CSV="$(IFS=,; echo "${R2S[*]}")"

log "RUN_ID=$RUN_ID"
log "GEX_FASTQ_DIR=$GEX_FASTQ_DIR"
log "REF_DIR=$REF_DIR"
log "OUT_DIR=$OUT_DIR"
log "THREADS=$THREADS"
log "WHITELIST=$WHITELIST"
log "Detected lanes: ${#R1S[@]}"

# ---------------------------------------
# Run STARsolo (GEX)
# Outputs: matrix.mtx + barcodes.tsv + features.tsv (10x-like)
# ---------------------------------------
log "Running STARsolo (GEX) ..."

STAR \
  --genomeDir "$REF_DIR/star" \
  --readFilesIn "$R2_CSV" "$R1_CSV" \
  --readFilesCommand pigz -dc \
  --runThreadN "$THREADS" \
  --outFileNamePrefix "$OUT_DIR/" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI AS nM CR UR CB UB GX GN sS sQ \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$WHITELIST" \
  --soloCBlen "$CB_LEN" \
  --soloUMIlen "$UMI_LEN" \
  --soloBarcodeReadLength 0 \
  --soloFeatures Gene \
  --soloUMIdedup 1MM_All \
  --soloCellFilter EmptyDrops_CR \
  --limitBAMsortRAM 12000000000

log "STARsolo finished."

# STARsolo writes under: Solo.out/
if [[ -d "$OUT_DIR/Solo.out" ]]; then
  log "STARsolo outputs:"
  find "$OUT_DIR/Solo.out" -maxdepth 4 -type f | sed 's#^#  - #' | head -n 80
else
  log "WARN: Solo.out not found. Check STAR logs under $OUT_DIR/Log.final.out"
fi

log "Done."
log "Primary outputs in: $OUT_DIR"
