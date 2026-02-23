#!/usr/bin/env bash
set -euo pipefail

# -------------------------
# Defaults (override via args/env)
# -------------------------
ROOT="${ROOT:-/work}"

# Input FASTQs (CRISPR library: combined R1/R2)
R1="${1:-$ROOT/data/fastq/SC3_v3_NextGem_DI_CRISPR_10K_crispr_S1_combined_R1_001.fastq.gz}"
R2="${2:-$ROOT/data/fastq/SC3_v3_NextGem_DI_CRISPR_10K_crispr_S1_combined_R2_001.fastq.gz}"

# Cell barcodes whitelist + feature reference (2-col TSV: feature_id \t sequence)
WHITELIST="${3:-$ROOT/data/barcode/filtered_feature_bc_matrix/barcodes.tsv.gz}"
FEATURE_TSV="${4:-$ROOT/data/ref/SC3_v3_NextGem_DI_CRISPR_10K_feature_ref.tsv}"

# QC params
R1_C="${5:-0,16}"          # cell barcode range on read1
NREADS="${6:-20000}"       # number of read pairs sampled by fba qc
THREADS="${7:-8}"
RUN_ID="${8:-qc_$(date +%Y%m%d_%H%M%S)}"

# Output layout
OUT_ROOT="${OUT_ROOT:-$ROOT/outputs}"
OUT_DIR="$OUT_ROOT/qc/$RUN_ID"
FBA_QC_DIR="$OUT_DIR/fba_qc"
TABLE_DIR="$OUT_DIR/tables"
FIG_DIR="$OUT_DIR/figures"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$FBA_QC_DIR" "$TABLE_DIR" "$FIG_DIR" "$LOG_DIR"

echo "==> RUN_ID=$RUN_ID"
echo "==> R1=$R1"
echo "==> R2=$R2"
echo "==> WHITELIST=$WHITELIST"
echo "==> FEATURE_TSV=$FEATURE_TSV"
echo "==> OUT_DIR=$OUT_DIR"

# -------------------------
# 1) Run fba qc
# -------------------------
echo "==> Running: fba qc"
fba qc \
  -1 "$R1" \
  -2 "$R2" \
  -w "$WHITELIST" \
  -f "$FEATURE_TSV" \
  -r1_c "$R1_C" \
  -n "$NREADS" \
  -t "$THREADS" \
  --output_directory "$FBA_QC_DIR" \
  2>&1 | tee "$LOG_DIR/fba_qc.log"

# Expected by tutorial: qc/feature_barcoding_output.tsv.gz inside output_directory
QC_TSV_GZ="$FBA_QC_DIR/qc/feature_barcoding_output.tsv.gz"
if [[ ! -f "$QC_TSV_GZ" ]]; then
  echo "❌ Expected QC TSV not found: $QC_TSV_GZ"
  echo "   Listing $FBA_QC_DIR:"
  find "$FBA_QC_DIR" -maxdepth 3 -type f | sed 's/^/   - /'
  exit 2
fi

# Copy key artifact into tables/
cp -f "$QC_TSV_GZ" "$TABLE_DIR/feature_barcoding_output.tsv.gz"

# -------------------------
# 2) Python plots + summary tables
# -------------------------
echo "==> Running: python QC plots"
python3 "$ROOT/qc/plot_qc.py" \
  --r1 "$R1" \
  --r2 "$R2" \
  --qc_tsv_gz "$QC_TSV_GZ" \
  --out_tables "$TABLE_DIR" \
  --out_figures "$FIG_DIR" \
  --nreads "$NREADS" \
  --run_id "$RUN_ID" \
  2>&1 | tee "$LOG_DIR/python_qc.log"

echo "Done."
echo "Outputs:"
echo " - $TABLE_DIR"
echo " - $FIG_DIR"
echo " - $LOG_DIR"
