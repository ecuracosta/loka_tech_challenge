#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Defaults
# ------------------------------------------------------------
ROOT="${ROOT:-/work}"
DATASET="${DATASET:-SC3_v3_NextGem_DI_CRISPR_10K}"
VERSION="${VERSION:-4.0.0}"

# Where your downloaded/unpacked inputs currently live (local)
FASTQ_BASE_LOCAL="${FASTQ_BASE_LOCAL:-$ROOT/data/fastq/SC3_v3_NextGem_DI_CRISPR_10K_fastqs}"
BARCODE_MATRIX_LOCAL="${BARCODE_MATRIX_LOCAL:-$ROOT/data/barcode/filtered_feature_bc_matrix}"
REF_LOCAL_DIR="${REF_LOCAL_DIR:-$ROOT/data/ref}"

# Run identity
RUN_ID="${RUN_ID:-RUN_$(date +%Y%m%dT%H%M%S)}"

# Stage controls
STAGE="${STAGE:-local}"          # local | s3
S3_PREFIX="${S3_PREFIX:-}"       # e.g. s3://loka-processed/customerX/runs
COMBINE_LANES="${COMBINE_LANES:-1}"  # 1 combine L00x into combined, 0 keep only raw
DRYRUN="${DRYRUN:-0}"

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------
usage() {
  cat <<'USAGE'
Usage:
  ingest_10x_run.sh [options]

Options:
  --run-id RUN_ID
  --dataset NAME
  --version VERSION
  --root PATH
  --stage local|s3
  --s3-prefix s3://bucket/prefix           (required if --stage s3)
  --combine-lanes 0|1
  --fastq-base-local PATH                  (where 10x fastq folders are)
  --barcode-matrix-local PATH              (filtered_feature_bc_matrix folder)
  --ref-local-dir PATH                     (feature_ref.tsv/csv folder)
  --dryrun 0|1
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-id) RUN_ID="$2"; shift 2 ;;
    --dataset) DATASET="$2"; shift 2 ;;
    --version) VERSION="$2"; shift 2 ;;
    --root) ROOT="$2"; shift 2 ;;
    --stage) STAGE="$2"; shift 2 ;;
    --s3-prefix) S3_PREFIX="$2"; shift 2 ;;
    --combine-lanes) COMBINE_LANES="$2"; shift 2 ;;
    --fastq-base-local) FASTQ_BASE_LOCAL="$2"; shift 2 ;;
    --barcode-matrix-local) BARCODE_MATRIX_LOCAL="$2"; shift 2 ;;
    --ref-local-dir) REF_LOCAL_DIR="$2"; shift 2 ;;
    --dryrun) DRYRUN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 2 ;;
  esac
done

if [[ "$STAGE" == "s3" && -z "$S3_PREFIX" ]]; then
  echo "ERROR: --s3-prefix is required when --stage s3"
  exit 2
fi

# ------------------------------------------------------------
# Layout
# ------------------------------------------------------------
RUN_DIR="$ROOT/data/runs/$RUN_ID"
FASTQ_DIR="$RUN_DIR/fastq"
GEX_RAW_DIR="$FASTQ_DIR/gex/raw"
CRISPR_RAW_DIR="$FASTQ_DIR/crispr/raw"
GEX_COMBINED_DIR="$FASTQ_DIR/gex/combined"
CRISPR_COMBINED_DIR="$FASTQ_DIR/crispr/combined"
REF_DIR="$RUN_DIR/ref"
MATRIX_DIR="$RUN_DIR/filtered_feature_bc_matrix"
META_DIR="$RUN_DIR/meta"
LOG_DIR="$RUN_DIR/logs"

MANIFEST_JSON="$RUN_DIR/manifest.json"
MANIFEST_TSV="$RUN_DIR/manifest.tsv"

mkdir -p \
  "$GEX_RAW_DIR" "$CRISPR_RAW_DIR" \
  "$GEX_COMBINED_DIR" "$CRISPR_COMBINED_DIR" \
  "$REF_DIR" "$MATRIX_DIR" "$META_DIR" "$LOG_DIR"

log() { echo "==> $*"; }

# ------------------------------------------------------------
# Helpers: link vs copy (materialize)
# ------------------------------------------------------------
materialize_tree() {
  # materialize_tree SRC_DIR DEST_DIR
  local src="$1"
  local dest="$2"

  if [[ "$STAGE" == "local" ]]; then
    # save disk: symlink files
    find "$src" -maxdepth 1 -type f -name "*.fastq.gz" -print0 | while IFS= read -r -d '' f; do
      ln -sfn "$f" "$dest/$(basename "$f")"
    done
  else
    # stage=s3: DO NOT symlink; copy files so sync to S3 is clean
    find "$src" -maxdepth 1 -type f -name "*.fastq.gz" -print0 | while IFS= read -r -d '' f; do
      cp -f "$f" "$dest/$(basename "$f")"
    done
  fi
}

combine_lanes() {
  # combine_lanes RAW_DIR OUT_DIR PREFIX
  local raw="$1"
  local out="$2"
  local prefix="$3"

  # 10x naming: *_L001_R1_001.fastq.gz etc
  for read in I1 I2 R1 R2; do
    local glob="${prefix}_S1_L00?_${read}_001.fastq.gz"
    local outfn="${prefix}_S1_combined_${read}_001.fastq.gz"

    # shellcheck disable=SC2086
    local files=( "$raw"/$glob )
    if [[ "${files[0]}" == "$raw/$glob" ]]; then
      # no match
      continue
    fi

    log "Combine ${read}: ${outfn}"
    if [[ "$DRYRUN" == "1" ]]; then
      printf "DRYRUN: cat %s > %s\n" "${files[*]}" "$out/$outfn"
    else
      # concatenating .gz is valid; gzip readers will stream all members
      cat "${files[@]}" > "$out/$outfn"
    fi
  done
}

# ------------------------------------------------------------
# 1) Materialize raw fastqs into run structure
# ------------------------------------------------------------
GEX_SRC="$FASTQ_BASE_LOCAL/${DATASET}_gex_fastqs"
CRISPR_SRC="$FASTQ_BASE_LOCAL/${DATASET}_crispr_fastqs"

log "ROOT=$ROOT"
log "RUN_ID=$RUN_ID"
log "STAGE=$STAGE"
log "RUN_DIR=$RUN_DIR"
log "GEX_SRC=$GEX_SRC"
log "CRISPR_SRC=$CRISPR_SRC"

if [[ ! -d "$GEX_SRC" || ! -d "$CRISPR_SRC" ]]; then
  echo "ERROR: expected GEX/CRISPR source dirs not found."
  exit 1
fi

log "Materialize GEX raw..."
materialize_tree "$GEX_SRC" "$GEX_RAW_DIR"

log "Materialize CRISPR raw..."
materialize_tree "$CRISPR_SRC" "$CRISPR_RAW_DIR"

# ------------------------------------------------------------
# 2) Copy ref + barcode matrix
# ------------------------------------------------------------
log "Copy ref files..."
if [[ -d "$REF_LOCAL_DIR" ]]; then
  cp -f "$REF_LOCAL_DIR"/* "$REF_DIR/" 2>/dev/null || true
else
  log "WARN: REF_LOCAL_DIR not found: $REF_LOCAL_DIR"
fi

log "Copy filtered_feature_bc_matrix (optional sanity check)..."
if [[ -d "$BARCODE_MATRIX_LOCAL" ]]; then
  cp -f "$BARCODE_MATRIX_LOCAL/"* "$MATRIX_DIR/" 2>/dev/null || true
else
  log "WARN: BARCODE_MATRIX_LOCAL not found: $BARCODE_MATRIX_LOCAL"
fi

# ------------------------------------------------------------
# 3) Combine lanes (optional)
# ------------------------------------------------------------
if [[ "$COMBINE_LANES" == "1" ]]; then
  log "Combine lanes (GEX)..."
  combine_lanes "$GEX_RAW_DIR" "$GEX_COMBINED_DIR" "${DATASET}_gex"

  log "Combine lanes (CRISPR)..."
  combine_lanes "$CRISPR_RAW_DIR" "$CRISPR_COMBINED_DIR" "${DATASET}_crispr"
else
  log "COMBINE_LANES=0 : skipping lane combination"
fi

# ------------------------------------------------------------
# 4) Write manifest (JSON + TSV)
# ------------------------------------------------------------
log "Write manifest.json + manifest.tsv"

cat > "$MANIFEST_JSON" <<EOF
{
  "run_id": "${RUN_ID}",
  "dataset": "${DATASET}",
  "version": "${VERSION}",
  "created_at": "$(date -Is)",
  "stage": "${STAGE}",
  "paths": {
    "run_dir": "${RUN_DIR}",
    "fastq": {
      "gex_raw": "${GEX_RAW_DIR}",
      "crispr_raw": "${CRISPR_RAW_DIR}",
      "gex_combined": "${GEX_COMBINED_DIR}",
      "crispr_combined": "${CRISPR_COMBINED_DIR}"
    },
    "ref": "${REF_DIR}",
    "filtered_feature_bc_matrix": "${MATRIX_DIR}"
  },
  "metadata_sources": {
    "benchling": {
      "status": "TODO",
      "notes": "Integrate Benchling metadata here (sample sheet, constructs, guide library, run info)"
    },
    "smartsheet": {
      "status": "TODO",
      "notes": "Integrate SmartSheet metadata here (run tracking, QC, pipeline params)"
    }
  }
}
EOF

cat > "$MANIFEST_TSV" <<EOF
run_id	${RUN_ID}
dataset	${DATASET}
version	${VERSION}
created_at	$(date -Is)
gex_raw	${GEX_RAW_DIR}
crispr_raw	${CRISPR_RAW_DIR}
gex_combined	${GEX_COMBINED_DIR}
crispr_combined	${CRISPR_COMBINED_DIR}
ref	${REF_DIR}
filtered_feature_bc_matrix	${MATRIX_DIR}
EOF

# ------------------------------------------------------------
# 5) Upload to S3 (optional), manifest LAST
# ------------------------------------------------------------
if [[ "$STAGE" == "s3" ]]; then
  DEST="${S3_PREFIX%/}/$RUN_ID"
  log "Upload run to S3: $DEST"

  if [[ "$DRYRUN" == "1" ]]; then
    echo "DRYRUN: aws s3 sync \"$RUN_DIR\" \"$DEST\" --exclude manifest.json"
    echo "DRYRUN: aws s3 cp  \"$MANIFEST_JSON\" \"$DEST/manifest.json\""
  else
    # Upload everything except manifest first
    aws s3 sync "$RUN_DIR" "$DEST" --exclude "manifest.json"
    # Upload manifest last (robust trigger point)
    aws s3 cp "$MANIFEST_JSON" "$DEST/manifest.json"
  fi

  log "S3 uploaded. Trigger can watch: $DEST/manifest.json"
fi

log "Done. Tree:"
find "$RUN_DIR" -maxdepth 4 -type f | sed 's#^#  #'
