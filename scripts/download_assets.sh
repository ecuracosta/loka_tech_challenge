#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------
# Cacheable asset downloader (local or S3)
# ----------------------------------------
# Examples:
#   # Local cache only (downloads from internet if missing)
#   ./scripts/download_assets.sh --asset grch38_refdata --stage local
#
#   # Use S3 as shared cache (download from S3 if exists, else fetch+upload)
#   ./scripts/download_assets.sh --asset grch38_refdata --stage s3 \
#       --s3-bucket my-loka-cache --s3-prefix assets/10x
#
# Outputs:
#   data/ref/refdata-gex-GRCh38-2020-A/  (extracted transcriptome)
#   data/cache/downloads/               (tarballs)
#   data/cache/markers/                 (idempotency markers)

ROOT="${ROOT:-/work}"
STAGE="local"                  # local | s3
ASSET="grch38_refdata"         # currently: grch38_refdata
CACHE_DIR="$ROOT/data/cache"
DOWNLOADS_DIR="$CACHE_DIR/downloads"
MARKERS_DIR="$CACHE_DIR/markers"
REF_DIR="$ROOT/data/ref"

S3_BUCKET=""
S3_PREFIX="assets"             # s3://bucket/prefix/<filename>

FORCE="0"                      # 1=redownload+reextract
QUIET="0"

# Stable (unsigned) URL for 10x GRCh38 refdata
REFDATA_URL_DEFAULT="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"

usage() {
  cat <<'USAGE'
Usage:
  download_assets.sh [options]

Options:
  --asset <name>         Asset key (default: grch38_refdata)
  --stage <local|s3>     Cache stage (default: local)
  --root <path>          Project root (default: /work or $ROOT)
  --cache-dir <path>     Cache dir (default: $ROOT/data/cache)
  --ref-dir <path>       Ref output dir (default: $ROOT/data/ref)
  --s3-bucket <name>     Required if --stage s3
  --s3-prefix <prefix>   S3 prefix (default: assets)
  --url <url>            Override download URL
  --force                Redownload and reextract
  --quiet                Less logging

Examples:
  ./scripts/download_assets.sh --asset grch38_refdata --stage local
  ./scripts/download_assets.sh --asset grch38_refdata --stage s3 --s3-bucket my-bucket --s3-prefix assets/10x
USAGE
}

URL_OVERRIDE=""

# -------------------------
# Parse args
# -------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --asset)      ASSET="$2"; shift 2 ;;
    --stage)      STAGE="$2"; shift 2 ;;
    --root)       ROOT="$2"; shift 2 ;;
    --cache-dir)  CACHE_DIR="$2"; shift 2 ;;
    --ref-dir)    REF_DIR="$2"; shift 2 ;;
    --s3-bucket)  S3_BUCKET="$2"; shift 2 ;;
    --s3-prefix)  S3_PREFIX="$2"; shift 2 ;;
    --url)        URL_OVERRIDE="$2"; shift 2 ;;
    --force)      FORCE="1"; shift ;;
    --quiet)      QUIET="1"; shift ;;
    -h|--help)    usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

DOWNLOADS_DIR="$CACHE_DIR/downloads"
MARKERS_DIR="$CACHE_DIR/markers"

mkdir -p "$DOWNLOADS_DIR" "$MARKERS_DIR" "$REF_DIR"

log() { [[ "$QUIET" == "1" ]] || echo "$@"; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing command: $1"; exit 127; }
}

# We accept curl OR wget; pick one.
pick_downloader() {
  if command -v curl >/dev/null 2>&1; then
    echo "curl"
  elif command -v wget >/dev/null 2>&1; then
    echo "wget"
  else
    echo "none"
  fi
}

download_http() {
  local url="$1"
  local out="$2"
  local tmp="${out}.partial"
  local dl
  dl="$(pick_downloader)"
  [[ "$dl" != "none" ]] || { echo "ERROR: need curl or wget to download over HTTP(S)"; exit 127; }

  log "==> HTTP download: $url"
  if [[ "$dl" == "curl" ]]; then
    # -L follow redirects, -f fail on HTTP errors, -C - resume
    curl -fL -C - -o "$tmp" "$url"
  else
    wget -c -O "$tmp" "$url"
  fi
  mv -f "$tmp" "$out"
}

s3_uri() {
  local filename="$1"
  echo "s3://${S3_BUCKET%/}/${S3_PREFIX%/}/${filename}"
}

have_aws() { command -v aws >/dev/null 2>&1; }

download_s3_if_exists() {
  local uri="$1"
  local out="$2"
  if ! have_aws; then
    echo "ERROR: aws cli not found but --stage s3 was requested."
    echo "Install awscli OR run with --stage local."
    exit 127
  fi

  # Check existence without failing the whole script
  set +e
  aws s3 ls "$uri" >/dev/null 2>&1
  local exists=$?
  set -e

  if [[ $exists -eq 0 ]]; then
    log "==> S3 cache hit: $uri"
    aws s3 cp "$uri" "$out"
    return 0
  fi
  return 1
}

upload_s3() {
  local file="$1"
  local uri="$2"
  if ! have_aws; then
    echo "ERROR: aws cli not found but upload to s3 requested."
    exit 127
  fi
  log "==> Upload to S3 cache: $uri"
  aws s3 cp "$file" "$uri"
}

extract_tar_gz_once() {
  local tarball="$1"
  local out_dir="$2"
  local marker="$3"

  if [[ "$FORCE" == "1" ]]; then
    rm -f "$marker"
    rm -rf "$out_dir"
  fi

  if [[ -f "$marker" ]]; then
    log "==> Extract: already done ($(basename "$marker"))"
    return 0
  fi

  mkdir -p "$out_dir"
  log "==> Extracting: $(basename "$tarball") -> $out_dir"
  tar -xzf "$tarball" -C "$out_dir" --strip-components=1

  # Create marker with some provenance
  {
    echo "tarball=$tarball"
    echo "extracted_at=$(date -Is)"
  } > "$marker"
}

# -------------------------
# Asset registry
# -------------------------
asset_url=""
asset_filename=""
asset_outdir=""
asset_marker=""

case "$ASSET" in
  grch38_refdata)
    asset_url="${URL_OVERRIDE:-$REFDATA_URL_DEFAULT}"
    asset_filename="$(basename "${asset_url%%\?*}")" # remove query string if any
    asset_outdir="$REF_DIR/refdata-gex-GRCh38-2020-A"
    asset_marker="$MARKERS_DIR/${ASSET}.extracted.ok"
    ;;
  *)
    echo "ERROR: unknown asset: $ASSET"
    echo "Supported: grch38_refdata"
    exit 1
    ;;
esac

tarball_path="$DOWNLOADS_DIR/$asset_filename"

# -------------------------
# Main flow
# -------------------------
log "==> ROOT=$ROOT"
log "==> ASSET=$ASSET  STAGE=$STAGE"
log "==> tarball=$tarball_path"
log "==> outdir=$asset_outdir"

if [[ "$FORCE" == "1" ]]; then
  rm -f "$tarball_path"
fi

if [[ "$STAGE" == "s3" ]]; then
  [[ -n "$S3_BUCKET" ]] || { echo "ERROR: --s3-bucket is required when --stage s3"; exit 1; }
  uri="$(s3_uri "$asset_filename")"

  if [[ ! -f "$tarball_path" ]]; then
    if ! download_s3_if_exists "$uri" "$tarball_path"; then
      # Cache miss -> download from internet, then upload
      download_http "$asset_url" "$tarball_path"
      upload_s3 "$tarball_path" "$uri"
    fi
  else
    log "==> Local cache hit: $tarball_path"
  fi
else
  # local stage
  if [[ ! -f "$tarball_path" ]]; then
    download_http "$asset_url" "$tarball_path"
  else
    log "==> Local cache hit: $tarball_path"
  fi
fi

# Minimal sanity check
if [[ ! -s "$tarball_path" ]]; then
  echo "ERROR: downloaded file is empty: $tarball_path"
  exit 1
fi

extract_tar_gz_once "$tarball_path" "$asset_outdir" "$asset_marker"

log "==> DONE"
log "Transcriptome ready at: $asset_outdir"
