#!/usr/bin/env bash
set -euo pipefail

PROFILE="${1:-local}"
RUN_ID="${2:-SC3_v3_NextGem_DI_CRISPR_10K}"
DO_INGEST="${3:-false}"
DO_QC="${4:-true}"

nextflow run workflow/main.nf \
  -profile "$PROFILE" \
  --run_id "$RUN_ID" \
  --do_ingest "$DO_INGEST" \
  --do_qc "$DO_QC" \
  -with-report outputs/nf_report_${RUN_ID}.html \
  -with-timeline outputs/nf_timeline_${RUN_ID}.html \
  -with-trace outputs/nf_trace_${RUN_ID}.tsv
