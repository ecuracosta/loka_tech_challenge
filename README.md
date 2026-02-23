# Perturb-seq Pipeline (Loka Tech Challenge)

This repository contains an end-to-end **Perturb-seq**-style workflow that emulates an AWS EventBridge-triggered pipeline locally, runs **primary processing** (GEX quantification + CRISPR guide extraction), and produces **secondary analysis-ready artifacts** (merged RNA+guides `AnnData` + summary plots/metrics).

The goal of this challenge deliverable is to demonstrate:
- A production-minded pipeline structure (clear stages, reproducible runs, debuggable interfaces)
- Proper separation of concerns: workflow orchestration (**Nextflow**) vs analysis utilities (**Python**)
- Practical handling of real-world issues (permissions/caches in containers, whitelists, lane handling, robust debug paths)
- A path to cloud execution (Docker + AWS cost estimate)

![Pipeline Architecture](architecture.png)

---

## High-level workflow

**Trigger ➜ Prepare inputs ➜ Run Nextflow (STARsolo) ➜ Extract guides (UMI-deduped) ➜ Integrate RNA+guides ➜ Secondary analysis outputs**

### Primary outputs
- **RNA (GEX)**: STARsolo matrices + BAMs (per sample)
- **Guides (CRISPR)**: guide count matrix (cells × guides), optionally **UMI-deduplicated**
- **Merged object**: `merged_rna_guides.h5ad` for downstream Perturb-seq analysis

### Secondary outputs (example)
- guide assignment per cell (`NEG`, `SINGLE`, `MULTI`)
- guide abundance plots / embeddings (lightweight summary)
- summary CSVs for quick inspection

---

## What I built / key engineering choices

### Orchestration: Nextflow (DSL2)
- A **Nextflow** workflow coordinates heavy compute steps (e.g. **STARsolo**).
- Parameters are exposed via CLI for explicit, testable runs (samplesheet path, STAR index, whitelists, threads, etc.).
- Supports incremental reruns using `-resume`.

### Trigger + glue layer: Python (`run_pipeline.py`)
- A Python entrypoint emulates an “AWS-like trigger” locally:
  - scans for input FASTQs
  - optionally down-samples FASTQs for ultra-fast debugging
  - writes a Nextflow samplesheet
  - launches Nextflow with explicit parameters
- Includes a **debug-fast mode** to minimize iteration time (single sample, low fraction downsample, low threads).

### Guide extraction: UMI-aware (Perturb-seq requirement)
- Implemented a guide extractor that counts **unique UMIs** per `(cell_barcode, guide, UMI)`:
  - avoids inflated guide counts due to PCR duplicates
  - produces a sparse MatrixMarket matrix compatible with Scanpy/AnnData tools
- Uses:
  - **cell barcode whitelist** (10x filtered barcodes)
  - **feature reference** (guide sequences)
  - configurable mismatch tolerances (CB and guide matching)

### Integration: RNA + guides
- A Python integration script merges:
  - STARsolo filtered RNA matrix (cells × genes)
  - guide matrix (cells × guides)
- Outputs a single `AnnData` `.h5ad` file ready for secondary Perturb-seq analysis.

### Reproducibility: Docker
- A Docker image was built to ensure consistent versions across environments:
  - STAR / STARsolo
  - Python scientific stack (numpy/scipy/pandas/anndata)
  - optional Scanpy-based plotting utilities
- Container runtime issues (Numba/Matplotlib cache permissions) are handled via environment variables:
  - `NUMBA_CACHE_DIR`
  - `MPLCONFIGDIR`

### AWS costs
- A cloud cost estimate is provided in **`cost_estimates.pdf`**:
  - includes a pragmatic view of compute + storage costs for running this pipeline on AWS
  - sized for typical runs and realistic service choices (Batch/ECS, S3, ECR, logs)

---

## Repository structure (important files)

> Paths may vary slightly depending on the environment; these are the logical components.

### Workflow (Nextflow)
- `workflow/main.nf`  
  Orchestrates primary processing steps (e.g. STARsolo).  
- `workflow/nextflow.config`  
  Local execution defaults (cpus/memory/time/workDir), configurable via CLI.

### Primary processing scripts
- `scripts/star_solo.sh`  
  Runs STARsolo for GEX quantification. Produces:
  - `Solo.out/Gene/{filtered,raw}/{matrix.mtx,barcodes.tsv,features.tsv}`
  - sorted BAM + STAR logs

### Trigger / pipeline entrypoint
- `workflow/run_pipeline.py`  
  Local trigger emulator + samplesheet generator + Nextflow launcher.  
  Key features:
  - `--debug-fast`: ultra-fast iteration (single pair + downsample)
  - `--downsample --fraction`: lightweight local debugging without full compute costs

### Guide extraction + integration
- `workflow/extract_guides_umi.py`  
  Extracts guides from CRISPR FASTQs and counts **unique UMIs** per cell/guide.
- `workflow/integrate_rna_guides.py`  
  Merges RNA + guides into a single `AnnData` `.h5ad` and assigns guides per cell.

### Secondary / plotting utilities
- `workflow/plot_guides_scanpy.py`  
  Lightweight guide-matrix visualization (guarded for tiny/degenerate matrices).

### AWS estimate
- `cost_estimates.pdf`  
  Estimated AWS costs for running the pipeline in cloud.

---

# Key Results Produced

---

## 🧬 STARsolo (GEX)

### RNA Count Matrices

- `Solo.out/Gene/filtered/matrix.mtx`
- `Solo.out/Gene/filtered/barcodes.tsv`
- `Solo.out/Gene/filtered/features.tsv`

---

## 🧪 Guides (CRISPR)

### Guide Count Matrix (UMI-deduplicated)

- `/work/outputs/guide_matrix_umi/matrix.mtx`
- `/work/outputs/guide_matrix_umi/barcodes.tsv`
- `/work/outputs/guide_matrix_umi/features.tsv`

---

## 🔗 Integrated Object

- `/work/outputs/merged_rna_guides.h5ad`

### Includes

- Guide calls as cell metadata:
  - `NEG` (no guide detected)
  - `single-guide`
  - `MULTI` (multiple guides)

---

# Notes on Assumptions / Open Extensions

---

## 🎯 Scope of This Prototype

- Focused on:
  - Producing correct **primary artifacts**
  - Providing a clean bridge to **secondary analysis**

---

## ☁️ AWS Deployment Mapping

In a full cloud deployment, the local trigger logic maps to:

- **S3 event**
  ➜ **EventBridge**
  ➜ **Step Functions**
  ➜ **Batch / ECS** *(or Nextflow Tower)*

---

## 🛠️ Production Hardening (Future Work)

- **Structured logging**
- **Richer QC report artifacts**
- **Formal unit tests** for parsing and matching logic
- **Per-stage provenance metadata**
  - Tool versions
  - Checksums