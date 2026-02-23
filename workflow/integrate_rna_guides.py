#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
import anndata as ad
ad.settings.allow_write_nullable_strings = True

def read_10x_mtx(folder: Path):
    """
    Reads MatrixMarket + barcodes.tsv + features.tsv (10x-style).
    Returns (X, barcodes, features).
    Assumes matrix.mtx has shape (features x cells).
    """
    mtx = mmread(folder / "matrix.mtx").tocsr()
    barcodes = pd.read_csv(folder / "barcodes.tsv", header=None, sep="\t")[0].astype(str).tolist()
    feats_df = pd.read_csv(folder / "features.tsv", header=None, sep="\t")
    # In STARsolo Solo.out/Gene, features.tsv is usually 2 columns (id, name) or 3 cols.
    # We'll keep first col as var_names and second as gene_symbols if present.
    feat_ids = feats_df.iloc[:, 0].astype(str).tolist()
    feat_names = feats_df.iloc[:, 1].astype(str).tolist() if feats_df.shape[1] > 1 else feat_ids
    return mtx, barcodes, feat_ids, feat_names

def strip_suffix(bc: str) -> str:
    # STARsolo barcodes often look like "AAAC...-1"
    return bc.split("-")[0]

def clr_normalize(X: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    # centered log-ratio per cell (row)
    X = X.astype(float)
    gm = np.exp(np.mean(np.log(X + eps), axis=1, keepdims=True))
    return np.log1p(X / gm)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", required=True, type=Path, help="STARsolo Solo.out/Gene/filtered (or raw) folder")
    ap.add_argument("--guides", required=True, type=Path, help="Guide matrix folder (matrix.mtx, barcodes.tsv, features.tsv)")
    ap.add_argument("--out", required=True, type=Path, help="Output .h5ad path")
    ap.add_argument("--min-guide-count", type=int, default=3, help="Min counts to call a guide positive (after merge)")
    ap.add_argument("--max-multi-fraction", type=float, default=0.1, help="If 2nd best guide >= frac*best -> mark as MULTI")
    args = ap.parse_args()

    # --- Load RNA ---
    rna = Path(args.rna)
    candidate = rna / "filtered"
    if not candidate.exists():
        candidate = rna / "raw"
    rna_dir = candidate
    rna_mtx, rna_bcs, rna_feat_ids, rna_feat_names = read_10x_mtx(args.rna)
    # RNA matrix is (genes x cells) -> AnnData expects (cells x genes)
    rna_adata = ad.AnnData(X=rna_mtx.T)
    rna_adata.obs_names = [strip_suffix(b) for b in rna_bcs]
    rna_adata.var_names = rna_feat_ids
    rna_adata.var["gene_symbol"] = rna_feat_names

    # --- Load guides ---
    g_mtx, g_bcs, g_feat_ids, g_feat_names = read_10x_mtx(args.guides)
    # guide matrix is (guides x cells) -> transpose to (cells x guides)
    gX = g_mtx.T.tocsr()
    g_bcs_clean = [strip_suffix(b) for b in g_bcs]  # yours are already 16bp, harmless
    guide_names = g_feat_ids  # like "RAB1A-2", "NON_TARGET-1"

    # --- Intersect cells ---
    rna_cells = pd.Index(rna_adata.obs_names)
    guide_cells = pd.Index(g_bcs_clean)
    common = rna_cells.intersection(guide_cells)

    if len(common) == 0:
        raise SystemExit("No overlapping cell barcodes between RNA and guide matrices. Check barcode formatting.")

    # subset in matching order
    rna_adata = rna_adata[common].copy()

    # build guide matrix aligned to common cells
    guide_pos = {bc: i for i, bc in enumerate(g_bcs_clean)}
    idx = np.array([guide_pos[bc] for bc in common], dtype=int)
    gX_common = gX[idx, :]

    # store guides in obsm (dense small matrix is ok here)
    g_dense = gX_common.toarray()
    rna_adata.obsm["guides_counts"] = g_dense
    rna_adata.uns["guide_names"] = guide_names

    # --- Demultiplex / assign guide per cell (simple but effective) ---
    # Strategy:
    # 1) choose top guide by raw counts
    # 2) call NEG if top < min-guide-count
    # 3) call MULTI if 2nd >= max_multi_fraction*top
    top = g_dense.argmax(axis=1)
    top_counts = g_dense[np.arange(g_dense.shape[0]), top]
    # second best
    sorted_counts = np.sort(g_dense, axis=1)
    second = sorted_counts[:, -2] if g_dense.shape[1] >= 2 else np.zeros_like(top_counts)

    labels = []
    for i in range(g_dense.shape[0]):
        if top_counts[i] < args.min_guide_count:
            labels.append("NEG")
        elif second[i] >= args.max_multi_fraction * top_counts[i]:
            labels.append("MULTI")
        else:
            labels.append(guide_names[top[i]])

    rna_adata.obs["guide_call"] = pd.Categorical(labels)

    # Also store CLR-normalized guide abundances (useful for plotting like the tutorial)
    rna_adata.obsm["guides_clr"] = clr_normalize(g_dense)

    # Save
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rna_adata.write_h5ad(args.out)
    print(f"Wrote merged AnnData: {args.out}")
    print(rna_adata.obs["guide_call"].value_counts())

if __name__ == "__main__":
    main()
