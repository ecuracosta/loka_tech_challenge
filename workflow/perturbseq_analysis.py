#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import scipy.sparse as sp

sc.settings.n_jobs = 1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--outdir", required=True)
    # opcional: te deja ajustar sin hardcode
    p.add_argument("--target-sum", type=float, default=1e4)
    return p.parse_args()


def normalize_total_safe(adata, target_sum=1e4):
    """
    Robust per-cell library size normalization:
      X <- X / libsize * target_sum
    Works for sparse CSR/CSC and dense.
    Avoids scanpy.pp.normalize_total internal edge-cases.
    """
    X = adata.X

    if sp.issparse(X):
        X = X.tocsr(copy=False)
        lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        # prevent division by zero
        lib[lib == 0] = 1.0
        scale = (target_sum / lib).astype(np.float64)
        # CSR row scaling (safe, no densify)
        X = sp.diags(scale).dot(X)
        adata.X = X
    else:
        lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        lib[lib == 0] = 1.0
        adata.X = (X / lib[:, None]) * float(target_sum)


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    h5ad_path = Path(args.h5ad)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"--h5ad not found: {h5ad_path}")

    try:
        adata = sc.read_h5ad(h5ad_path)
    except Exception:
        print(f"[ERROR] Failed to read h5ad: {h5ad_path}", file=sys.stderr)
        raise

    # --- robust: drop zero-count cells (prevents downstream crashes/segfaults) ---
    X = adata.X
    row_sum = np.asarray(X.sum(axis=1)).ravel() if sp.issparse(X) else np.asarray(X.sum(axis=1)).ravel()
    keep = row_sum > 0
    if (~keep).any():
        print(f"[WARN] Dropping {(~keep).sum()} cells with zero counts", file=sys.stderr)
        adata = adata[keep].copy()

    # if fast-mode leaves nothing, exit gracefully (not a crash)
    if adata.n_obs == 0 or adata.n_vars == 0:
        (outdir / "qc_metrics.csv").write_text("Empty AnnData after filtering.\n")
        (outdir / "umap_by_guide.png").write_text("Empty AnnData after filtering.\n")
        print("[WARN] Empty AnnData after filtering; wrote placeholders.", file=sys.stderr)
        return

    # Basic QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # write qc even if later steps fail
    qc_cols = [c for c in ["total_counts", "n_genes_by_counts", "pct_counts_in_top_50_genes"] if c in adata.obs.columns]
    qc = adata.obs[qc_cols].copy() if qc_cols else pd.DataFrame(index=adata.obs_names)
    qc["guide_call"] = adata.obs.get("guide_call", "NA").astype(str)
    qc.to_csv(outdir / "qc_metrics.csv", index=True)

    # RNA normalization / HVG / embedding
    # PATCH: avoid scanpy.pp.normalize_total bug
    normalize_total_safe(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)

    # HVG can fail if too few genes/cells; guard it
    n_top = min(2000, adata.n_vars)
    if n_top < 2 or adata.n_obs < 3:
        (outdir / "umap_by_guide.png").write_text("Not enough cells/genes for HVG/PCA/UMAP in fast mode.\n")
        print("[WARN] Too few cells/genes for UMAP; QC written.", file=sys.stderr)
        return

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat_v3")
    if "highly_variable" not in adata.var.columns or adata.var["highly_variable"].sum() < 2:
        (outdir / "umap_by_guide.png").write_text("No HVGs found in fast mode.\n")
        print("[WARN] No HVGs found; QC written.", file=sys.stderr)
        return

    ad_hvg = adata[:, adata.var["highly_variable"]].copy()

    n_comps = min(30, ad_hvg.n_obs - 1, ad_hvg.n_vars - 1)
    if n_comps < 2:
        (outdir / "umap_by_guide.png").write_text("Not enough data for PCA/UMAP in fast mode.\n")
        print("[WARN] Not enough data for PCA; QC written.", file=sys.stderr)
        return

    sc.pp.pca(ad_hvg, n_comps=n_comps)
    sc.pp.neighbors(ad_hvg)
    sc.tl.umap(ad_hvg)

    # UMAP colored by guide_call
    if "guide_call" in adata.obs.columns:
        ad_hvg.obs["guide_call"] = adata.obs["guide_call"].astype(str).values
        sc.pl.umap(ad_hvg, color="guide_call", show=False)
        plt.savefig(outdir / "umap_by_guide.png", dpi=200, bbox_inches="tight")
        plt.close()
    else:
        (outdir / "umap_by_guide.png").write_text("No guide_call column found.\n")

    # Differential expression: RAB1A-2 vs NON_TARGET-1 (if present)
    if "guide_call" in adata.obs.columns:
        groups = set(adata.obs["guide_call"].astype(str).values)
        if ("RAB1A-2" in groups) and ("NON_TARGET-1" in groups):
            ad_de = adata.copy()
            ad_de.obs["guide_call"] = ad_de.obs["guide_call"].astype(str)

            sc.tl.rank_genes_groups(
                ad_de,
                groupby="guide_call",
                groups=["RAB1A-2"],
                reference="NON_TARGET-1",
                method="wilcoxon"
            )

            df = sc.get.rank_genes_groups_df(ad_de, group="RAB1A-2")
            df.to_csv(outdir / "de_RAB1A2_vs_NON_TARGET1.csv", index=False)

            x = df["logfoldchanges"].values
            p = df["pvals_adj"].fillna(1.0).values
            y = -np.log10(np.maximum(p, 1e-300))

            plt.figure()
            plt.scatter(x, y, s=6)
            plt.xlabel("logFC (RAB1A-2 vs NON_TARGET-1)")
            plt.ylabel("-log10(p_adj)")
            plt.title("DE volcano")
            plt.savefig(outdir / "volcano_RAB1A2_vs_NON_TARGET1.png", dpi=200, bbox_inches="tight")
            plt.close()
        else:
            (outdir / "de_RAB1A2_vs_NON_TARGET1.csv").write_text("Missing required guide groups.\n")
            (outdir / "volcano_RAB1A2_vs_NON_TARGET1.png").write_text("Missing required guide groups.\n")

    print("Wrote analysis to:", outdir)


if __name__ == "__main__":
    main()
