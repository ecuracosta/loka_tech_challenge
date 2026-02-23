#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import scipy.sparse as sp
import scanpy as sc
import os
os.environ["NUMBA_DISABLE_CACHING"] = "1"
import matplotlib.pyplot as plt

def read_mtx(mtx_path: Path, n_feat: int, n_cell: int):
    # MatrixMarket: rows=features, cols=cells, 1-based
    rows, cols, data = [], [], []
    with open(mtx_path, "rt") as f:
        header = f.readline()
        dims = f.readline().strip().split()
        # if file has comments, you'd need to skip; here assume 2-line header like we wrote
        for line in f:
            i, j, v = line.strip().split()
            rows.append(int(i)-1)
            cols.append(int(j)-1)
            data.append(int(v))
    X = sp.coo_matrix((data, (rows, cols)), shape=(n_feat, n_cell)).tocsr()
    return X

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, type=Path, help="Directory with matrix.mtx, barcodes.tsv, features.tsv")
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    barcodes = [l.strip() for l in open(args.indir/"barcodes.tsv")]
    feats = [l.strip().split("\t")[0] for l in open(args.indir/"features.tsv")]
    X = read_mtx(args.indir/"matrix.mtx", n_feat=len(feats), n_cell=len(barcodes))

    # scanpy expects cells × features, so transpose
    ad = sc.AnnData(X=X.T)
    ad.obs_names = barcodes
    ad.var_names = feats

    # Normalize + log
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    # --- Robust guard for tiny / degenerate matrices ---
    X = ad.X
    try:
        var = np.array(X.toarray()).var(axis=0) if hasattr(X, "toarray") else np.array(X).var(axis=0)
        ok = np.isfinite(var).all() and (var > 0).sum() >= 2 and ad.n_obs >= 10
    except Exception:
        ok = False

    if ok:
        sc.pp.pca(ad, n_comps=min(10, ad.n_obs - 1, ad.n_vars - 1))
        sc.pp.neighbors(ad)
        sc.tl.umap(ad)
    else:
        print(f"[WARN] Skipping PCA/UMAP (n_obs={ad.n_obs}, n_vars={ad.n_vars}) — matrix too small/degenerate.")

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Heatmap: abundancia por guía
    # (simple: mean log1p counts por feature)
    means = np.asarray(ad.X.mean(axis=0)).ravel()
    df = pd.DataFrame({"feature": feats, "mean_log1p": means}).sort_values("mean_log1p", ascending=False)
    df.to_csv(args.outdir/"feature_mean_log1p.csv", index=False)

    # UMAP colored by the top feature (hard assignment)
    # (si querés más parecido al tutorial, luego hacemos demultiplex thresholds)
    Xdense = ad.X.toarray() if sp.issparse(ad.X) else ad.X
    top_idx = Xdense.argmax(axis=1)
    ad.obs["top_feature"] = [feats[i] for i in top_idx]

    # 1) guides per cell (cuántas guías con count > 0 por célula)
    guides_per_cell = (Xdense > 0).sum(axis=1)
    plt.figure()
    plt.hist(guides_per_cell, bins=np.arange(guides_per_cell.max() + 2) - 0.5)
    plt.xlabel("Guides detected per cell")
    plt.ylabel("Number of cells")
    plt.tight_layout()
    plt.savefig(args.outdir / "guides_per_cell_hist.png", dpi=200)
    plt.close()

    # 2) scatter guide1 vs guide2 (solo si hay exactamente 2 guías)
    if Xdense.shape[1] == 2:
        plt.figure()
        plt.scatter(Xdense[:, 0], Xdense[:, 1])
        plt.xlabel(feats[0])
        plt.ylabel(feats[1])
        plt.tight_layout()
        plt.savefig(args.outdir / "guide_scatter_2guides.png", dpi=200)
        plt.close()

    if "X_umap" in ad.obsm:
        sc.pl.umap(ad, color="top_feature", show=False)
        plt.savefig(args.outdir / "umap_top_feature.png", dpi=200, bbox_inches="tight")
        plt.close()
    else:
        print("[WARN] No UMAP embedding found; skipping umap_top_feature.png")

    print("Wrote:", args.outdir)

if __name__ == "__main__":
    main()
