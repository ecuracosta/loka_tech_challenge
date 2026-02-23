#!/usr/bin/env python3
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse
import matplotlib.pyplot as plt


def read_tsv_gz(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        return pd.read_csv(f, sep="\t", header=None)


def load_10x_mtx(matrix_dir: Path):
    """Load 10x MTX (matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz)."""
    mtx_path = matrix_dir / "matrix.mtx.gz"
    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    features_path = matrix_dir / "features.tsv.gz"

    if not (mtx_path.exists() and barcodes_path.exists() and features_path.exists()):
        raise FileNotFoundError(f"Missing 10x matrix files in {matrix_dir}")

    # matrix is features x barcodes
    with gzip.open(mtx_path, "rb") as f:
        mtx = spio.mmread(f).tocsr()

    barcodes = read_tsv_gz(barcodes_path)[0].astype(str).tolist()

    feats = read_tsv_gz(features_path)
    # 10x features.tsv.gz usually: id, name, feature_type (v3+)
    if feats.shape[1] >= 3:
        feature_id = feats[0].astype(str).tolist()
        feature_name = feats[1].astype(str).tolist()
        feature_type = feats[2].astype(str).tolist()
    else:
        feature_id = feats[0].astype(str).tolist()
        feature_name = feats[1].astype(str).tolist() if feats.shape[1] > 1 else feature_id
        feature_type = ["Unknown"] * len(feature_id)

    features = pd.DataFrame(
        {"feature_id": feature_id, "feature_name": feature_name, "feature_type": feature_type}
    )
    return mtx, barcodes, features


def split_modalities(mtx: sparse.csr_matrix, features: pd.DataFrame):
    """Return indices for GEX and CRISPR-like features."""
    ft = features["feature_type"].str.lower()

    gex_idx = ft.isin(["gene expression", "gex", "rna"])
    # common labels for guide capture
    crispr_idx = ft.str.contains("crispr") | ft.str.contains("guide") | ft.str.contains("feature barcode")

    return np.where(gex_idx.values)[0], np.where(crispr_idx.values)[0]


def guide_calling_simple(guide_counts: sparse.csr_matrix, guide_names: list[str],
                         barcodes: list[str], min_umis: int, ratio_top2: float) -> pd.DataFrame:
    """
    Simple calling:
      - compute top1, top2 counts per cell
      - call if top1 >= min_umis and (top2==0 or top1/top2 >= ratio_top2)
    """
    # guide_counts: guides x cells
    gc = guide_counts.tocsr()

    top1_name = []
    top1 = []
    top2 = []
    n_nonzero = []
    call = []
    status = []

    for j, bc in enumerate(barcodes):
        col = gc[:, j]
        if col.nnz == 0:
            top1_name.append("")
            top1.append(0)
            top2.append(0)
            n_nonzero.append(0)
            call.append("")
            status.append("no_guide_counts")
            continue

        vals = col.data
        idxs = col.indices  # guide indices

        # sort descending
        order = np.argsort(vals)[::-1]
        vals_sorted = vals[order]
        idxs_sorted = idxs[order]

        t1 = int(vals_sorted[0])
        g1 = guide_names[int(idxs_sorted[0])]
        t2 = int(vals_sorted[1]) if len(vals_sorted) > 1 else 0

        top1_name.append(g1)
        top1.append(t1)
        top2.append(t2)
        n_nonzero.append(int(len(vals_sorted)))

        if t1 < min_umis:
            call.append("")
            status.append("below_min_umis")
        else:
            if t2 == 0:
                call.append(g1)
                status.append("called_single")
            else:
                if (t1 / max(t2, 1)) >= ratio_top2:
                    call.append(g1)
                    status.append("called_ratio_ok")
                else:
                    call.append("")
                    status.append("ambiguous_multi")

    df = pd.DataFrame({
        "barcode": barcodes,
        "top1_guide": top1_name,
        "top1_umis": top1,
        "top2_umis": top2,
        "n_guides_nonzero": n_nonzero,
        "called_guide": call,
        "status": status,
    })
    return df


def plot_hist(values, title, outpath):
    plt.figure()
    plt.hist(values, bins=60)
    plt.title(title)
    plt.xlabel("Value")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def plot_top_guides(guide_sums: pd.Series, outpath: Path, topn=30):
    top = guide_sums.sort_values(ascending=False).head(topn)
    plt.figure(figsize=(10, 5))
    plt.bar(top.index.astype(str), top.values)
    plt.xticks(rotation=75, ha="right")
    plt.title(f"Top {topn} guides by total UMIs")
    plt.xlabel("Guide")
    plt.ylabel("UMIs")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--cellranger-outs", required=True, help="Path to cellranger multi outs/ directory")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-guide-umis", type=int, default=3)
    ap.add_argument("--ratio-top2", type=float, default=3.0)
    args = ap.parse_args()

    run_id = args.run_id
    outs = Path(args.cellranger_outs)
    outdir = Path(args.outdir)

    # Prefer filtered matrix for guide calling; also load raw for QC if needed later
    filtered = outs / "filtered_feature_bc_matrix"
    raw = outs / "raw_feature_bc_matrix"

    if not filtered.exists():
        # fallback: maybe you only have the tarball matrix you downloaded, but after cellranger it should exist
        raise FileNotFoundError(f"Expected {filtered} to exist. Did cellranger multi finish?")

    mtx_f, bcs, feats = load_10x_mtx(filtered)
    gex_idx, crispr_idx = split_modalities(mtx_f, feats)

    qc_tables = outdir / "qc" / "tables"
    qc_figs = outdir / "qc" / "figures"
    sec_dir = outdir / "secondary"
    qc_tables.mkdir(parents=True, exist_ok=True)
    qc_figs.mkdir(parents=True, exist_ok=True)
    sec_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "run_id": run_id,
        "cellranger_outs": str(outs),
        "n_barcodes_filtered": len(bcs),
        "n_features_total_filtered": int(mtx_f.shape[0]),
        "n_gex_features_filtered": int(len(gex_idx)),
        "n_crispr_features_filtered": int(len(crispr_idx)),
    }

    # If no crispr features found, bail with an informative message
    if len(crispr_idx) == 0:
        summary["error"] = "No CRISPR/Feature Barcode features detected in filtered_feature_bc_matrix."
        (qc_tables / f"{run_id}__qc_summary.json").write_text(json.dumps(summary, indent=2))
        raise RuntimeError(summary["error"])

    # Extract guide matrix (guides x cells)
    guide_mtx = mtx_f[crispr_idx, :]
    guide_feats = feats.iloc[crispr_idx].copy()
    guide_names = guide_feats["feature_name"].astype(str).tolist()

    # Guide totals + per-cell totals
    guide_totals = np.asarray(guide_mtx.sum(axis=1)).ravel()
    cell_guide_umis = np.asarray(guide_mtx.sum(axis=0)).ravel()

    guide_sums = pd.Series(guide_totals, index=guide_names)
    guide_sums.to_csv(qc_tables / f"{run_id}__guide_totals.tsv", sep="\t", header=False)

    # Save per-cell guide counts summary
    pd.DataFrame({"barcode": bcs, "guide_umis": cell_guide_umis}).to_csv(
        qc_tables / f"{run_id}__guide_counts_per_cell.tsv", sep="\t", index=False
    )

    # Calling
    calls = guide_calling_simple(
        guide_counts=guide_mtx,
        guide_names=guide_names,
        barcodes=bcs,
        min_umis=args.min_guide_umis,
        ratio_top2=args.ratio_top2,
    )
    calls.to_csv(qc_tables / f"{run_id}__guide_calls.tsv", sep="\t", index=False)

    # QC summary stats
    summary.update({
        "guide_calling": {
            "min_guide_umis": args.min_guide_umis,
            "ratio_top2": args.ratio_top2,
            "n_called": int((calls["called_guide"] != "").sum()),
            "n_ambiguous_multi": int((calls["status"] == "ambiguous_multi").sum()),
            "n_no_guide_counts": int((calls["status"] == "no_guide_counts").sum()),
        }
    })

    # Plots
    plot_hist(
        cell_guide_umis,
        title=f"{run_id}: guide UMIs per cell (filtered)",
        outpath=qc_figs / f"{run_id}__guide_umis_per_cell_hist.png",
    )
    plot_top_guides(
        guide_sums=guide_sums,
        outpath=qc_figs / f"{run_id}__top_guides_by_umis.png",
        topn=30,
    )

    (qc_tables / f"{run_id}__qc_summary.json").write_text(json.dumps(summary, indent=2))

    print("✅ Done")
    print(f"- tables:  {qc_tables}")
    print(f"- figures: {qc_figs}")


if __name__ == "__main__":
    main()
