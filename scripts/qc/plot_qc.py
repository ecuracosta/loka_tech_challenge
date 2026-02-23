#!/usr/bin/env python3
import argparse
import gzip
import os
from collections import Counter, defaultdict

import pandas as pd
import matplotlib.pyplot as plt


def read_fastq_per_base_counts_gz(path: str, max_reads: int) -> pd.DataFrame:
    """
    Stream a gzipped FASTQ and compute per-base nucleotide fractions for the first max_reads reads.
    Returns a DataFrame with columns: pos, A, C, G, T, N, other, total
    """
    counts_by_pos = defaultdict(Counter)
    n_reads = 0

    with gzip.open(path, "rt") as fh:
        while n_reads < max_reads:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                break

            n_reads += 1
            for i, base in enumerate(seq):
                b = base.upper()
                if b in ("A", "C", "G", "T", "N"):
                    counts_by_pos[i][b] += 1
                else:
                    counts_by_pos[i]["other"] += 1

    if n_reads == 0:
        raise RuntimeError(f"No reads found in {path}")

    # Build table
    max_len = max(counts_by_pos.keys())
    rows = []
    for pos in range(max_len + 1):
        c = counts_by_pos[pos]
        total = sum(c.values())
        rows.append(
            {
                "pos": pos,
                "A": c.get("A", 0),
                "C": c.get("C", 0),
                "G": c.get("G", 0),
                "T": c.get("T", 0),
                "N": c.get("N", 0),
                "other": c.get("other", 0),
                "total": total,
            }
        )

    df = pd.DataFrame(rows)
    # Convert to fractions
    for col in ["A", "C", "G", "T", "N", "other"]:
        df[col] = df[col] / df["total"].replace(0, pd.NA)

    return df


def plot_per_base_content(df: pd.DataFrame, title: str, out_png: str) -> None:
    plt.figure()
    plt.plot(df["pos"], df["A"], label="A")
    plt.plot(df["pos"], df["C"], label="C")
    plt.plot(df["pos"], df["G"], label="G")
    plt.plot(df["pos"], df["T"], label="T")
    plt.plot(df["pos"], df["N"], label="N")
    if "other" in df.columns:
        plt.plot(df["pos"], df["other"], label="other")
    plt.xlabel("Position (0-based)")
    plt.ylabel("Fraction")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def parse_matching_startpos(field: str):
    """
    fba qc writes positions like '0:16' or '31:51' (strings).
    Return start position as int, else None.
    """
    if not isinstance(field, str):
        return None
    if field in ("NA", "", "None"):
        return None
    if ":" not in field:
        return None
    try:
        start = int(field.split(":")[0])
        return start
    except Exception:
        return None


def plot_hist(values, title: str, out_png: str, xlabel: str) -> None:
    plt.figure()
    # default bins: auto-ish but stable
    plt.hist(values, bins=50)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--qc_tsv_gz", required=True)
    ap.add_argument("--out_tables", required=True)
    ap.add_argument("--out_figures", required=True)
    ap.add_argument("--nreads", type=int, default=20000)
    ap.add_argument("--run_id", default="qc_run")
    args = ap.parse_args()

    os.makedirs(args.out_tables, exist_ok=True)
    os.makedirs(args.out_figures, exist_ok=True)

    # -------------------------
    # Per-base content (FASTQ-derived)
    # -------------------------
    r1_df = read_fastq_per_base_counts_gz(args.r1, args.nreads)
    r2_df = read_fastq_per_base_counts_gz(args.r2, args.nreads)

    r1_tsv = os.path.join(args.out_tables, "per_base_content_read1.tsv")
    r2_tsv = os.path.join(args.out_tables, "per_base_content_read2.tsv")
    r1_df.to_csv(r1_tsv, sep="\t", index=False)
    r2_df.to_csv(r2_tsv, sep="\t", index=False)

    plot_per_base_content(
        r1_df,
        title=f"Reads per base sequence content (R1) [{args.run_id}]",
        out_png=os.path.join(args.out_figures, "qc_reads_per_base_R1.png"),
    )
    plot_per_base_content(
        r2_df,
        title=f"Reads per base sequence content (R2) [{args.run_id}]",
        out_png=os.path.join(args.out_figures, "qc_reads_per_base_R2.png"),
    )

    # -------------------------
    # Feature barcode positions (TSV-derived)
    # -------------------------
    qc = pd.read_csv(args.qc_tsv_gz, sep="\t", compression="gzip")

    # Columns per tutorial example:
    # cb_matching_pos, fb_matching_pos, plus barcodes
    cb_start = qc.get("cb_matching_pos", pd.Series([None] * len(qc))).map(parse_matching_startpos)
    fb_start = qc.get("fb_matching_pos", pd.Series([None] * len(qc))).map(parse_matching_startpos)

    cb_vals = [v for v in cb_start.dropna().astype(int).tolist()]
    fb_vals = [v for v in fb_start.dropna().astype(int).tolist()]

    # Summary table
    summary_rows = []
    summary_rows.append({"metric": "rows_total", "value": len(qc)})
    summary_rows.append({"metric": "cb_matches", "value": len(cb_vals)})
    summary_rows.append({"metric": "fb_matches", "value": len(fb_vals)})
    if cb_vals:
        summary_rows += [
            {"metric": "cb_startpos_min", "value": min(cb_vals)},
            {"metric": "cb_startpos_median", "value": float(pd.Series(cb_vals).median())},
            {"metric": "cb_startpos_max", "value": max(cb_vals)},
        ]
    if fb_vals:
        summary_rows += [
            {"metric": "fb_startpos_min", "value": min(fb_vals)},
            {"metric": "fb_startpos_median", "value": float(pd.Series(fb_vals).median())},
            {"metric": "fb_startpos_max", "value": max(fb_vals)},
        ]

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(args.out_tables, "matching_positions_summary.tsv"), sep="\t", index=False)

    # Hist plots
    if cb_vals:
        plot_hist(
            cb_vals,
            title=f"Distribution of cell-barcode match start positions (R1) [{args.run_id}]",
            out_png=os.path.join(args.out_figures, "qc_cb_matching_startpos.png"),
            xlabel="Start position (0-based)",
        )
    if fb_vals:
        plot_hist(
            fb_vals,
            title=f"Distribution of feature-barcode match start positions (R2) [{args.run_id}]",
            out_png=os.path.join(args.out_figures, "qc_fb_matching_startpos.png"),
            xlabel="Start position (0-based)",
        )

    print("OK")
    print(f"Tables: {args.out_tables}")
    print(f"Figures: {args.out_figures}")


if __name__ == "__main__":
    main()
