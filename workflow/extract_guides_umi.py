#!/usr/bin/env python3
import argparse
import gzip
import json
from pathlib import Path
from collections import defaultdict

def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return 10**9
    return sum(x != y for x, y in zip(a, b))

def parse_range(s: str):
    # "0,16" -> (0,16)
    a, b = s.split(",")
    return int(a), int(b)

def open_maybe_gz(path: str):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

def load_cb_whitelist(cb_whitelist_gz: Path):
    # barcodes.tsv.gz has "AAAAAAAAAAAAAAAA-1"
    cbs = []
    with gzip.open(cb_whitelist_gz, "rt") as f:
        for line in f:
            cb = line.strip().split("-")[0]
            if cb:
                cbs.append(cb)
    return cbs

def build_seed_index(cbs, seed_len=8):
    idx = defaultdict(list)
    for cb in cbs:
        idx[cb[:seed_len]].append(cb)
    return idx

def correct_cb(cb_obs: str, cb_set: set, seed_index, max_mm: int, seed_len=8):
    if cb_obs in cb_set:
        return cb_obs, 0, True

    if max_mm <= 0:
        return None, None, False

    # candidates from seed bucket
    bucket = seed_index.get(cb_obs[:seed_len], [])
    if not bucket:
        return None, None, False

    best = None
    best_d = 10**9
    tie = False
    for cb in bucket:
        d = hamming(cb_obs, cb)
        if d < best_d:
            best, best_d, tie = cb, d, False
        elif d == best_d:
            tie = True

    if best_d <= max_mm and (not tie):
        return best, best_d, True
    return None, None, False

def load_feature_ref(tsv: Path):
    # each line: name \t sequence
    feats = []
    with open(tsv, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, seq = parts[0], parts[1].upper()
            feats.append((name, seq))
    return feats

def best_guide_match(obs: str, feats, max_mm: int):
    best = None
    best_d = 10**9
    tie = False
    for name, seq in feats:
        if len(obs) != len(seq):
            continue
        d = hamming(obs, seq)
        if d < best_d:
            best, best_d, tie = name, d, False
        elif d == best_d:
            tie = True
    if best is None:
        return None, None, False
    if best_d <= max_mm and (not tie):
        return best, best_d, True
    return None, None, False

def iter_fastq_pairs(r1_path: Path, r2_path: Path):
    r1_fh = gzip.open(r1_path, "rt") if str(r1_path).endswith(".gz") else open(r1_path, "rt")
    r2_fh = gzip.open(r2_path, "rt") if str(r2_path).endswith(".gz") else open(r2_path, "rt")
    try:
        while True:
            h1 = r1_fh.readline()
            h2 = r2_fh.readline()
            if not h1 or not h2:
                break
            s1 = r1_fh.readline().strip()
            s2 = r2_fh.readline().strip()
            r1_fh.readline(); r2_fh.readline()  # +
            r1_fh.readline(); r2_fh.readline()  # qual
            yield s1, s2
    finally:
        r1_fh.close()
        r2_fh.close()

def write_mtx(outdir: Path, counts, barcodes, features):
    # counts: dict[(feat_idx, cell_idx)] -> int
    # MatrixMarket coordinate integer general
    nnz = len(counts)
    with open(outdir/"matrix.mtx", "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{len(features)} {len(barcodes)} {nnz}\n")
        for (i, j), v in sorted(counts.items()):
            f.write(f"{i+1} {j+1} {v}\n")

    with open(outdir/"barcodes.tsv", "wt") as f:
        for cb in barcodes:
            f.write(cb + "\n")

    with open(outdir/"features.tsv", "wt") as f:
        for feat in features:
            f.write(f"{feat}\t{feat}\tCRISPR\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-1", dest="r1", required=True, type=Path)
    ap.add_argument("-2", dest="r2", required=True, type=Path)
    ap.add_argument("-w", dest="cb_whitelist", required=True, type=Path)   # barcodes.tsv.gz
    ap.add_argument("-f", dest="feature_ref", required=True, type=Path)    # edited.tsv
    ap.add_argument("-o", dest="outdir", required=True, type=Path)

    ap.add_argument("--max-cb-mm", type=int, default=2)
    ap.add_argument("--max-guide-mm", type=int, default=2)

    ap.add_argument("--r1-cb", type=str, default="0,16")
    ap.add_argument("--umi-start", type=int, default=16)
    ap.add_argument("--umi-len", type=int, default=12)
    ap.add_argument("--r2-guide", type=str, default="31,51")

    ap.add_argument("--seed-len", type=int, default=8)

    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cb0, cb1 = parse_range(args.r1_cb)
    g0, g1 = parse_range(args.r2_guide)

    whitelist = load_cb_whitelist(args.cb_whitelist)
    cb_set = set(whitelist)
    seed_index = build_seed_index(whitelist, seed_len=args.seed_len)

    feats = load_feature_ref(args.feature_ref)
    feat_names = [n for n, _ in feats]

    total = 0
    pass_cb = 0
    pass_guide = 0

    # Dedup UMI per (cell, guide, umi)
    seen = set()

    # Count unique UMIs per (cell, guide)
    per_cell_guide = defaultdict(int)

    for r1_seq, r2_seq in iter_fastq_pairs(args.r1, args.r2):
        total += 1

        cb_obs = r1_seq[cb0:cb1].upper()
        umi = r1_seq[args.umi_start:args.umi_start+args.umi_len].upper()

        cb_corr, cb_mm, ok_cb = correct_cb(cb_obs, cb_set, seed_index, args.max_cb_mm, seed_len=args.seed_len)
        if not ok_cb:
            continue
        pass_cb += 1

        guide_obs = r2_seq[g0:g1].upper()
        guide, gmm, ok_g = best_guide_match(guide_obs, feats, args.max_guide_mm)
        if not ok_g:
            continue
        pass_guide += 1

        key = (cb_corr, guide, umi)
        if key in seen:
            continue
        seen.add(key)
        per_cell_guide[(cb_corr, guide)] += 1

    # Build sparse counts matrix (features x cells)
    # Keep only cells that appear
    cells = sorted({cb for (cb, g) in per_cell_guide.keys()})
    cell_to_idx = {cb:i for i, cb in enumerate(cells)}
    feat_to_idx = {f:i for i, f in enumerate(feat_names)}

    counts = {}
    for (cb, g), v in per_cell_guide.items():
        i = feat_to_idx[g]
        j = cell_to_idx[cb]
        counts[(i, j)] = v

    write_mtx(args.outdir, counts, cells, feat_names)

    summary = {
        "total_read_pairs": total,
        "whitelist_passing_pairs": pass_cb,
        "guide_matched_pairs": pass_guide,
        "cells_with_ge1_guide_umi": len(cells),
        "guides": len(feat_names),
        "max_cb_mm": args.max_cb_mm,
        "max_guide_mm": args.max_guide_mm,
        "r1_cb": args.r1_cb,
        "r2_guide": args.r2_guide,
        "umi_start": args.umi_start,
        "umi_len": args.umi_len,
        "seed_len": args.seed_len,
    }
    (args.outdir/"summary.json").write_text(json.dumps(summary, indent=2))

    print("Done.")
    print(f"Total read pairs processed: {total:,}")
    print(f"Whitelist-passing pairs:   {pass_cb:,}")
    print(f"Guide-matched pairs:       {pass_guide:,}")
    print(f"Cells with ≥1 guide UMI:   {len(cells):,}")
    print(f"Guides:                    {len(feat_names)}")
    print(f"Wrote: {args.outdir}")

if __name__ == "__main__":
    main()
    