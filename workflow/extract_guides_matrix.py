#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict
from pathlib import Path

def read_whitelist(path: Path) -> set[str]:
    wl = set()
    with open(path, "rt") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            wl.add(s)
    return wl

def read_feature_ref(path: Path) -> dict[str, str]:
    # TSV: feature_id \t sequence
    feats = {}
    with open(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fid, seq = line.split("\t")
            feats[seq] = fid
    return feats

def open_maybe_gz(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def iter_fastq_seqs(path: Path):
    # yields sequence lines (NR%4==2)
    with open_maybe_gz(path) as fh:
        i = 0
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip("\n")
            _plus = fh.readline()
            _qual = fh.readline()
            i += 1
            yield seq

def write_mtx(out_dir: Path, counts: dict[tuple[int,int], int], n_cells: int, n_feats: int):
    # Matrix Market format (1-based indices)
    out_dir.mkdir(parents=True, exist_ok=True)
    mtx_path = out_dir / "matrix.mtx"
    with open(mtx_path, "wt") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n")
        f.write(f"{n_feats} {n_cells} {len(counts)}\n")
        # convention: features as rows, cells as columns
        for (cell_j, feat_i), v in counts.items():
            f.write(f"{feat_i+1} {cell_j+1} {v}\n")

def main():
    ap = argparse.ArgumentParser(description="Extract 10x cell barcode + guide barcode and build cell×guide matrix.")
    ap.add_argument("--r1", required=True, type=Path, help="R1 FASTQ(.gz) (contains cell barcode at 0:16)")
    ap.add_argument("--r2", required=True, type=Path, help="R2 FASTQ(.gz) (contains feature at 31:51)")
    ap.add_argument("--whitelist", required=True, type=Path, help="One barcode per line, 16bp (cell-associated)")
    ap.add_argument("--feature-ref", required=True, type=Path, help="TSV: feature_id<TAB>sequence")
    ap.add_argument("--outdir", required=True, type=Path, help="Output directory")
    ap.add_argument("--cb-start", type=int, default=0)
    ap.add_argument("--cb-end", type=int, default=16)
    ap.add_argument("--fb-start", type=int, default=31)
    ap.add_argument("--fb-end", type=int, default=51)
    ap.add_argument("--max-reads", type=int, default=0, help="If >0, process only first N reads (debug)")
    args = ap.parse_args()

    wl = read_whitelist(args.whitelist)
    feat_seq_to_id = read_feature_ref(args.feature_ref)

    # Assign integer ids
    feat_ids = sorted(set(feat_seq_to_id.values()))
    feat_id_to_i = {fid: i for i, fid in enumerate(feat_ids)}

    cell_to_j = {}
    cell_list = []

    counts = defaultdict(int)
    kept = 0
    total = 0

    r1_it = iter_fastq_seqs(args.r1)
    r2_it = iter_fastq_seqs(args.r2)

    for r1_seq, r2_seq in zip(r1_it, r2_it):
        total += 1
        if args.max_reads and total > args.max_reads:
            break

        cb = r1_seq[args.cb_start:args.cb_end]
        if cb not in wl:
            continue

        fb = r2_seq[args.fb_start:args.fb_end]
        fid = feat_seq_to_id.get(fb)
        if fid is None:
            continue

        if cb not in cell_to_j:
            cell_to_j[cb] = len(cell_list)
            cell_list.append(cb)

        j = cell_to_j[cb]
        i = feat_id_to_i[fid]
        counts[(j, i)] += 1
        kept += 1

        if kept % 1_000_000 == 0:
            print(f"... kept {kept:,} reads (total {total:,})")

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    # write barcodes/features like 10x
    with open(out / "barcodes.tsv", "wt") as f:
        for cb in cell_list:
            f.write(cb + "\n")

    with open(out / "features.tsv", "wt") as f:
        # features.tsv: id \t name \t type (simple)
        for fid in feat_ids:
            f.write(f"{fid}\t{fid}\tCRISPR\n")

    write_mtx(out, counts, n_cells=len(cell_list), n_feats=len(feat_ids))

    print(f"Done. total_reads={total:,} kept_reads={kept:,} cells={len(cell_list)} features={len(feat_ids)}")
    print(f"Outputs: {out}/matrix.mtx, features.tsv, barcodes.tsv")

if __name__ == "__main__":
    main()
