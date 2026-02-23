#!/usr/bin/env python3
import argparse
import gzip
import random
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Downsample paired FASTQ files by record.")
    p.add_argument("--r1", required=True, type=Path)
    p.add_argument("--r2", required=True, type=Path)
    p.add_argument("--out1", required=True, type=Path)
    p.add_argument("--out2", required=True, type=Path)
    p.add_argument("--fraction", type=float, default=0.1)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def open_fastq(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def open_out(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "wt")


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    with open_fastq(args.r1) as f1, open_fastq(args.r2) as f2, open_out(args.out1) as o1, open_out(args.out2) as o2:
        total = kept = 0
        while True:
            rec1 = [f1.readline() for _ in range(4)]
            rec2 = [f2.readline() for _ in range(4)]
            if not rec1[0] or not rec2[0]:
                break
            total += 1
            if rng.random() < args.fraction:
                o1.writelines(rec1)
                o2.writelines(rec2)
                kept += 1

    print(f"Downsampled: total_pairs={total} kept_pairs={kept} fraction={args.fraction} seed={args.seed}")


if __name__ == "__main__":
    main()