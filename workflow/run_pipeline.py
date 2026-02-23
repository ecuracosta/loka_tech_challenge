#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
import random
import shlex
import subprocess
from pathlib import Path
from typing import Optional


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Local trigger + optional FASTQ downsampling + run Nextflow pipeline."
    )

    p.add_argument("--local", action="store_true", help="Run in local mode (required for now).")
    p.add_argument("--data-dir", default="data", help="Base data directory (default: data).")

    p.add_argument("--downsample", action="store_true", help="Downsample FASTQs used in samplesheet.")
    p.add_argument("--fraction", type=float, default=0.1, help="Downsample fraction per read (0-1).")
    p.add_argument("--seed", type=int, default=0, help="Random seed for downsampling.")

    p.add_argument("--debug-fast", action="store_true", help="Ultra-fast debug mode (very small fraction, 1 sample).")

    p.add_argument("--nextflow-workflow", default="workflow/main.nf")
    p.add_argument("--nextflow-config", default="workflow/nextflow.config")

    p.add_argument("--dry-run", action="store_true", help="Print what would run, but don't execute nextflow.")

    p.add_argument("--star-index", default=None, help="Path to STAR index (required for --do-solo).")
    p.add_argument(
        "--gex-cb-whitelist",
        default="None",
        help='GEX cell-barcode whitelist for STARsolo. Use "None" to disable whitelist (debug). '
             'Example: /work/ref/3M-february-2018.txt',
    )

    p.add_argument("--threads", type=int, default=8, help="Threads passed to pipeline.")

    # Stages (mirror params in Nextflow)
    p.add_argument("--do-solo", action="store_true", help="Run STARsolo stage.")
    p.add_argument("--do-guides", action="store_true", help="Run guide assignment stage.")
    p.add_argument("--do-integrate", action="store_true", help="Run integration stage.")
    p.add_argument("--do-analysis", action="store_true", help="Run downstream analysis stage.")

    # Guide-related inputs
    p.add_argument("--crispr-r1", default=None)
    p.add_argument("--crispr-r2", default=None)
    p.add_argument("--cb-whitelist", default=None)
    p.add_argument("--feature-ref", default=None)

    p.add_argument("--max-cb-mm", type=int, default=2)
    p.add_argument("--max-guide-mm", type=int, default=2)

    p.add_argument(
        "--library",
        choices=["gex", "crispr", "combined"],
        default="gex",
        help="Library type to use for samplesheet inference (default: gex).",
    )
    p.add_argument(
        "--out-fastq-dir",
        default="data/downsampled_fastq",
        help="Where to write downsampled FASTQs (default: data/downsampled_fastq).",
    )

    # Integration params
    p.add_argument("--min-guide-count", type=int, default=3)
    p.add_argument("--max-multi-fraction", type=float, default=0.1)

    p.add_argument(
        "--nf-args",
        default="",
        help='Extra args passed verbatim to nextflow run (e.g. "--foo 1 --bar baz").',
    )

    return p.parse_args()

from pathlib import Path

def resolve_data_dir(data_dir: str) -> Path:
    """
    Make data_dir robust to different working directories.

    Priority:
      1) If user passes an absolute/relative path that exists -> use it.
      2) If running from /work/workflow and data is in /work/data -> auto-fallback.
      3) If running from anywhere and /work/data exists -> use it.
    """
    p = Path(data_dir)

    # 1) as provided (relative to cwd)
    if p.exists():
        return p.resolve()

    # 2) common container layout: code in /work/workflow, data in /work/data
    if data_dir == "data":
        p2 = Path("/work/data")
        if p2.exists():
            return p2.resolve()

    # 3) fallback: try interpreting as /work/<data_dir>
    p3 = Path("/work") / data_dir
    if p3.exists():
        return p3.resolve()

    # If nothing works, raise with helpful message
    raise FileNotFoundError(
        f"data_dir not found: '{data_dir}' (tried: {p}, /work/data, {p3})"
    )

def downsample_fastq(
    input_file: Path,
    output_file: Path,
    fraction: float,
    rng: random.Random,
) -> tuple[int, int]:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    in_fh = gzip.open(input_file, "rt") if input_file.name.endswith(".gz") else open(input_file, "rt")
    out_fh = gzip.open(output_file, "wt") if output_file.name.endswith(".gz") else open(output_file, "wt")

    kept = 0
    total = 0
    try:
        while True:
            rec = [in_fh.readline() for _ in range(4)]
            if not rec[0]:
                break
            total += 1
            if rng.random() < fraction:
                out_fh.writelines(rec)
                kept += 1
    finally:
        in_fh.close()
        out_fh.close()

    return total, kept


def infer_sample_pairs(fastq_dir: Path) -> list[tuple[str, Path, Path]]:
    r1_files = sorted(fastq_dir.glob("*_R1_*.fastq*"))
    pairs: list[tuple[str, Path, Path]] = []

    for r1 in r1_files:
        r2_name = r1.name.replace("_R1_", "_R2_", 1)
        r2 = fastq_dir / r2_name
        if not r2.exists():
            print(f"WARNING: could not find R2 pair for {r1.name} (expected {r2_name}); skipping")
            continue

        sample = r1.name.split("_R1_", 1)[0]
        pairs.append((sample, r1, r2))

    return pairs


def write_samplesheet(pairs: list[tuple[str, Path, Path]], out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing samplesheet to {out_path}")
    with open(out_path, "w") as f:
        f.write("sample,fastq_1,fastq_2\n")
        for sample, r1, r2 in pairs:
            f.write(f"{sample},{r1.resolve().as_posix()},{r2.resolve().as_posix()}\n")
    return out_path


def _require_path(flag_name: str, value: Optional[str]) -> Path:
    if not value:
        raise SystemExit(f"ERROR: {flag_name} is required for this mode.")
    p = Path(value)
    if not p.exists():
        raise SystemExit(f"ERROR: {flag_name} does not exist: {p}")
    return p


def run_pipeline(
    nextflow_workflow: str,
    nextflow_config: str,
    samplesheet: Path,
    threads: int,
    *,
    do_solo: bool,
    star_index: Optional[str],
    gex_cb_whitelist: str,
    do_guides: bool,
    do_integrate: bool,
    do_analysis: bool,
    crispr_r1: Optional[str],
    crispr_r2: Optional[str],
    cb_whitelist: Optional[str],
    feature_ref: Optional[str],
    max_cb_mm: int,
    max_guide_mm: int,
    min_guide_count: int,
    max_multi_fraction: float,
    extra_env: Optional[dict[str, str]] = None,
    dry_run: bool = False,
    nf_args: str = "",
) -> int:

    cmd: list[str] = [
        "nextflow", "run", nextflow_workflow,
        "-c", nextflow_config,
        "--samplesheet", str(samplesheet),
        "--threads", str(threads),
        "--gex_cb_whitelist", str(gex_cb_whitelist),
        "-with-report",
    ]

    if do_solo:
        if not star_index:
            raise SystemExit("ERROR: --do-solo requires --star-index")
        cmd += ["--do_solo", "--star_index", star_index]

    if do_guides:
        _require_path("--crispr-r1", crispr_r1)
        _require_path("--crispr-r2", crispr_r2)
        _require_path("--cb-whitelist", cb_whitelist)
        _require_path("--feature-ref", feature_ref)

        cmd += [
            "--do_guides",
            "--crispr_r1", str(Path(crispr_r1).resolve()),
            "--crispr_r2", str(Path(crispr_r2).resolve()),
            "--cb_whitelist", str(Path(cb_whitelist).resolve()),
            "--feature_ref", str(Path(feature_ref).resolve()),
            "--max_cb_mm", str(max_cb_mm),
            "--max_guide_mm", str(max_guide_mm),
        ]

    if do_integrate:
        cmd += [
            "--do_integrate",
            "--min_guide_count", str(min_guide_count),
            "--max_multi_fraction", str(max_multi_fraction),
        ]

    if do_analysis:
        cmd += ["--do_analysis"]

    if nf_args:
        cmd += shlex.split(nf_args)

    print("Running pipeline:\n  " + " ".join(cmd))

    if dry_run:
        print("(dry-run) Not executing Nextflow.")
        return 0

    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)

    res = subprocess.run(cmd, env=env)
    return res.returncode


def main() -> int:
    print(">>> ENTERED main()")
    print(f">>> run_pipeline.py = {Path(__file__).resolve()}")

    args = parse_args()

    if not args.local:
        print("Nothing to do: only --local is implemented right now.")
        return 0

    if args.debug_fast:
        print(">>> DEBUG FAST MODE ENABLED")
        args.downsample = True
        args.fraction = 0.0001
        args.library = "gex"
        args.threads = min(args.threads, 2)

    if not (0.0 < args.fraction <= 1.0):
        raise SystemExit(f"ERROR: --fraction must be in (0,1], got {args.fraction}")

    # 1) Decide base FASTQ dir to use
    data_dir = resolve_data_dir(args.data_dir)
    if args.library == "combined":
        fastq_dir = data_dir / "fastq"
    elif args.library == "gex":
        fastq_dir = data_dir / "fastq" / "SC3_v3_NextGem_DI_CRISPR_10K_fastqs" / "SC3_v3_NextGem_DI_CRISPR_10K_gex_fastqs"
    else:  # crispr
        fastq_dir = data_dir / "fastq" / "SC3_v3_NextGem_DI_CRISPR_10K_fastqs" / "SC3_v3_NextGem_DI_CRISPR_10K_crispr_fastqs"

    if not fastq_dir.exists():
        raise SystemExit(f"FASTQ dir not found: {fastq_dir} (exists={fastq_dir.exists()})")

    print(f">>> chosen fastq_dir = {fastq_dir}")

    # 2) infer pairs
    pairs = infer_sample_pairs(fastq_dir)
    pairs = [(sample, r1.resolve(), r2.resolve()) for sample, r1, r2 in pairs]
    if args.debug_fast:
        pairs = pairs[:1]
    print(f">>> inferred pairs = {len(pairs)}")
    if not pairs:
        raise SystemExit(f"No R1/R2 pairs inferred in {fastq_dir}")

    # 3) optional downsample
    if args.downsample:
        out_dir = Path(args.out_fastq_dir) / args.library
        out_dir.mkdir(parents=True, exist_ok=True)

        rng = random.Random(args.seed)
        new_pairs = []

        for sample, r1, r2 in pairs:
            r1_out = out_dir / r1.name
            r2_out = out_dir / r2.name

            # REUSE if already exists and non-empty
            if r1_out.exists() and r2_out.exists() and r1_out.stat().st_size > 0 and r2_out.stat().st_size > 0:
                print(f"Reusing existing downsampled pair for {sample}:")
                print(f"  - {r1_out}")
                print(f"  - {r2_out}")
            else:
                t1, k1 = downsample_fastq(r1, r1_out, fraction=args.fraction, rng=rng)
                t2, k2 = downsample_fastq(r2, r2_out, fraction=args.fraction, rng=rng)
                print(f"  - {sample}: R1 kept {k1}/{t1}, R2 kept {k2}/{t2}")

            new_pairs.append((sample, r1_out.resolve(), r2_out.resolve()))

        pairs = new_pairs
        fastq_dir = out_dir
        print(f">>> using downsampled fastq_dir = {fastq_dir}")

    # 4) write samplesheet
    repo_root = Path(__file__).resolve().parents[0]  # /work/workflow
    samplesheet_path = write_samplesheet(pairs, repo_root / "outputs" / "samplesheet.csv")

    # 5) run Nextflow
    rc = run_pipeline(
        args.nextflow_workflow,
        args.nextflow_config,
        samplesheet=samplesheet_path,
        threads=args.threads,
        do_solo=args.do_solo,
        star_index=args.star_index,
        gex_cb_whitelist=args.gex_cb_whitelist,
        do_guides=args.do_guides,
        do_integrate=args.do_integrate,
        do_analysis=args.do_analysis,
        crispr_r1=args.crispr_r1,
        crispr_r2=args.crispr_r2,
        cb_whitelist=args.cb_whitelist,
        feature_ref=args.feature_ref,
        max_cb_mm=args.max_cb_mm,
        max_guide_mm=args.max_guide_mm,
        min_guide_count=args.min_guide_count,
        max_multi_fraction=args.max_multi_fraction,
        extra_env={"FASTQ_DIR": str(fastq_dir)},
        dry_run=args.dry_run,
        nf_args=args.nf_args,
    )

    if rc != 0:
        print(f"Pipeline failed (exit code {rc})")
    else:
        print("Pipeline finished OK")

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
