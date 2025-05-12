#!/usr/bin/env python3
"""
Extract PolyPhen-damaging variants from VCF files for **multiple**
individuals, compare sums in ‘real’ regions vs. 100 random replicates,
and report the difference.

Usage:
    python extract_pp_dam_batch.py individuals.txt
"""

from __future__ import annotations
import os
import sys
import subprocess
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd


# ─────────────────────────────  CONSTANT PATHS  ────────────────────────────── #
VCF_DIR           = Path("/lisc/scratch/admixlab/aigerim/ensembl-vep")
VCF_FILE_TEMPLATE = "sstar_subset_filter1_{chr}.vcf.gz"
CHROMOSOMES       = [str(i) for i in range(1, 23)]

BED_FILE = Path("Polyphen_dam_sites.txt")
GENOME_FILE   = Path("genome.txt")                      # BEDTools genome file

BASE_DIR = (Path("/lisc/scratch/admixlab/aigerim/sstar-analysis-main")
            / "African" / "San21")                      # one sub-dir per ID

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)


# ─────────────────────────────  HELPER FUNCTIONS  ──────────────────────────── #
def extract_genotypes(bed_file: Path, individual: str) -> List[pd.DataFrame]:
    """Return one DataFrame per chromosome with genotypes for *individual*."""
    dfs: List[pd.DataFrame] = []

    for chr_ in CHROMOSOMES:
        print(f"[INFO] {individual}: processing chr{chr_}", flush=True)
        vcf_file = VCF_DIR / VCF_FILE_TEMPLATE.format(chr=chr_)

        # temporary filenames share a prefix → easy clean-up later
        prefix       = f"{individual}_chr{chr_}_dam"
        chr_real_bed = Path(f"{prefix}_regions.bed")
        chr_bed = Path(f"{prefix}.bed")
        chr_inter    = Path(f"{prefix}_intersect.bed")

        # ---  filter original BEDs by chromosome  --------------------------
        with bed_file.open() as src, chr_real_bed.open("w") as dst:
            for line in src:
                if line.startswith(f"chr{chr_}\t") or line.startswith(f"{chr_}\t"):
                    dst.write(line)

        with BED_FILE.open() as src, chr_bed.open("w") as dst:
            for line in src:
                if line.startswith(f"chr{chr_}\t") or line.startswith(f"{chr_}\t"):
                    dst.write(line)

        # ---  sort & intersect  --------------------------------------------
        def _sort(in_file: Path, out_file: Path) -> None:
            subprocess.run(
                ["bedtools", "sort", "-i", in_file, "-g", GENOME_FILE],
                stdout=out_file.open("w"), check=True
            )

        sorted_real = Path(f"{prefix}_sorted_regions.bed")
        sorted = Path(f"{prefix}_sorted.bed")
        _sort(chr_real_bed, sorted_real)
        _sort(chr_bed, sorted)

        subprocess.run(
            ["bedtools", "intersect", "-a", sorted_real, "-b", sorted,
             "-wa", "-wb", "-sorted"],
            stdout=chr_inter.open("w"), check=True
        )

        if chr_inter.stat().st_size == 0:          # no overlaps → skip chr
            _cleanup(prefix)
            continue

        # ---  pull the genotypes  ------------------------------------------
        intersect = pd.read_csv(chr_inter, sep="\t", header=None).drop_duplicates()
        temp_pos  = Path(f"{prefix}_positions.txt")
        intersect[[3, 4, 5]].to_csv(temp_pos, sep="\t", index=False, header=False)

        temp_out  = Path(f"{prefix}_genotypes.txt")
        view_cmd  = ["bcftools", "view", "-m2", "-M2", "-v", "snps",
                     "-R", temp_pos, "-s", individual, vcf_file]
        query_cmd = ["bcftools", "query", "-f",
                     '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n']

        with temp_out.open("w") as fout:
            p1 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p2 = subprocess.Popen(query_cmd, stdin=p1.stdout, stdout=fout, stderr=subprocess.PIPE)
            p1.stdout.close()
            p1.communicate()
            p2.communicate()

        if temp_out.stat().st_size == 0:
            _cleanup(prefix)
            continue

        df = pd.read_csv(
            temp_out, sep="\t",
            names=["CHROM", "POS", "REF", "ALT", "GT"],
            dtype={"CHROM": str, "POS": int}
        ).drop_duplicates(subset=["CHROM", "POS"])

        dfs.append(df)
        _cleanup(prefix)

    return dfs


def _cleanup(prefix: str) -> None:
    """Remove all temporary files whose names start with *prefix*."""
    for f in Path.cwd().glob(f"{prefix}*"):
        try:
            f.unlink()
        except FileNotFoundError:
            pass


def process_genotypes(dfs: List[pd.DataFrame], csv_out: Path | None = None) -> int:
    """Merge chromosome DataFrames, recode GT, optionally write CSV, return sum."""
    if not dfs:
        return 0

    df = pd.concat(dfs, ignore_index=True)
    df["Recoded_GT"] = df["GT"].map(
        lambda gt: 0 if gt in {"./.", ".|.", "0/0", "0|0"} else 1
    )

    if csv_out:
        df.to_csv(csv_out, index=False)

    return int(df["Recoded_GT"].sum(skipna=True))


def run_for_individual(individual: str) -> None:
    print("=" * 72)
    print(f"[INFO] Starting {individual}", flush=True)

    indiv_dir = BASE_DIR / individual
    real_bed  = indiv_dir / f"{individual}_overlap.bed"

    if not real_bed.exists():
        print(f"[ERROR] Missing BED: {real_bed} – skipping.", flush=True)
        return

    # ---  real regions  -----------------------------------------------------
    real_dfs = extract_genotypes(real_bed, individual)
    real_sum = process_genotypes(real_dfs)

    # ---  100 random replicates  -------------------------------------------
    rand_sums, missing = [], []
    for rep in range(1, 101):
        rand_bed = indiv_dir / f"{individual}.random.replicate{rep}_overlap.bed"
        if not rand_bed.exists():
            missing.append(rep)
            continue
        dfs = extract_genotypes(rand_bed, individual)
        rand_sums.append(process_genotypes(dfs))

    if not rand_sums:
        print(f"[ERROR] No random replicates processed for {individual}")
        return

    # ---  summary  ----------------------------------------------------------
    mean_rand = float(np.mean(rand_sums))
    sd_rand   = float(np.std(rand_sums, ddof=1))

    print(f"  • real-region sum        : {real_sum}")
    print(f"  • mean random (n={len(rand_sums)}) : {mean_rand:.2f} ± {sd_rand:.2f}")
    print(f"  • difference (real-mean) : {real_sum - mean_rand:.2f}")
    if missing:
        print(f"    (missing replicates: {missing})")

    out_csv = RESULTS_DIR / f"{individual}_pp_dam_100.csv"
    rows = [{"type": "real", "sum": real_sum}] + \
           [{"type": "random", "sum": s} for s in rand_sums]
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"[INFO] Wrote {out_csv}\n")


# ────────────────────────────────  MAIN  ───────────────────────────────────── #
def main() -> None:
    if len(sys.argv) != 2:
        print("Usage: python extract_pp_dam_batch.py <INDIVIDUAL_LIST_FILE>")
        sys.exit(1)

    list_file = Path(sys.argv[1])
    if not list_file.exists():
        sys.exit(f"[ERROR] File not found: {list_file}")

    ids = [line.strip() for line in list_file.open() if line.strip()]
    if not ids:
        sys.exit("[ERROR] The list file is empty.")

    for ind in ids:
        run_for_individual(ind)


if __name__ == "__main__":
    main()
