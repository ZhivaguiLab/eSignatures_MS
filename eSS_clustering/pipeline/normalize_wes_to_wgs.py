#!/usr/bin/env python3
"""
normalize_wes_to_wgs.py

A4 (manuscript revision) Pass 1: correct WES samples onto their own species'
WGS trinucleotide-opportunity basis, before any cross-species comparison is
attempted. Cross-species opportunity normalization (Pass 2) is deliberately
NOT done here -- genome builds for human/mouse are confirmed (GRCh38, mm10),
but rat/chicken/celegans builds are not yet confirmed and chicken/celegans
opportunity tables don't exist yet. Rat, chicken, and celegans samples in
this dataset are already 100% WGS, so they need no correction in this pass
regardless.

For each WES sample, the raw counts are rescaled per SBS-96 context by
(that genome's WGS trinucleotide relative frequency) / (that genome's WES
[exome-capture] trinucleotide relative frequency), then renormalized to sum
to the sample's original total mutation count -- i.e. this only reshapes the
profile, it does not change the sample's total mutation burden. WGS samples
are left untouched (ratio == 1 in every context that matters).

Produces a complete new data/input/SBS/-equivalent directory: human and
mouse files are corrected, all other species' files are copied unchanged
(they're required for the pipeline's per-species discovery either way, and
must stay byte-identical since they're not being corrected in this pass).

Usage
-----
    python pipeline/normalize_wes_to_wgs.py \\
        --input_dir data/input/SBS \\
        --output_dir data/input_wes_to_wgs_normalized/SBS \\
        --context_dist_dir /path/to/context-distributions
"""

import argparse
import os
import re
import shutil

import pandas as pd

BUILD = {"human": "GRCh38", "mouse": "mm10"}
CORRECT_SPECIES = set(BUILD.keys())
ALL_SPECIES = ["human", "mouse", "rat", "chicken", "celegans"]


def trinuc_rel_freq(context_dist_dir, genome, exome):
    suffix = "_96_exome.csv" if exome else "_96.csv"
    path = os.path.join(context_dist_dir, f"context_counts_{genome}{suffix}")
    df = pd.read_csv(path, index_col=0)
    if "Y" in df.columns:
        df = df.drop(columns=["Y"])
    totals = df.sum(axis=1)
    return totals / totals.sum()


def context_to_trinuc(label):
    m = re.match(r"^(.)\[(.)>(.)\](.)$", label)
    five, ref, _alt, three = m.groups()
    return five + ref + three


def wes_to_wgs_ratio(context_dist_dir, contexts, genome):
    wes = trinuc_rel_freq(context_dist_dir, genome, exome=True)
    wgs = trinuc_rel_freq(context_dist_dir, genome, exome=False)
    ratio_by_trinuc = wgs / wes
    return pd.Series({c: ratio_by_trinuc[context_to_trinuc(c)] for c in contexts})


def tech_of(sample_name):
    if re.search("exome", sample_name, re.I):
        return "WES"
    if re.search("genome", sample_name, re.I):
        return "WGS"
    return None


def correct_species_file(input_dir, output_dir, species, context_dist_dir):
    raw_path = os.path.join(input_dir, f"filtered_{species}_307.txt")
    raw = pd.read_csv(raw_path, sep='\t', index_col=0)
    contexts = list(raw.index)

    ratio = wes_to_wgs_ratio(context_dist_dir, contexts, BUILD[species])

    n_wes = n_wgs = n_unknown = 0
    corrected = raw.copy().astype(float)
    for col in raw.columns:
        tech = tech_of(col)
        if tech == "WES":
            n_wes += 1
            counts = raw[col].astype(float)
            total = counts.sum()
            rescaled = counts * ratio.reindex(contexts).values
            corrected[col] = rescaled / rescaled.sum() * total  # reshape only; preserve total burden
        elif tech == "WGS":
            n_wgs += 1
        else:
            n_unknown += 1

    print(f"  {species}: {n_wgs} WGS (unchanged), {n_wes} WES (corrected to {BUILD[species]} WGS basis), "
          f"{n_unknown} unrecognized tech (left unchanged)")

    corrected.to_csv(os.path.join(output_dir, f"filtered_{species}_307.txt"), sep='\t')
    col_sums = corrected.sum(axis=0)
    normalized = corrected.div(col_sums, axis=1)
    normalized.to_csv(os.path.join(output_dir, f"normalized_filtered_{species}_307.tsv"), sep='\t')


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input_dir", default="data/input/SBS")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--context_dist_dir", required=True,
                        help="Directory containing context_counts_<genome>_96[_exome].csv files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for species in ALL_SPECIES:
        raw_path = os.path.join(args.input_dir, f"filtered_{species}_307.txt")
        if not os.path.exists(raw_path):
            continue
        if species in CORRECT_SPECIES:
            correct_species_file(args.input_dir, args.output_dir, species, args.context_dist_dir)
        else:
            # Already WGS-only in this dataset -- copy raw + normalized unchanged.
            shutil.copy(raw_path, os.path.join(args.output_dir, f"filtered_{species}_307.txt"))
            norm_src = os.path.join(args.input_dir, f"normalized_filtered_{species}_307.tsv")
            if os.path.exists(norm_src):
                shutil.copy(norm_src, os.path.join(args.output_dir, f"normalized_filtered_{species}_307.tsv"))
            print(f"  {species}: copied unchanged (already 100% WGS)")

    print(f"\nDone. Output -> {args.output_dir}")


if __name__ == "__main__":
    main()
