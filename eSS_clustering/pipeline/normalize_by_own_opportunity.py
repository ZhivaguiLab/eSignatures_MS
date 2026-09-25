#!/usr/bin/env python3
"""
normalize_by_own_opportunity.py

A4 alternative to the liftover/translation approach (normalize_wes_to_wgs.py):
instead of rescaling WES samples onto a chosen target genome's (WGS) basis,
normalize EVERY sample -- WGS included -- by its own species+build+sequencing
-technology trinucleotide opportunity, independently. No common reference
genome is chosen or needed; each sample only ever uses its own opportunity
table.

    rate[c]    = raw_count[c] / opportunity[c]   (own species/build/tech table)
    profile[c] = rate[c] / sum(rate)              (renormalize to sum 1)

This produces the pipeline's per-sample "pre-normalized" file directly.
Raw counts (filtered_<species>_307.txt) are left completely unchanged --
only the normalized/proportional representation used for clustering and
consensus averaging is affected. This is a mathematically different result
from normalize_wes_to_wgs.py, not a refinement of it -- WGS samples are
corrected here (never in the liftover approach), and WES samples are
corrected differently (normalized by their own opportunity only, not
rescaled onto WGS's).

Species without a confirmed opportunity table (chicken, celegans, as of
this writing) are left on the original raw/raw.sum() normalization,
unchanged, and clearly reported as such -- this is a partial-coverage run,
not a full cross-species correction.

Usage
-----
    python pipeline/normalize_by_own_opportunity.py \\
        --input_dir data/input/SBS \\
        --output_dir data/input_opportunity_normalized/SBS \\
        --context_dist_dir /path/to/context-distributions
"""

import argparse
import os
import re
import shutil

import pandas as pd

BUILD = {"human": "GRCh38", "mouse": "mm10", "rat": "rn7"}
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


def own_opportunity(context_dist_dir, contexts, genome, exome):
    """Per-context relative opportunity for this exact (genome, tech), no ratio to anything else."""
    rel_freq = trinuc_rel_freq(context_dist_dir, genome, exome)
    return pd.Series({c: rel_freq[context_to_trinuc(c)] for c in contexts})


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
    genome = BUILD[species]

    n_wes = n_wgs = n_unknown = 0
    norm_cols = {}
    for col in raw.columns:
        tech = tech_of(col)
        counts = raw[col].astype(float)
        if tech in ("WES", "WGS"):
            opp = own_opportunity(context_dist_dir, contexts, genome, exome=(tech == "WES"))
            rate = counts / opp.reindex(contexts).values
            norm_cols[col] = rate / rate.sum()
            if tech == "WES":
                n_wes += 1
            else:
                n_wgs += 1
        else:
            # Unrecognized tech -- fall back to plain proportion, can't pick an opportunity table.
            norm_cols[col] = counts / counts.sum()
            n_unknown += 1

    print(f"  {species}: {n_wgs} WGS (own {genome}-WGS opportunity), "
          f"{n_wes} WES (own {genome}-exome opportunity), "
          f"{n_unknown} unrecognized tech (plain proportion, no correction)")

    # Raw counts unchanged.
    raw.to_csv(os.path.join(output_dir, f"filtered_{species}_307.txt"), sep='\t')
    norm_df = pd.DataFrame(norm_cols).reindex(contexts)
    norm_df.to_csv(os.path.join(output_dir, f"normalized_filtered_{species}_307.tsv"), sep='\t')


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
            # No opportunity table for this species yet -- copy raw + normalized unchanged.
            shutil.copy(raw_path, os.path.join(args.output_dir, f"filtered_{species}_307.txt"))
            norm_src = os.path.join(args.input_dir, f"normalized_filtered_{species}_307.tsv")
            if os.path.exists(norm_src):
                shutil.copy(norm_src, os.path.join(args.output_dir, f"normalized_filtered_{species}_307.tsv"))
            print(f"  {species}: copied unchanged (no opportunity table available yet)")

    print(f"\nDone. Output -> {args.output_dir}")


if __name__ == "__main__":
    main()
