#!/usr/bin/env python3
"""
build_translation_verification_matrices.py

For every WES sample that was corrected by normalize_wes_to_wgs.py, build a
before/after SBS-96 profile matrix pair (suitable for sigProfilerPlotting)
and a per-sample cosine-similarity table, so the effect of the correction
can be checked sample by sample -- both visually and numerically.

Usage
-----
    python pipeline/build_translation_verification_matrices.py \\
        --orig_dir data/input/SBS \\
        --corrected_dir data/input_wes_to_wgs_normalized/SBS \\
        --output_dir results_translation_verification/SBS

Then plot both matrices, e.g.:
    python -c "
import sigProfilerPlotting as sigPlt
sigPlt.plotSBS('results_translation_verification/SBS/translated_samples_BEFORE.tsv',
                'results_translation_verification/SBS/plots', 'translated_BEFORE', '96',
                percentage=True)
sigPlt.plotSBS('results_translation_verification/SBS/translated_samples_AFTER.tsv',
                'results_translation_verification/SBS/plots', 'translated_AFTER', '96',
                percentage=True)
"
"""

import argparse
import os
import re

import numpy as np
import pandas as pd


def cosine(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--orig_dir", required=True)
    parser.add_argument("--corrected_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--species", nargs="+", default=["human", "mouse"],
                        help="Species to include (must have both orig and corrected files).")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    before_cols, after_cols, meta_rows = {}, {}, []

    for species in args.species:
        orig = pd.read_csv(os.path.join(args.orig_dir, f"filtered_{species}_307.txt"),
                           sep='\t', index_col=0)
        corr = pd.read_csv(os.path.join(args.corrected_dir, f"filtered_{species}_307.txt"),
                           sep='\t', index_col=0)
        contexts = list(orig.index)

        wes_cols = [c for c in orig.columns if re.search("exome", c, re.I)]
        for col in wes_cols:
            tag = f"{species}__{col}"
            o = orig[col].reindex(contexts).astype(float)
            c = corr[col].reindex(contexts).astype(float)
            before_cols[tag] = o / o.sum()
            after_cols[tag] = c / c.sum()
            meta_rows.append({"sample": tag, "species": species, "original_name": col})

    if not before_cols:
        raise ValueError("No WES samples found for the given species -- nothing to compare.")

    common_idx = list(next(iter(before_cols.values())).index)
    before_df = pd.DataFrame(before_cols).reindex(common_idx)
    after_df = pd.DataFrame(after_cols).reindex(common_idx)
    before_df.index.name = "MutationType"
    after_df.index.name = "MutationType"

    before_df.to_csv(os.path.join(args.output_dir, "translated_samples_BEFORE.tsv"), sep='\t')
    after_df.to_csv(os.path.join(args.output_dir, "translated_samples_AFTER.tsv"), sep='\t')

    meta = pd.DataFrame(meta_rows)
    meta.to_csv(os.path.join(args.output_dir, "translated_samples_metadata.tsv"), sep='\t', index=False)

    rows = [
        {"sample": tag, "species": tag.split("__", 1)[0], "original_name": tag.split("__", 1)[1],
         "cosine_before_vs_after": cosine(before_df[tag].values, after_df[tag].values)}
        for tag in before_cols
    ]
    sim_df = pd.DataFrame(rows).sort_values("cosine_before_vs_after")
    sim_df.to_csv(os.path.join(args.output_dir, "cosine_similarity_before_vs_after.tsv"),
                 sep='\t', index=False)

    print(f"Total translated samples: {len(before_cols)} "
          f"({meta.groupby('species').size().to_dict()})")
    print(f"\nCosine similarity (before vs after) distribution:")
    print(sim_df["cosine_before_vs_after"].describe().to_string())
    print(f"\nMost affected (lowest cosine) samples:")
    print(sim_df.head(10).to_string(index=False))
    print(f"\nOutputs -> {args.output_dir}")


if __name__ == "__main__":
    main()
