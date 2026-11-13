#!/usr/bin/env python3
"""
collect_final_outputs.py — copy a run's figure PDFs and PNGs into one folder.

    python pipeline/collect_final_outputs.py --run_dir results/min307/SBS

Writes <run_dir>/final_output/ (replaced on every run), with copies of:

    cluster_summaries/        cluster_summary, main_clusters_with_dendrogram,
                              small_clusters_summary (.pdf + .png; Step 2b)
    species_compound_matrix/  combined figure (.pdf + .png; Step 6b) and the
                              human / mouse / other_species matrices (Step 6)
    cosine_similarity/        the three eSS vs COSMIC heatmaps (Step 5)
    dendrogram/               the clustering dendrogram (.pdf + .png; Step 1)
    signature_profiles/       eSS and small-cluster profile plots (Step 1)
    decomposition/            COSMIC decomposition plots (Step 8, SBS)
    run_parameters.txt        the run's settings and results (Step 7)
    MANIFEST.txt              every file, where it was copied from, and any
                              expected file that is missing

The originals stay where they are; these are copies. A missing file is listed
in MANIFEST.txt and reported as a warning, and the rest are still copied.
"""

import argparse
import datetime
import os
import shutil
import sys

# (folder in final_output, source path relative to the run folder, required)
# required=False: only made in some runs (e.g. Chrome missing, or not SBS).
FILES = [
    ("cluster_summaries", "reports/cluster_summary.pdf", False),
    ("cluster_summaries", "reports/cluster_summary.png", False),
    ("cluster_summaries", "reports/main_clusters_with_dendrogram.pdf", False),
    ("cluster_summaries", "reports/main_clusters_with_dendrogram.png", False),
    ("cluster_summaries", "reports/small_clusters_summary.pdf", False),
    ("cluster_summaries", "reports/small_clusters_summary.png", False),
    ("species_compound_matrix", "species_compound_matrix_combined.pdf", False),
    ("species_compound_matrix", "species_compound_matrix_combined.png", False),
    ("species_compound_matrix", "species_compound_matrix_human.pdf", True),
    ("species_compound_matrix", "species_compound_matrix_mouse.pdf", True),
    ("species_compound_matrix", "species_compound_matrix_other_species.pdf", True),
    ("cosine_similarity", "cosine_similarity_heatmap.pdf", True),
    ("cosine_similarity", "cosine_similarity_heatmap_high_similarity_only.pdf", True),
    ("cosine_similarity", "cosine_similarity_heatmap_no_artifacts.pdf", True),
    ("dendrogram", "interactive_clustering/all_esignature_models_dendrogram.pdf", True),
    ("dendrogram", "interactive_clustering/all_esignature_models_dendrogram.png", True),
    ("signature_profiles", "main_clusters/main_clusters_hierarchical_plots/{T}_*_plots_main_clusters_hierarchical.pdf", True),
    ("signature_profiles", "small_clusters/small_clusters_hierarchical_plots/{T}_*_plots_small_clusters_hierarchical.pdf", False),
    ("decomposition", "decomposition/Decompose_Solution/SBS96_Decomposition_Plots.pdf", False),
    ("", "run_parameters.txt", True),
]


def resolve(run_dir, pattern, mutation_type):
    import glob
    pattern = pattern.replace("{T}", mutation_type)
    matches = sorted(glob.glob(os.path.join(run_dir, pattern)))
    return matches[0] if matches else None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run_dir", required=True, help="<run folder>/<TYPE>, e.g. results/min307/SBS")
    parser.add_argument("--mutation_type", default=None,
                        help="SBS, DBS or ID (default: the run folder's name)")
    args = parser.parse_args()

    run_dir = os.path.abspath(args.run_dir)
    mutation_type = (args.mutation_type or os.path.basename(run_dir)).upper()
    out = os.path.join(run_dir, "final_output")
    if os.path.isdir(out):
        shutil.rmtree(out)          # only this generated folder, so no stale copies remain
    os.makedirs(out)

    copied, missing = [], []
    for folder, pattern, required in FILES:
        src = resolve(run_dir, pattern, mutation_type)
        if src is None:
            if required or mutation_type == "SBS":
                missing.append(pattern.replace("{T}", mutation_type))
            continue
        dest_dir = os.path.join(out, folder)
        os.makedirs(dest_dir, exist_ok=True)
        dest = os.path.join(dest_dir, os.path.basename(src))
        shutil.copy2(src, dest)
        copied.append((os.path.relpath(dest, out), os.path.relpath(src, run_dir)))

    lines = [f"Final outputs of {run_dir}",
             f"Collected: {datetime.datetime.now().isoformat(timespec='seconds')}",
             "Copies; the originals are in the run folder at the paths shown.", "",
             f"{len(copied)} files:"]
    lines += [f"  {dest:70s} <- {src}" for dest, src in copied]
    if missing:
        lines += ["", f"{len(missing)} expected file(s) not found (not made in this run):"]
        lines += [f"  {m}" for m in missing]
    with open(os.path.join(out, "MANIFEST.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Collected {len(copied)} files into {out}")
    for m in missing:
        print(f"WARNING: not found, not collected: {m}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
