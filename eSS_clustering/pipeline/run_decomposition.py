#!/usr/bin/env python3
"""
run_decomposition.py — decompose COSMIC signatures into the eSS of a run
with SigProfilerAssignment (decompose_fit).

Each COSMIC signature is reconstructed from the run's eSS consensus profiles.
A COSMIC signature whose best reconstruction has cosine similarity below the
novelty threshold (default 0.8, SigProfilerAssignment's own default) is kept
as itself (not explained by the eSS). The decomposed signatures are then
assigned to the samples. The COSMIC-to-eSS map
(Decompose_Solution/De_Novo_map_to_COSMIC_SBS96.csv, used for Figure 5) does
not depend on the samples; only the activities do.

Usage
-----
    # eSS from results/min307_cos0.90/SBS/; COSMIC is also used as the samples
    python pipeline/run_decomposition.py --run_dir results/min307_cos0.90/SBS

    # assign to the samples clustered in that run, or any other matrix;
    # a different novelty threshold
    python pipeline/run_decomposition.py --run_dir results/min307_cos0.90/SBS \\
        --samples clustered --threshold 0.85

Output (default <run_dir>/decomposition/):
    Decompose_Solution/...           SigProfilerAssignment output
    decomposition_parameters.txt/.json  settings, versions, input hashes and
                                     how many COSMIC signatures the eSS explain

Requires SigProfilerAssignment 1.1.4 (the version the results were verified
with); another version stops the script unless --allow_other_version is given,
since outputs differ between versions.
"""

import argparse
import datetime
import glob
import hashlib
import json
import os
import re
import sys

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_COSMIC = os.path.join(REPO_ROOT, "data", "references", "COSMIC_v3.6_SBS_GRCh38.txt")
EXPECTED_SPA_VERSION = "1.1.4"
ESS_PROFILES = os.path.join("main_clusters", "main_clusters_NORMALIZED_weighted_avg_profiles.tsv")


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def package_version(name):
    from importlib.metadata import version, PackageNotFoundError
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def clustered_samples(run_dir, out_path):
    """
    Write the counts of the samples clustered in run_dir to one matrix.

    The counts come from the run's preprocessed data (data/input_cleaned/...),
    which must still match the input files recorded for the run.
    """
    with open(os.path.join(run_dir, "run_parameters.json")) as f:
        record = json.load(f)
    data_dir = os.path.join(REPO_ROOT, record["inputs"]["data_dir"])
    stamp_path = os.path.join(data_dir, "preprocessing_stamp.json")
    if not os.path.exists(stamp_path):
        sys.exit(f"Preprocessed data for this run not found in {data_dir}; rerun the pipeline.")
    with open(stamp_path) as f:
        stamp = json.load(f)
    if stamp.get("input_files_sha256") != record["inputs"]["input_files_sha256"]:
        sys.exit(f"{data_dir} was rebuilt since this run (different inputs); rerun the pipeline.")

    files = sorted(f for f in glob.glob(os.path.join(data_dir, "*.txt"))
                   if "normaliz" not in os.path.basename(f).lower())
    counts = pd.concat([pd.read_csv(f, sep="\t", index_col=0) for f in files], axis=1)
    expected = record["results"]["samples_clustered"]
    if counts.shape[1] != expected:
        sys.exit(f"Found {counts.shape[1]} samples in {data_dir}, but the run clustered "
                 f"{expected}; rerun the pipeline.")
    counts.index.name = "MutationType"
    counts.to_csv(out_path, sep="\t")
    return out_path, counts.shape[1]


def summarize(output_dir, threshold):
    """How many COSMIC signatures the eSS reconstruct (>= threshold) vs keep as novel."""
    path = os.path.join(output_dir, "Decompose_Solution", "De_Novo_map_to_COSMIC_SBS96.csv")
    d = pd.read_csv(path, skipinitialspace=True)
    d.columns = [c.strip() for c in d.columns]
    explained = d["Global NMF Signatures"].str.contains(r"eSS\d+")
    return {
        "cosmic_signatures": int(len(d)),
        "explained_by_eSS": int(explained.sum()),
        "not_explained": int((~explained).sum()),
        "explained_min_cosine": float(d.loc[explained, "Cosine Similarity"].astype(float).min())
        if explained.any() else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run_dir", help="Pipeline output for one run, e.g. results/min307_cos0.90/SBS")
    parser.add_argument("--esignatures", help="eSS profiles TSV (default: the run's main-cluster profiles)")
    parser.add_argument("--cosmic", default=DEFAULT_COSMIC, help="COSMIC signatures to decompose")
    parser.add_argument("--samples", default="cosmic",
                        help="Samples to assign activities to: 'cosmic' (default; the COSMIC "
                             "file itself, as in SPA_1000_shuffled_libraries/null_library_fit.py), "
                             "'clustered' (the samples clustered in --run_dir) or an SBS96 matrix path")
    parser.add_argument("--output", help="Output directory (default: <run_dir>/decomposition)")
    parser.add_argument("--threshold", type=float, default=0.8,
                        help="Novelty threshold (new_signature_thresh_hold): COSMIC signatures "
                             "whose best eSS reconstruction has cosine below this stay as themselves.")
    parser.add_argument("--no_plots", action="store_true", help="Skip the decomposition plots")
    parser.add_argument("--allow_other_version", action="store_true",
                        help=f"Run with a SigProfilerAssignment version other than {EXPECTED_SPA_VERSION}")
    args = parser.parse_args()

    spa_version = package_version("SigProfilerAssignment")
    if spa_version is None:
        sys.exit("SigProfilerAssignment is not installed (pip install SigProfilerAssignment==1.1.4).")
    if spa_version != EXPECTED_SPA_VERSION and not args.allow_other_version:
        sys.exit(f"SigProfilerAssignment {spa_version} is installed; results were verified with "
                 f"{EXPECTED_SPA_VERSION}, and outputs differ between versions. Install "
                 f"{EXPECTED_SPA_VERSION} or pass --allow_other_version.")

    if not args.run_dir and not args.esignatures:
        parser.error("give --run_dir or --esignatures")
    esignatures = args.esignatures or os.path.join(args.run_dir, ESS_PROFILES)
    output = args.output or os.path.join(args.run_dir or ".", "decomposition")
    os.makedirs(output, exist_ok=True)

    if args.samples == "clustered":
        if not args.run_dir:
            parser.error("--samples clustered needs --run_dir")
        samples, n_samples = clustered_samples(
            args.run_dir, os.path.join(output, "decomposition_input_samples.tsv"))
    elif args.samples == "cosmic":
        # SigProfilerAssignment needs the samples' first column to be named
        # MutationType; the COSMIC file calls it Type.
        cosmic = pd.read_csv(args.cosmic, sep="\t", index_col=0)
        cosmic.index.name = "MutationType"
        samples = os.path.join(output, "decomposition_input_samples.tsv")
        cosmic.to_csv(samples, sep="\t")
        n_samples = cosmic.shape[1]
    else:
        samples = args.samples
        n_samples = pd.read_csv(samples, sep="\t", index_col=0).shape[1]

    settings = dict(
        input_type="matrix",
        genome_build="GRCh38",          # matches the COSMIC GRCh38 reference; no effect on matrix input
        collapse_to_SBS96=False,
        new_signature_thresh_hold=args.threshold,
        nnls_add_penalty=0.05,          # SigProfilerAssignment 1.1.4 defaults, stated explicitly
        nnls_remove_penalty=0.01,
        initial_remove_penalty=0.05,
        add_background_signatures=True,
        connected_sigs=False,
        make_plots=not args.no_plots,
        export_probabilities=False,
        verbose=False,
    )
    print(f"Decomposing {os.path.relpath(args.cosmic, REPO_ROOT)} into "
          f"{os.path.relpath(esignatures, REPO_ROOT)}; assigning to {n_samples} samples")

    from SigProfilerAssignment import Analyzer
    Analyzer.decompose_fit(samples=samples, output=output, signatures=args.cosmic,
                           signature_database=esignatures, **settings)

    summary = summarize(output, args.threshold)
    record = {
        "run_date": datetime.datetime.now().isoformat(timespec="seconds"),
        "run_dir": os.path.relpath(args.run_dir, REPO_ROOT) if args.run_dir else None,
        "inputs": {
            "signatures_decomposed": os.path.relpath(args.cosmic, REPO_ROOT),
            "signature_database": os.path.relpath(esignatures, REPO_ROOT),
            "samples": {"clustered": "clustered samples of run_dir",
                        "cosmic": "COSMIC signatures (same file as decomposed)"}.get(
                            args.samples, os.path.relpath(samples, REPO_ROOT)),
            "n_samples": n_samples,
            "sha256": {os.path.basename(p): sha256(p) for p in [args.cosmic, esignatures, samples]},
        },
        "settings": settings,
        "summary": summary,
        "software": {name: package_version(name) for name in [
            "SigProfilerAssignment", "SigProfilerMatrixGenerator", "SigProfilerPlotting",
            "numpy", "pandas", "scipy"]},
    }
    with open(os.path.join(output, "decomposition_parameters.json"), "w") as f:
        json.dump(record, f, indent=2)
    lines = [
        "COSMIC decomposition into eSS (SigProfilerAssignment decompose_fit)",
        f"Date:     {record['run_date']}",
        f"Run:      {record['run_dir']}",
        "",
        f"Decomposed:        {record['inputs']['signatures_decomposed']}",
        f"Into (database):   {record['inputs']['signature_database']}",
        f"Samples assigned:  {record['inputs']['samples']} ({n_samples})",
        f"Novelty threshold: {args.threshold} (cosine; below it a COSMIC signature stays as itself)",
        "",
        "SETTINGS",
        *[f"  {k}: {v}" for k, v in settings.items()],
        "",
        "SUMMARY",
        f"  COSMIC signatures:           {summary['cosmic_signatures']}",
        f"  Reconstructed from eSS:      {summary['explained_by_eSS']}"
        + (f" (lowest cosine {summary['explained_min_cosine']:.2f})" if summary['explained_min_cosine'] else ""),
        f"  Not explained by eSS:        {summary['not_explained']}",
        "",
        "INPUTS (sha256)",
        *[f"  {k}: {v}" for k, v in record["inputs"]["sha256"].items()],
        "",
        "SOFTWARE",
        *[f"  {k}: {v}" for k, v in record["software"].items()],
    ]
    with open(os.path.join(output, "decomposition_parameters.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Done. {summary['explained_by_eSS']} of {summary['cosmic_signatures']} COSMIC signatures "
          f"reconstructed from eSS (threshold {args.threshold}). Output → {output}")


if __name__ == "__main__":
    main()
