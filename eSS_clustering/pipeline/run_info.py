#!/usr/bin/env python3
"""
run_info.py — run naming and the run parameter summary for run_pipeline.sh.

    python pipeline/run_info.py name --mutation_type SBS \\
        [--min_mutations 250 | per-species] [--custom_thresholds aai-split] \\
        [--normalization wes-to-wgs]
        → prints e.g. "min307", "min250", "min-per-species" or
          "min307_aai-split_wes-to-wgs"

The name records the minimum mutation count and the optional extras. It does
not include the clustering cosine threshold (0.9 by default), which is
recorded in run_parameters.txt; runs that differ only in that threshold would
share a folder.

    python pipeline/run_info.py summary --output_dir results/min307/SBS \\
        --cosmic_profiles data/references/COSMIC_v3.6_SBS_GRCh38.txt --cosmic_threshold 0.85
        → adds the COSMIC comparison to run_parameters.json and writes a
          readable run_parameters.txt next to it
"""

import argparse
import hashlib
import json
import os

import numpy as np
import pandas as pd
import yaml

from utils import min_mutations as cutoff_config

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_PREPROCESSING = os.path.join(REPO_ROOT, "config", "preprocessing.yaml")


def run_name(mutation_type, preprocessing_config=DEFAULT_PREPROCESSING,
             custom_thresholds="none", normalization="none", min_mutations=None):
    with open(preprocessing_config) as f:
        config = yaml.safe_load(f) or {}
    cutoffs = cutoff_config.resolve(config.get(mutation_type), min_mutations)
    parts = [cutoff_config.label(cutoffs)]
    if custom_thresholds and custom_thresholds != "none":
        parts.append(custom_thresholds if ":" not in custom_thresholds else "custom")
    if normalization and normalization != "none":
        parts.append(normalization)
    return "_".join(parts)


def sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def cosmic_comparison(output_dir, cosmic_profiles, threshold):
    profiles = pd.read_csv(os.path.join(output_dir, "main_clusters",
                                        "main_clusters_NORMALIZED_weighted_avg_profiles.tsv"),
                           sep="\t", index_col=0)
    cosmic = pd.read_csv(cosmic_profiles, sep="\t", index_col=0).reindex(profiles.index)
    a = profiles.values / np.linalg.norm(profiles.values, axis=0)
    b = cosmic.values / np.linalg.norm(cosmic.values, axis=0)
    best = (a.T @ b).max(axis=1)
    return {
        "reference": os.path.relpath(cosmic_profiles, REPO_ROOT),
        "reference_sha256": sha256(cosmic_profiles),
        "signatures_compared": int(cosmic.shape[1]),
        "artifact_signatures_included": True,
        "match_threshold": threshold,
        "matched": int((best >= threshold).sum()),
        "unmatched": int((best < threshold).sum()),
    }


def readable(record):
    pre, clu, inp, res = (record["preprocessing"], record["clustering"],
                          record["inputs"], record["results"])
    lines = [
        f"eSignatures clustering run — {record['mutation_type']}",
        f"Run name:    {record.get('run_name', '')}",
        f"Date:        {record['run_date']}",
        f"Code:        commit {record['code']['commit']}"
        + (" (with uncommitted changes)" if record['code']['uncommitted_changes'] else ""),
        "",
        "PREPROCESSING",
        f"  Minimum mutations per sample: {cutoff_config.describe(pre['min_mutations'] or {})}"
        f" ({'config default' if pre.get('min_mutations_setting', 'default') == 'default' else 'option: ' + pre['min_mutations_setting']})",
        f"  Excluded by name:             {pre['exclude'] or 'none'}",
        f"  Opportunity normalization:    {pre['normalization']}",
        f"  Config:                       {pre['config']}" + (" (skipped)" if pre['skipped'] else ""),
        "",
        "CLUSTERING",
        f"  Clustering cosine similarity: {clu['cosine_similarity']} (distance {1 - clu['cosine_similarity']:.4g}), {clu['linkage']} linkage",
        f"  Custom thresholds:            {clu['custom_thresholds'] or 'none'}",
        f"  Main cluster:                 >= {clu['main_cluster_min_samples']} samples",
        f"  Consensus profile:            {clu['averaging_method']}",
        "",
    ]
    if "cosmic" in record:
        c = record["cosmic"]
        lines += [
            "COSMIC COMPARISON",
            f"  Reference:                    {c['reference']} ({c['signatures_compared']} signatures, artifacts included)",
            f"  Match threshold:              max cosine similarity >= {c['match_threshold']}",
            "",
        ]
    lines += [
        "RESULTS",
        f"  Samples clustered:            {res['samples_clustered']}",
        f"  Main clusters (eSS):          {res['main_clusters']}",
        f"  Small clusters:               {res['small_clusters']}",
        f"  Singletons:                   {res['singletons']}",
    ]
    if "cosmic" in record:
        lines.append(f"  COSMIC matched / unmatched:   {record['cosmic']['matched']} / {record['cosmic']['unmatched']}")
    lines += ["", "INPUTS (sha256)"]
    for name, h in (inp.get("input_files_sha256") or {}).items():
        lines.append(f"  {name}: {h}")
    for name, h in (inp.get("context_tables_sha256") or {}).items():
        lines.append(f"  {name}: {h}")
    if "cosmic" in record:
        lines.append(f"  {os.path.basename(record['cosmic']['reference'])}: {record['cosmic']['reference_sha256']}")
    lines += ["", "SOFTWARE"]
    lines += [f"  {k}: {v}" for k, v in record["software"].items()]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    n = sub.add_parser("name")
    n.add_argument("--mutation_type", required=True)
    n.add_argument("--preprocessing_config", default=DEFAULT_PREPROCESSING)
    n.add_argument("--custom_thresholds", default="none")
    n.add_argument("--normalization", default="none")
    n.add_argument("--min_mutations", default=None)

    s = sub.add_parser("summary")
    s.add_argument("--output_dir", required=True, help="results/<run name>/<TYPE>")
    s.add_argument("--run_name", default="")
    s.add_argument("--cosmic_profiles", required=True)
    s.add_argument("--cosmic_threshold", type=float, required=True)

    args = parser.parse_args()
    if args.command == "name":
        print(run_name(args.mutation_type.upper(), args.preprocessing_config,
                       args.custom_thresholds, args.normalization, args.min_mutations))
        return

    json_path = os.path.join(args.output_dir, "run_parameters.json")
    with open(json_path) as f:
        record = json.load(f)
    record["run_name"] = args.run_name
    record["cosmic"] = cosmic_comparison(args.output_dir, args.cosmic_profiles, args.cosmic_threshold)
    with open(json_path, "w") as f:
        json.dump(record, f, indent=2)
    txt_path = os.path.join(args.output_dir, "run_parameters.txt")
    with open(txt_path, "w") as f:
        f.write(readable(record))
    print(f"Run parameters → {txt_path}")


if __name__ == "__main__":
    main()
