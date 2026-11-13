#!/usr/bin/env python3
"""
make_expected.py — register a new version of the SBS input data and write the
results the tests expect for it.

    python tests/make_expected.py --version 2026-11-13 \
        --description "Collaborator's new human file: +70 samples"

Fingerprints data/input/SBS (or --input_dir), runs the same preprocessing and
clustering as the tests, and writes:

    tests/expected/<version>/SBS_cluster_membership.tsv               (default)
    tests/expected/<version>/SBS_cluster_membership_aai_split.tsv     (aai-split)
    tests/expected/<version>/SBS_cluster_membership_per_species.tsv   (per-species)

and adds the version to tests/expected/data_versions.json with its file
fingerprints and expected numbers (profiles, main / small / singletons, COSMIC
matched / unmatched). The version is a date, YYYY-MM-DD: the day the inputs
changed.

These become what the tests check from then on, so review them before
committing: compare the numbers with a full pipeline run on the same inputs,
and diff the membership files against the previous version's.
"""

import argparse
import datetime
import os
import re
import shutil
import sys
import tempfile

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TESTS_DIR)

import data_versions as dv  # noqa: E402
import test_clustering_reproducibility as t  # noqa: E402  (same code as the tests)
from test_clustering_reproducibility import pc, get_config, cutoff_config  # noqa: E402


def cluster(data_dir, cfg, cutoffs, thresholds, tmp):
    counts, normalized, mappings = t.pc.load_data(
        data_dir, cfg, os.path.join(t.REPO_ROOT, "config", "sample_mapping.tsv"))
    pc.check_loaded_data(counts, normalized)
    pc.check_min_mutations(counts, mappings, cutoffs)
    counts = pc.filter_zero_columns(counts, tmp)
    normalized = normalized[counts.columns]
    result = t.run_clustering(normalized, tmp, dict(thresholds))
    return result, normalized, counts, mappings


def compute(input_dir, out_dir):
    """Expected numbers for the inputs in input_dir; membership files into out_dir."""
    cfg = get_config("SBS")
    pre_config = pc.load_preprocessing_config(t.PREPROCESSING_CONFIG)["SBS"]
    expected = {}
    tmp = tempfile.mkdtemp(prefix="ess_make_expected_")
    try:
        base = os.path.join(tmp, "input")
        shutil.copytree(input_dir, os.path.join(base, "SBS"))
        default_dir = pc.determine_data_directory(base, "SBS", t.PREPROCESSING_CONFIG)
        per_species_dir = pc.determine_data_directory(base, "SBS", t.PREPROCESSING_CONFIG,
                                                      min_mutations_setting="per-species")
        runs = [
            ("default", default_dir, cutoff_config.resolve(pre_config), cfg.default_custom_thresholds),
            ("aai_split", default_dir, cutoff_config.resolve(pre_config),
             cfg.custom_threshold_presets["aai-split"]),
            ("per_species", per_species_dir, cutoff_config.resolve(pre_config, "per-species"),
             cfg.default_custom_thresholds),
        ]
        for key, data_dir, cutoffs, thresholds in runs:
            result, normalized, counts, mappings = cluster(data_dir, cfg, cutoffs, thresholds, tmp)
            main, small, *_, singletons = result
            entry = {"clusters": [len(main), len(small), len(singletons)]}
            if key == "per_species":
                entry["profiles"] = int(normalized.shape[1])
            else:
                expected["profiles"] = int(normalized.shape[1])
                entry["cosmic"] = list(t.cosmic_counts(main, normalized, counts, mappings,
                                                       os.path.join(tmp, f"summary_{key}")))
            expected[key] = entry
            table = t.membership_table(main, small, singletons)
            table.to_csv(os.path.join(out_dir, t.MEMBERSHIP_FILES[key]), sep="\t", index=False)
            print(f"  {key:12s} {normalized.shape[1]} profiles, main/small/singletons "
                  f"{entry['clusters']}" + (f", COSMIC matched/unmatched {entry['cosmic']}"
                                            if "cosmic" in entry else ""))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return {"profiles": expected.pop("profiles"), **expected}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", required=True, help="date the inputs changed, YYYY-MM-DD")
    parser.add_argument("--description", required=True, help="what changed in these inputs")
    parser.add_argument("--input_dir", default=dv.INPUT_DIR,
                        help="unfiltered SBS96 files (default: data/input/SBS)")
    parser.add_argument("--force", action="store_true", help="replace an existing version")
    args = parser.parse_args()

    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", args.version):
        parser.error("--version must be a date, YYYY-MM-DD")
    datetime.date.fromisoformat(args.version)

    registry = dv.load_registry() if os.path.exists(dv.REGISTRY) else {}
    files = dv.describe(args.input_dir)
    for version, entry in registry.items():
        if entry["files"] == files and version != args.version:
            sys.exit(f"These inputs are already registered as {version}.")
    if args.version in registry and not args.force:
        sys.exit(f"{args.version} is already registered; use --force to replace it.")

    out_dir = dv.expected_dir(args.version)
    os.makedirs(out_dir, exist_ok=True)
    print(f"Expected results for input data version {args.version} ({args.input_dir})")
    expected = compute(args.input_dir, out_dir)
    registry[args.version] = {"description": args.description, "files": files,
                              "expected": expected}
    dv.save_registry(registry)
    print(f"Registered {args.version} in {dv.REGISTRY}; membership files in {out_dir}")
    print("Review the numbers and diff the membership files before committing.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
