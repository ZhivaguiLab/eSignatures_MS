"""
data_versions.py — which release of the SBS input data is in data/input/SBS.

Each release of the input matrices is registered in
tests/expected/data_versions.json under its date (YYYY-MM-DD), with:

    description   what changed in that release
    files         for each unfiltered_*_SBS96.txt: number of samples and a
                  fingerprint of its contents (sample names, contexts, counts)
    expected      the results the tests expect for it (profile and cluster
                  counts, COSMIC matched/unmatched)

and its expected cluster membership files are in tests/expected/<date>/.

The fingerprint is a SHA-256 of the matrix itself, not of the file's bytes, so
line endings or number formatting don't matter, only the data.

New releases so far only add samples after the existing ones, so an earlier
release can be rebuilt from the current files by keeping each file's first N
samples (build_inputs); the fingerprint confirms the rebuild is exact.
"""

import hashlib
import json
import os

import numpy as np
import pandas as pd

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
EXPECTED_DIR = os.path.join(TESTS_DIR, "expected")
REGISTRY = os.path.join(EXPECTED_DIR, "data_versions.json")
INPUT_DIR = os.path.join(REPO_ROOT, "data", "input", "SBS")

# The release the published clustering on main (and the old filtered_*_307.txt
# inputs) was made from.
PUBLISHED_VERSION = "2026-09-28"


def read_matrix(path):
    return pd.read_csv(path, sep="\t", index_col=0)


def fingerprint(df):
    h = hashlib.sha256()
    h.update("\t".join(map(str, df.index)).encode())
    h.update(b"\n")
    h.update("\t".join(map(str, df.columns)).encode())
    h.update(b"\n")
    h.update(np.ascontiguousarray(df.values.astype(np.float64)).tobytes())
    return h.hexdigest()


def input_files(input_dir=INPUT_DIR):
    return sorted(f for f in os.listdir(input_dir) if f.endswith(".txt"))


def describe(input_dir=INPUT_DIR):
    """{file: {"samples": n, "fingerprint": sha256}} for the inputs in input_dir."""
    out = {}
    for f in input_files(input_dir):
        df = read_matrix(os.path.join(input_dir, f))
        out[f] = {"samples": int(df.shape[1]), "fingerprint": fingerprint(df)}
    return out


def load_registry():
    if not os.path.exists(REGISTRY):
        return {}
    with open(REGISTRY) as f:
        return json.load(f)


def save_registry(registry):
    with open(REGISTRY, "w") as f:
        json.dump(dict(sorted(registry.items())), f, indent=2)
        f.write("\n")


def current_version(input_dir=INPUT_DIR):
    """The registered release that data/input/SBS matches, or None."""
    files = describe(input_dir)
    for version, entry in load_registry().items():
        if entry["files"] == files:
            return version
    return None


def expected_dir(version):
    return os.path.join(EXPECTED_DIR, version)


def build_inputs(version, dest, input_dir=INPUT_DIR):
    """
    Write the inputs of an earlier (or the same) release into dest, from the
    current files, keeping each file's first N samples. Returns dest, or None
    if the current files can't reproduce that release exactly.
    """
    entry = load_registry()[version]
    os.makedirs(dest, exist_ok=True)
    for f, spec in entry["files"].items():
        src = os.path.join(input_dir, f)
        if not os.path.exists(src):
            return None
        df = read_matrix(src)
        if df.shape[1] < spec["samples"]:
            return None
        df = df.iloc[:, :spec["samples"]]
        if fingerprint(df) != spec["fingerprint"]:
            return None
        df.to_csv(os.path.join(dest, f), sep="\t")
    return dest


UNKNOWN_DATA_MESSAGE = (
    "The input data in data/input/SBS is not a registered data version "
    "(tests/expected/data_versions.json), so there are no expected results for it. "
    "If the new inputs are intended, register them and generate their expected results:\n"
    "    python tests/make_expected.py --version YYYY-MM-DD --description \"what changed\"\n"
    "then review the new numbers and membership files before committing them.")
