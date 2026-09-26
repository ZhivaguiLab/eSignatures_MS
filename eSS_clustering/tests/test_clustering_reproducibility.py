#!/usr/bin/env python3
"""
Reproducibility checks for the SBS eSignature clustering.

Runs the clustering on data/input/SBS (no plotting of signature profiles) and
checks it against the published result frozen in
tests/expected/SBS_cluster_membership.tsv: 49 main clusters, 16 small
clusters, 131 singletons from 671 profiles, with identical membership and
eSS numbering.

Run from the eSS_clustering directory:

    python -m unittest discover tests -v
"""

import os
import shutil
import sys
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "pipeline"))

# pdf2image is only used to convert plots to PNG, which these tests skip.
try:
    import pdf2image  # noqa: F401
except ImportError:
    import types
    sys.modules["pdf2image"] = types.SimpleNamespace(convert_from_path=None)

import perform_clustering as pc  # noqa: E402
from utils.mutation_type import get_config  # noqa: E402

EXPECTED_MEMBERSHIP = os.path.join(REPO_ROOT, "tests", "expected",
                                   "SBS_cluster_membership.tsv")
COSMIC_SBS = os.path.join(REPO_ROOT, "data", "references",
                          "COSMIC_v3.6_SBS_GRCh38.txt")
COSINE_SIMILARITY = 0.9
COSMIC_MATCH_THRESHOLD = 0.85

# Samples removed by the SBS exclusion patterns in config/preprocessing.yaml.
EXCLUDED_MOUSE = [
    'MEF_Deoxynivalenol_Patulin__Genome_1', 'MEF_Deoxynivalenol_Patulin__Genome_2',
    'MEF_Deoxynivalenol_Patulin__Genome_3', 'MEF_Deoxynivalenol_Patulin__Genome_4',
    'MEF_Deoxynivalenol__Genome_1', 'MEF_Deoxynivalenol__Genome_2',
    'MEF_Xenon__Genome_1', 'MEF_Xenon__Genome_2',
    'MEF_Xenon_XPA-/-_Genome_1', 'MEF_Xenon_XPA-/-_Genome_2',
]


def run_clustering(normalized_df, output_dir, custom_thresholds):
    return pc.cluster_signatures_with_custom_thresholds(
        X=normalized_df,
        output_dir=output_dir,
        output_name="dendrogram",
        cos_dist_threshold=1 - COSINE_SIMILARITY,
        custom_thresholds=custom_thresholds,
    )


class ClusteringReproducibilityTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="ess_clustering_test_")
        cls.cfg = get_config("SBS")

        # Work on a copy so the preprocessing cache isn't written into data/.
        data_base = os.path.join(cls.tmp, "input")
        shutil.copytree(os.path.join(REPO_ROOT, "data", "input", "SBS"),
                        os.path.join(data_base, "SBS"))
        cls.data_base = data_base
        data_dir = pc.determine_data_directory(
            data_base, "SBS",
            os.path.join(REPO_ROOT, "config", "preprocessing.yaml"))

        counts, normalized, cls.mappings = pc.load_data(
            data_dir, cls.cfg,
            os.path.join(REPO_ROOT, "config", "sample_mapping.tsv"))
        pc.check_loaded_data(counts, normalized)
        cls.counts = pc.filter_zero_columns(counts, cls.tmp)
        cls.normalized = normalized[cls.counts.columns]

        cls.default = run_clustering(cls.normalized, cls.tmp,
                                     dict(cls.cfg.default_custom_thresholds))
        cls.no_custom = run_clustering(cls.normalized, cls.tmp, {})
        cls.expected = pd.read_csv(EXPECTED_MEMBERSHIP, sep="\t",
                                   keep_default_na=False)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # ── Input ────────────────────────────────────────────────────────────────

    def test_input_profile_count(self):
        self.assertEqual(self.normalized.shape, (96, 671))

    def test_mouse_exclusions_applied(self):
        mouse_originals = set(self.mappings["mouse"])
        self.assertFalse(mouse_originals & set(EXCLUDED_MOUSE))

    # ── Cluster counts ───────────────────────────────────────────────────────

    def test_default_cluster_counts(self):
        main, small, *_, singletons = self.default
        self.assertEqual((len(main), len(small), len(singletons)), (49, 16, 131))

    def test_cluster_counts_without_custom_thresholds(self):
        main, small, *_, singletons = self.no_custom
        self.assertEqual((len(main), len(small), len(singletons)), (48, 16, 131))

    def test_every_sample_assigned_once(self):
        main, small, *_, singletons = self.default
        assigned = ([s for v in main.values() for s in v] +
                    [s for v in small.values() for s in v] + singletons)
        self.assertEqual(len(assigned), len(set(assigned)))
        self.assertEqual(set(assigned), set(self.normalized.columns))

    # ── Membership and numbering ─────────────────────────────────────────────

    def test_membership_matches_published(self):
        main, small, *_, singletons = self.default
        rows = [(s, "main", f"eSS{i + 1}")
                for i, v in enumerate(main.values()) for s in v]
        rows += [(s, "small", f"small{i + 1}")
                 for i, v in enumerate(small.values()) for s in v]
        rows += [(s, "singleton", "") for s in singletons]
        got = pd.DataFrame(rows, columns=["sample", "group", "cluster"])

        merged = self.expected.merge(got, on="sample", how="outer",
                                     suffixes=("_expected", "_got"),
                                     indicator=True)
        diff = merged[(merged["_merge"] != "both") |
                      (merged["group_expected"] != merged["group_got"]) |
                      (merged["cluster_expected"] != merged["cluster_got"])]
        self.assertTrue(diff.empty,
                        f"{len(diff)} samples differ from the published "
                        f"clustering:\n{diff.head(20).to_string()}")

    def test_tree_cut_matches_dendrogram_colours(self):
        # Rebuild clusters the old way, from dendrogram leaf colours, and
        # check they match the fcluster tree cut used by the pipeline.
        X = self.normalized
        Z = sch.linkage(pdist(X.T, "cosine"), "average")
        R = sch.dendrogram(Z, color_threshold=1 - COSINE_SIMILARITY,
                           above_threshold_color="gray", no_plot=True)
        by_colour = {}
        for colour, leaf in zip(R["leaves_color_list"], R["leaves"]):
            by_colour.setdefault(colour, set()).add(X.columns[leaf])
        by_colour.pop("gray")

        main, small, *_ = self.no_custom
        pipeline = [set(v) for v in main.values()] + [set(v) for v in small.values()]
        self.assertEqual(sorted(map(sorted, by_colour.values())),
                         sorted(map(sorted, pipeline)))

    def test_custom_split_only_changes_aai_dbp_cluster(self):
        main_d, small_d, *_ = self.default
        main_n, small_n, *_ = self.no_custom
        default = {frozenset(v) for v in main_d.values()} | \
                  {frozenset(v) for v in small_d.values()}
        no_custom = {frozenset(v) for v in main_n.values()} | \
                    {frozenset(v) for v in small_n.values()}
        split_from = no_custom - default
        split_into = default - no_custom
        self.assertEqual(len(split_from), 1)
        self.assertEqual(len(split_into), 2)
        self.assertEqual(set().union(*split_into), set(next(iter(split_from))))

    # ── Consensus profiles and COSMIC matching ───────────────────────────────

    def test_cosmic_match_count(self):
        main, *_ = self.default
        out = os.path.join(self.tmp, "summary")
        os.makedirs(out, exist_ok=True)
        path = pc.summarize_and_average_clusters(
            main, self.normalized, self.counts, out, "main", self.mappings,
            "eSS", averaging_method="equal_replicate")
        profiles = pd.read_csv(path, sep="\t", index_col=0)
        cosmic = pd.read_csv(COSMIC_SBS, sep="\t", index_col=0)
        self.assertTrue(profiles.index.equals(cosmic.index))

        A = profiles.values / np.linalg.norm(profiles.values, axis=0)
        B = cosmic.values / np.linalg.norm(cosmic.values, axis=0)
        best = (A.T @ B).max(axis=1)
        self.assertEqual(int((best >= COSMIC_MATCH_THRESHOLD).sum()), 26)
        self.assertEqual(int((best < COSMIC_MATCH_THRESHOLD).sum()), 23)

    # ── Safety checks fail loudly ────────────────────────────────────────────

    def test_missing_cluster_sample_raises(self):
        main, *_ = self.default
        first = next(iter(main))
        broken = dict(main)
        broken[first] = main[first] + ["not_a_real_sample"]
        with self.assertRaises(AssertionError):
            pc.summarize_and_average_clusters(
                broken, self.normalized, self.counts, self.tmp, "main",
                self.mappings, "eSS")

    def test_mismatched_normalized_profile_raises(self):
        normalized = self.normalized.copy()
        a, b = normalized.columns[:2]
        normalized[[a, b]] = normalized[[b, a]].values
        with self.assertRaises(AssertionError):
            pc.check_loaded_data(self.counts, normalized)

    def test_sample_missing_from_normalized_raises(self):
        with self.assertRaises(AssertionError):
            pc.check_loaded_data(self.counts, self.normalized.iloc[:, 1:])

    def test_preprocessing_cache_rebuilt_when_patterns_change(self):
        config = os.path.join(self.tmp, "preprocessing.yaml")
        with open(config, "w") as f:
            f.write("SBS:\n  exclude:\n    - Xenon\ncase_sensitive: false\n")
        data_dir = pc.determine_data_directory(self.data_base, "SBS", config)
        mouse = pd.read_csv(os.path.join(data_dir, "filtered_mouse_307.txt"),
                            sep="\t", index_col=0)
        self.assertTrue(any("Deoxynivalenol" in c for c in mouse.columns))
        self.assertFalse(any("Xenon" in c for c in mouse.columns))


if __name__ == "__main__":
    unittest.main(verbosity=2)
