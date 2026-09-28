#!/usr/bin/env python3
"""
Reproducibility checks for the SBS eSignature clustering.

Runs preprocessing and clustering on data/input/SBS (no plotting of signature
profiles) and checks them against:

- tests/expected/SBS_cluster_membership.tsv: 48 main clusters, 16 small
  clusters, 123 singletons from 653 profiles (min_mutations 307 for every
  species, no custom thresholds), with identical membership and eSS numbering;
- tests/expected/SBS_cluster_membership_aai_split.tsv: the same with the
  optional 'aai-split' custom thresholds (49 / 16 / 123);
- tests/expected/published_307_input_samples.tsv: the samples in the
  previously published filtered_*_307.txt inputs, which preprocessing must
  reproduce when every species' cutoff is 307 (apart from the 18 C. elegans
  CX-5461 samples below 307 that the old filtering script let through);
- tests/expected/main_SBS_cluster_membership.tsv and main_SBS_eSS_profiles.tsv:
  the published clustering on main (671 profiles, 49 / 16 / 131). With 307 for
  every species, preprocessing plus the 18 samples the old script never tested
  must reproduce it exactly; without them, it must equal main minus those 18.

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
import run_info  # noqa: E402
from utils.mutation_type import get_config  # noqa: E402

EXPECTED_MEMBERSHIP = os.path.join(REPO_ROOT, "tests", "expected",
                                   "SBS_cluster_membership.tsv")
EXPECTED_MEMBERSHIP_AAI_SPLIT = os.path.join(REPO_ROOT, "tests", "expected",
                                             "SBS_cluster_membership_aai_split.tsv")
PUBLISHED_307_INPUTS = os.path.join(REPO_ROOT, "tests", "expected",
                                    "published_307_input_samples.tsv")
PREPROCESSING_CONFIG = os.path.join(REPO_ROOT, "config", "preprocessing.yaml")
MAIN_MEMBERSHIP = os.path.join(REPO_ROOT, "tests", "expected",
                               "main_SBS_cluster_membership.tsv")
MAIN_PROFILES = os.path.join(REPO_ROOT, "tests", "expected",
                             "main_SBS_eSS_profiles.tsv")
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
        cls.raw_dir = os.path.join(data_base, "SBS")
        cls.pre_config = pc.load_preprocessing_config(PREPROCESSING_CONFIG)["SBS"]
        data_dir = pc.determine_data_directory(data_base, "SBS", PREPROCESSING_CONFIG)
        cls.data_dir = data_dir

        counts, normalized, cls.mappings = pc.load_data(
            data_dir, cls.cfg,
            os.path.join(REPO_ROOT, "config", "sample_mapping.tsv"))
        pc.check_loaded_data(counts, normalized)
        pc.check_min_mutations(counts, cls.mappings, cls.pre_config["min_mutations"])
        cls.counts = pc.filter_zero_columns(counts, cls.tmp)
        cls.normalized = normalized[cls.counts.columns]

        cls.default = run_clustering(cls.normalized, cls.tmp,
                                     dict(cls.cfg.default_custom_thresholds))
        cls.aai_split = run_clustering(cls.normalized, cls.tmp,
                                       dict(cls.cfg.custom_threshold_presets["aai-split"]))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    # ── Input ────────────────────────────────────────────────────────────────

    def test_input_profile_count(self):
        self.assertEqual(self.normalized.shape, (96, 653))

    def test_mouse_exclusions_applied(self):
        mouse_originals = set(self.mappings["mouse"])
        self.assertFalse(mouse_originals & set(EXCLUDED_MOUSE))

    # ── Cluster counts ───────────────────────────────────────────────────────

    def test_default_cluster_counts(self):
        main, small, *_, singletons = self.default
        self.assertEqual((len(main), len(small), len(singletons)), (48, 16, 123))

    def test_cluster_counts_with_aai_split(self):
        main, small, *_, singletons = self.aai_split
        self.assertEqual((len(main), len(small), len(singletons)), (49, 16, 123))

    def test_no_custom_thresholds_by_default(self):
        self.assertEqual(self.cfg.default_custom_thresholds, {})
        self.assertEqual(pc.parse_custom_thresholds(None, self.cfg), {})
        self.assertEqual(pc.parse_custom_thresholds("aai-split", self.cfg),
                         {"Aristolochic_acid_I": 0.095, "Dibenzo[a,l]pyrene": 0.095})
        with self.assertRaises(ValueError):
            pc.parse_custom_thresholds("aai_split", self.cfg)

    def test_every_sample_assigned_once(self):
        main, small, *_, singletons = self.default
        assigned = ([s for v in main.values() for s in v] +
                    [s for v in small.values() for s in v] + singletons)
        self.assertEqual(len(assigned), len(set(assigned)))
        self.assertEqual(set(assigned), set(self.normalized.columns))

    # ── Membership and numbering ─────────────────────────────────────────────

    def _assert_membership(self, result, expected_path):
        main, small, *_, singletons = result
        got = membership_table(main, small, singletons)
        expected = pd.read_csv(expected_path, sep="\t", keep_default_na=False)
        merged = expected.merge(got, on="sample", how="outer",
                                     suffixes=("_expected", "_got"),
                                     indicator=True)
        diff = merged[(merged["_merge"] != "both") |
                      (merged["group_expected"] != merged["group_got"]) |
                      (merged["cluster_expected"] != merged["cluster_got"])]
        self.assertTrue(diff.empty,
                        f"{len(diff)} samples differ from {os.path.basename(expected_path)}:"
                        f"\n{diff.head(20).to_string()}")

    def test_membership_matches_expected(self):
        self._assert_membership(self.default, EXPECTED_MEMBERSHIP)

    def test_membership_matches_expected_with_aai_split(self):
        self._assert_membership(self.aai_split, EXPECTED_MEMBERSHIP_AAI_SPLIT)

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

        main, small, *_ = self.default
        pipeline = [set(v) for v in main.values()] + [set(v) for v in small.values()]
        self.assertEqual(sorted(map(sorted, by_colour.values())),
                         sorted(map(sorted, pipeline)))

    def test_aai_split_only_changes_aai_dbp_cluster(self):
        main_s, small_s, *_ = self.aai_split
        main_d, small_d, *_ = self.default
        split = {frozenset(v) for v in main_s.values()} | \
                {frozenset(v) for v in small_s.values()}
        default = {frozenset(v) for v in main_d.values()} | \
                  {frozenset(v) for v in small_d.values()}
        split_from = default - split
        split_into = split - default
        self.assertEqual(len(split_from), 1)
        self.assertEqual(len(split_into), 2)
        self.assertEqual(set().union(*split_into), set(next(iter(split_from))))

    # ── Consensus profiles and COSMIC matching ───────────────────────────────

    def _cosmic_counts(self, result, label):
        main, *_ = result
        out = os.path.join(self.tmp, f"summary_{label}")
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
        return int((best >= COSMIC_MATCH_THRESHOLD).sum()), int((best < COSMIC_MATCH_THRESHOLD).sum())

    def test_cosmic_match_count(self):
        self.assertEqual(self._cosmic_counts(self.default, "default"), (25, 23))

    def test_cosmic_match_count_with_aai_split(self):
        self.assertEqual(self._cosmic_counts(self.aai_split, "aai_split"), (26, 23))

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
        # Own copy of the data, so the shared cache used by other tests is
        # left alone. Build it with the real config, then change the config.
        base = os.path.join(self.tmp, "input_rebuild")
        shutil.copytree(self.raw_dir, os.path.join(base, "SBS"))
        pc.determine_data_directory(base, "SBS", PREPROCESSING_CONFIG)
        config = os.path.join(self.tmp, "preprocessing_xenon_only.yaml")
        with open(config, "w") as f:
            f.write("SBS:\n  exclude:\n    - Xenon\ncase_sensitive: false\n")
        data_dir = pc.determine_data_directory(base, "SBS", config)
        mouse = pd.read_csv(os.path.join(data_dir, "mouse_SBS96.txt"),
                            sep="\t", index_col=0)
        self.assertTrue(any("Deoxynivalenol" in c for c in mouse.columns))
        self.assertFalse(any("Xenon" in c for c in mouse.columns))
        # No min_mutations in this config, so low-count samples come back.
        self.assertLess(mouse.sum().min(), self.pre_config["min_mutations"]["mouse"])

    # ── Preprocessing: per-species minimum mutation counts ──────────────────

    def _raw_files(self):
        return {pc.detect_species(f)[0]: os.path.join(self.raw_dir, f)
                for f in sorted(os.listdir(self.raw_dir)) if f.endswith(".txt")}

    def test_every_sample_meets_species_cutoff(self):
        cutoffs = self.pre_config["min_mutations"]
        totals = self.counts.sum(axis=0)
        for species, mapping in self.mappings.items():
            names = [n for n in mapping.values() if n in totals]
            self.assertTrue(names, species)
            self.assertGreaterEqual(totals[names].min(), cutoffs[species], species)

    def test_preprocessing_matches_independent_filter(self):
        # Re-filter each unfiltered file from scratch and compare with the
        # cleaned files the pipeline clustered.
        cutoffs = self.pre_config["min_mutations"]
        patterns = [p.lower() for p in self.pre_config["exclude"]]
        for species, path in self._raw_files().items():
            raw = pd.read_csv(path, sep="\t", index_col=0)
            totals = raw.sum(axis=0)
            keep = [c for c in raw.columns
                    if totals[c] >= cutoffs[species]
                    and not any(p in c.lower() for p in patterns)]
            name = pc.cleaned_file_name(os.path.basename(path))
            cleaned = pd.read_csv(os.path.join(self.data_dir, name), sep="\t", index_col=0)
            self.assertEqual(list(cleaned.columns), keep, species)
            np.testing.assert_array_equal(cleaned.values, raw[keep].values)

            norm = pd.read_csv(os.path.join(self.data_dir, "normalized_" + name[:-4] + ".tsv"),
                               sep="\t", index_col=0)
            self.assertEqual(list(norm.columns), keep, species)
            np.testing.assert_allclose(norm.values, (raw[keep] / raw[keep].sum()).values,
                                       rtol=0, atol=1e-12)

    def test_removed_log_accounts_for_every_sample(self):
        log = pd.read_csv(os.path.join(self.data_dir, pc.PREPROCESSING_REMOVED_LOG))
        cutoffs = self.pre_config["min_mutations"]
        for species, path in self._raw_files().items():
            raw = pd.read_csv(path, sep="\t", index_col=0)
            removed = log[log.species == species]
            kept = pd.read_csv(os.path.join(self.data_dir,
                                            pc.cleaned_file_name(os.path.basename(path))),
                               sep="\t", index_col=0).columns
            self.assertEqual(sorted(set(removed["sample"]) | set(kept)), sorted(raw.columns))
            self.assertFalse(set(removed["sample"]) & set(kept))
            below = removed[removed.reason.str.startswith("below")]
            self.assertTrue((below.total_mutations < cutoffs[species]).all(), species)

    def test_cx5461_celegans_filtered(self):
        # The old filtering script only tested names containing "exome" or
        # "genome", so these 21 samples were never filtered. Only the three
        # with >= 307 SBSs may remain.
        clustered = set(self.normalized.columns)
        self.assertEqual(sorted(s for s in clustered if s.startswith("celegans_CX-5461")),
                         ["celegans_CX-5461_UVA14", "celegans_CX-5461_UVA17",
                          "celegans_CX-5461_UVA18"])

    def test_307_cutoff_reproduces_published_inputs(self):
        # With 307 for every species, preprocessing must keep exactly the
        # samples in the previously published filtered_*_307.txt files,
        # except the 18 C. elegans CX-5461 samples below 307.
        config = os.path.join(self.tmp, "preprocessing_307.yaml")
        with open(config, "w") as f:
            f.write("SBS:\n  min_mutations:\n" +
                    "".join(f"    {sp}: 307\n" for sp in self._raw_files()) +
                    "case_sensitive: false\n")
        base = os.path.join(self.tmp, "input307")
        shutil.copytree(self.raw_dir, os.path.join(base, "SBS"))
        data_dir = pc.determine_data_directory(base, "SBS", config)

        published = pd.read_csv(PUBLISHED_307_INPUTS, sep="\t")
        for species, path in self._raw_files().items():
            cleaned = pd.read_csv(os.path.join(data_dir, pc.cleaned_file_name(os.path.basename(path))),
                                  sep="\t", index_col=0)
            pub = published[published.species == species].set_index("sample").total_mutations
            extra_published = set(pub.index) - set(cleaned.columns)
            self.assertFalse(set(cleaned.columns) - set(pub.index), species)
            if species == "celegans":
                self.assertEqual(len(extra_published), 18)
                self.assertTrue(all(s.startswith("CX-5461") for s in extra_published))
                self.assertTrue((pub[sorted(extra_published)] < 307).all())
            else:
                self.assertFalse(extra_published, species)
            np.testing.assert_allclose(cleaned.sum()[cleaned.columns].values,
                                       pub[cleaned.columns].values, rtol=0, atol=1e-5)

    def test_cached_preprocessing_is_identical(self):
        before = pc.preprocessing_stamp(self.data_dir, [], False)
        again = pc.determine_data_directory(self.data_base, "SBS", PREPROCESSING_CONFIG)
        self.assertEqual(again, self.data_dir)
        self.assertEqual(pc.preprocessing_stamp(self.data_dir, [], False), before)

    def test_missing_species_cutoff_raises(self):
        config = os.path.join(self.tmp, "preprocessing_missing.yaml")
        with open(config, "w") as f:
            f.write("SBS:\n  min_mutations:\n    mouse: 235\ncase_sensitive: false\n")
        base = os.path.join(self.tmp, "input_missing")
        shutil.copytree(self.raw_dir, os.path.join(base, "SBS"))
        with self.assertRaises(ValueError):
            pc.determine_data_directory(base, "SBS", config)

    def test_sample_below_cutoff_raises(self):
        with self.assertRaises(AssertionError):
            pc.check_min_mutations(self.counts, self.mappings,
                                   {sp: 10**9 for sp in self.mappings})



def membership_table(main, small, singletons):
    rows = [(s, "main", f"eSS{i + 1}") for i, v in enumerate(main.values()) for s in v]
    rows += [(s, "small", f"small{i + 1}") for i, v in enumerate(small.values()) for s in v]
    rows += [(s, "singleton", "") for s in singletons]
    return pd.DataFrame(rows, columns=["sample", "group", "cluster"])


def old_notebook_tested(name):
    """filter_based_on_tmb.ipynb only applied its cutoff to these names."""
    return "exome" in name.lower() or "genome" in name.lower()


class MatchesMainAt307Test(unittest.TestCase):
    """
    Preprocessing with 307 for every species against the published
    clustering on main, which was built from the old notebook's
    filtered_*_307.txt files with the AAI/DBP split (so these tests use the
    'aai-split' preset).
    """

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="ess_clustering_307_")
        cls.cfg = get_config("SBS")
        raw_dir = os.path.join(REPO_ROOT, "data", "input", "SBS")
        base = os.path.join(cls.tmp, "input")
        shutil.copytree(raw_dir, os.path.join(base, "SBS"))
        cls.raw_files = {pc.detect_species(f)[0]: os.path.join(raw_dir, f)
                         for f in sorted(os.listdir(raw_dir)) if f.endswith(".txt")}

        # Same exclusions as the real config, 307 for every species.
        cls.exclude = pc.load_preprocessing_config(PREPROCESSING_CONFIG)["SBS"]["exclude"]
        config = os.path.join(cls.tmp, "preprocessing_307.yaml")
        with open(config, "w") as f:
            f.write("SBS:\n  exclude:\n" + "".join(f"    - {p}\n" for p in cls.exclude) +
                    "  min_mutations:\n" + "".join(f"    {sp}: 307\n" for sp in cls.raw_files) +
                    "case_sensitive: false\n")
        cls.pure_dir = pc.determine_data_directory(base, "SBS", config)

        # Recreate the old notebook's gap: add back the samples it never
        # tested (names without "exome"/"genome"), in their original order.
        cls.legacy_dir = os.path.join(cls.tmp, "legacy", "SBS")
        os.makedirs(cls.legacy_dir)
        cls.added_back = []
        for species, path in cls.raw_files.items():
            raw = pd.read_csv(path, sep="\t", index_col=0)
            name = pc.cleaned_file_name(os.path.basename(path))
            kept = set(pd.read_csv(os.path.join(cls.pure_dir, name), sep="\t", index_col=0).columns)
            legacy = [c for c in raw.columns
                      if c in kept or (not old_notebook_tested(c)
                                       and not any(p.lower() in c.lower() for p in cls.exclude))]
            cls.added_back += [(species, c) for c in legacy if c not in kept]
            counts = raw[legacy]
            counts.to_csv(os.path.join(cls.legacy_dir, name), sep="\t")
            counts.div(counts.sum(axis=0), axis=1).to_csv(
                os.path.join(cls.legacy_dir, "normalized_" + name[:-4] + ".tsv"), sep="\t")

        cls.pure = cls._cluster(cls.pure_dir, "pure")
        cls.legacy = cls._cluster(cls.legacy_dir, "legacy")
        cls.main_membership = pd.read_csv(MAIN_MEMBERSHIP, sep="\t", keep_default_na=False)
        cls.main_profiles = pd.read_csv(MAIN_PROFILES, sep="\t", index_col=0)

    @classmethod
    def _cluster(cls, data_dir, label):
        counts, normalized, mappings = pc.load_data(
            data_dir, cls.cfg, os.path.join(REPO_ROOT, "config", "sample_mapping.tsv"))
        pc.check_loaded_data(counts, normalized)
        counts = pc.filter_zero_columns(counts, cls.tmp)
        normalized = normalized[counts.columns]
        main, small, *_, singletons = run_clustering(
            normalized, cls.tmp, dict(cls.cfg.custom_threshold_presets["aai-split"]))
        out = os.path.join(cls.tmp, f"summary_{label}")
        os.makedirs(out)
        profiles = pd.read_csv(pc.summarize_and_average_clusters(
            main, normalized, counts, out, "main", mappings, "eSS"), sep="\t", index_col=0)
        return dict(counts=counts, main=main, small=small, singletons=singletons,
                    profiles=profiles, table=membership_table(main, small, singletons))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_only_cx5461_celegans_were_untested(self):
        self.assertEqual(len(self.added_back), 18)
        self.assertTrue(all(sp == "celegans" and s.startswith("CX-5461")
                            for sp, s in self.added_back))

    def test_legacy_inputs_match_published_inputs(self):
        # Same samples, in the same order, as main's filtered_*_307.txt files
        # (after main's Xenon/Deoxynivalenol exclusions).
        published = pd.read_csv(PUBLISHED_307_INPUTS, sep="\t")
        for species, path in self.raw_files.items():
            name = pc.cleaned_file_name(os.path.basename(path))
            legacy = pd.read_csv(os.path.join(self.legacy_dir, name), sep="\t", index_col=0)
            pub = published[published.species == species]
            pub = pub[~pub["sample"].str.lower().str.contains("|".join(p.lower() for p in self.exclude))]
            self.assertEqual(list(legacy.columns), list(pub["sample"]), species)
            np.testing.assert_allclose(legacy.sum().values, pub.total_mutations.values,
                                       rtol=0, atol=1e-5)

    def test_legacy_reproduces_main_exactly(self):
        t = self.legacy["table"]
        self.assertEqual(len(t), 671)
        self.assertEqual((len(self.legacy["main"]), len(self.legacy["small"]),
                          len(self.legacy["singletons"])), (49, 16, 131))
        pd.testing.assert_frame_equal(t.reset_index(drop=True),
                                      self.main_membership.reset_index(drop=True))

    def test_legacy_reproduces_main_eSS_profiles(self):
        got, expected = self.legacy["profiles"], self.main_profiles
        self.assertEqual(list(got.columns), list(expected.columns))
        self.assertTrue(got.index.equals(expected.index))
        np.testing.assert_allclose(got.values, expected.values, rtol=0, atol=1e-15)

    def test_pure_307_equals_main_minus_untested_samples(self):
        removed = {f"celegans_{s}" for _, s in self.added_back}
        main_t = self.main_membership
        self.assertEqual(set(self.pure["table"]["sample"]), set(main_t["sample"]) - removed)

        def groups(t, group):
            return {frozenset(g["sample"]) for _, g in t[t.group == group].groupby("cluster")}
        main_main = {frozenset(c - removed) for c in groups(main_t, "main")}
        self.assertEqual(main_main, groups(self.pure["table"], "main"))
        self.assertEqual(groups(main_t, "small"), groups(self.pure["table"], "small"))
        self.assertEqual(set(self.pure["singletons"]),
                         set(main_t[main_t.group == "singleton"]["sample"]) - removed)
        self.assertEqual((len(self.pure["main"]), len(self.pure["small"]),
                          len(self.pure["singletons"])), (49, 16, 123))

    def test_pure_307_unchanged_eSS_have_identical_profiles(self):
        removed = {f"celegans_{s}" for _, s in self.added_back}
        main_t = self.main_membership[self.main_membership.group == "main"]
        pure_ids = {frozenset(v): f"eSS{i + 1}" for i, v in enumerate(self.pure["main"].values())}
        unchanged = 0
        for cid, g in main_t.groupby("cluster"):
            members = frozenset(g["sample"])
            if members & removed:
                continue
            new_id = pure_ids[members]
            np.testing.assert_allclose(self.pure["profiles"][new_id].values,
                                       self.main_profiles[cid].values, rtol=0, atol=1e-15)
            unchanged += 1
        self.assertEqual(unchanged, 48)

    def test_pure_307_cosmic_match_count(self):
        cosmic = pd.read_csv(COSMIC_SBS, sep="\t", index_col=0)
        p = self.pure["profiles"].reindex(cosmic.index)
        best = ((p.values / np.linalg.norm(p.values, axis=0)).T @
                (cosmic.values / np.linalg.norm(cosmic.values, axis=0))).max(axis=1)
        self.assertEqual((int((best >= COSMIC_MATCH_THRESHOLD).sum()),
                          int((best < COSMIC_MATCH_THRESHOLD).sum())), (26, 23))



CONTEXT_DIR = os.path.join(REPO_ROOT, "data", "references", "context_distributions")


def opportunity(genome, exome, contexts):
    """Independent re-implementation: per-context relative trinucleotide frequency."""
    df = pd.read_csv(os.path.join(CONTEXT_DIR, f"context_counts_{genome}_96{'_exome' if exome else ''}.csv"),
                     index_col=0).drop(columns=["Y"], errors="ignore")
    freq = df.sum(axis=1) / df.sum(axis=1).sum()
    return np.array([freq[c[0] + c[2] + c[6]] for c in contexts])   # "A[C>A]G" -> "ACG"


class OpportunityNormalizationTest(unittest.TestCase):
    """The optional wes-to-wgs / own-opportunity normalization in preprocessing."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="ess_clustering_norm_")
        base = os.path.join(cls.tmp, "input")
        shutil.copytree(os.path.join(REPO_ROOT, "data", "input", "SBS"), os.path.join(base, "SBS"))
        cls.plain = pc.determine_data_directory(base, "SBS", PREPROCESSING_CONFIG)
        cls.dirs = {m: pc.determine_data_directory(base, "SBS", PREPROCESSING_CONFIG, normalization=m)
                    for m in ("wes-to-wgs", "own-opportunity")}
        cls.species = ["celegans", "chicken", "human", "mouse", "rat"]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def _read(self, d, sp, normalized=False):
        name = f"normalized_{sp}_SBS96.tsv" if normalized else f"{sp}_SBS96.txt"
        return pd.read_csv(os.path.join(d, name), sep="\t", index_col=0)

    def test_separate_caches_same_samples(self):
        self.assertEqual(len({self.plain, *self.dirs.values()}), 3)
        for d in self.dirs.values():
            for sp in self.species:
                self.assertEqual(list(self._read(d, sp).columns), list(self._read(self.plain, sp).columns))

    def test_wes_to_wgs_matches_independent_calculation(self):
        d = self.dirs["wes-to-wgs"]
        builds = {"human": "GRCh38", "mouse": "mm10"}
        for sp in self.species:
            before, after = self._read(self.plain, sp), self._read(d, sp)
            np.testing.assert_allclose(after.sum().values, before.sum().values, rtol=1e-12)
            for col in before.columns:
                if sp in builds and "exome" in col.lower():
                    ratio = (opportunity(builds[sp], False, before.index) /
                             opportunity(builds[sp], True, before.index))
                    expected = before[col].values * ratio
                    expected = expected / expected.sum() * before[col].sum()
                else:
                    expected = before[col].values
                np.testing.assert_allclose(after[col].values, expected, rtol=1e-10, err_msg=f"{sp} {col}")
            norm = self._read(d, sp, normalized=True)
            np.testing.assert_allclose(norm.values, (after / after.sum()).values, atol=1e-15)

    def test_own_opportunity_matches_independent_calculation(self):
        d = self.dirs["own-opportunity"]
        builds = {"human": "GRCh38", "mouse": "mm10", "rat": "rn7"}
        for sp in self.species:
            before = self._read(self.plain, sp)
            np.testing.assert_array_equal(self._read(d, sp).values, before.values)   # counts unchanged
            norm = self._read(d, sp, normalized=True)
            for col in before.columns:
                tech = "WES" if "exome" in col.lower() else "WGS" if "genome" in col.lower() else None
                if sp in builds and tech:
                    rate = before[col].values / opportunity(builds[sp], tech == "WES", before.index)
                    expected = rate / rate.sum()
                else:
                    expected = before[col].values / before[col].sum()
                np.testing.assert_allclose(norm[col].values, expected, rtol=1e-10, err_msg=f"{sp} {col}")

    def test_uncorrected_species_unchanged(self):
        for d in self.dirs.values():
            for sp in ["celegans", "chicken"]:
                np.testing.assert_allclose(self._read(d, sp, True).values,
                                           self._read(self.plain, sp, True).values, atol=1e-15)

    def test_own_opportunity_profiles_pass_profile_check(self):
        d = self.dirs["own-opportunity"]
        counts, norm, _ = pc.load_data(d, get_config("SBS"),
                                       os.path.join(REPO_ROOT, "config", "sample_mapping.tsv"))
        with self.assertRaises(AssertionError):
            pc.check_loaded_data(counts, norm)                 # not counts / total
        pc.check_loaded_data(counts, norm, profiles_are_proportions=False)

    def test_normalization_not_allowed_for_dbs(self):
        with self.assertRaises(ValueError):
            pc.determine_data_directory(os.path.join(self.tmp, "input"), "DBS",
                                        PREPROCESSING_CONFIG, normalization="wes-to-wgs")


class RunNameTest(unittest.TestCase):

    def test_run_names(self):
        name = lambda **k: run_info.run_name("SBS", 0.9, PREPROCESSING_CONFIG, **k)
        self.assertEqual(name(), "min307_cos0.90")
        self.assertEqual(name(custom_thresholds="aai-split"), "min307_cos0.90_aai-split")
        self.assertEqual(name(normalization="wes-to-wgs"), "min307_cos0.90_wes-to-wgs")
        self.assertEqual(name(normalization="own-opportunity"), "min307_cos0.90_own-opportunity")


if __name__ == "__main__":
    unittest.main(verbosity=2)
