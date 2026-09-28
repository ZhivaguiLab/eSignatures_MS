"""
utils/opportunity_normalization.py

Optional trinucleotide-opportunity normalization of SBS96 profiles (manuscript
revision item A4). Two methods, both applied per species to a samples x
contexts count table, with the sequencing technology of each sample taken from
its name ("exome" -> WES, "genome" -> WGS):

wes-to-wgs
    Human and mouse WES samples are rescaled onto their own genome's WGS
    opportunity: counts[c] * WGS_freq[c] / WES_freq[c], then scaled back to
    the sample's original total. WGS samples, other species and samples of
    unknown technology are unchanged. Normalized profile = counts / total.

own-opportunity
    Every human, mouse and rat sample (WGS and WES) is normalized by its own
    genome + technology opportunity: rate[c] = counts[c] / freq[c], profile =
    rate / sum(rate). Counts are unchanged. Chicken, C. elegans and samples
    of unknown technology keep the plain profile counts / total.

Opportunity tables are the SigProfilerMatrixGenerator context_counts files in
data/references/context_distributions/. See NORMALIZATION_APPROACH.md.
"""

import os
import re

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_CONTEXT_DIR = os.path.join(REPO_ROOT, "data", "references", "context_distributions")

# Genome build used for each species' opportunity tables, per method.
WES_TO_WGS_BUILDS = {"human": "GRCh38", "mouse": "mm10"}
OWN_OPPORTUNITY_BUILDS = {"human": "GRCh38", "mouse": "mm10", "rat": "rn7"}

METHODS = ("wes-to-wgs", "own-opportunity")


def trinuc_rel_freq(context_dir, genome, exome):
    suffix = "_96_exome.csv" if exome else "_96.csv"
    path = os.path.join(context_dir, f"context_counts_{genome}{suffix}")
    df = pd.read_csv(path, index_col=0)
    if "Y" in df.columns:
        df = df.drop(columns=["Y"])
    totals = df.sum(axis=1)
    return totals / totals.sum()


def context_to_trinuc(label):
    m = re.match(r"^(.)\[(.)>(.)\](.)$", label)
    five, ref, _alt, three = m.groups()
    return five + ref + three


def per_context(rel_freq, contexts):
    """Map trinucleotide frequencies onto SBS96 context labels (by label, not position)."""
    return pd.Series({c: rel_freq[context_to_trinuc(c)] for c in contexts})


def tech_of(sample_name):
    if re.search("exome", sample_name, re.I):
        return "WES"
    if re.search("genome", sample_name, re.I):
        return "WGS"
    return None


def tables_used(method, species):
    """Opportunity table files a method reads for a species (for provenance)."""
    builds = WES_TO_WGS_BUILDS if method == "wes-to-wgs" else OWN_OPPORTUNITY_BUILDS
    if species not in builds:
        return []
    g = builds[species]
    return [f"context_counts_{g}_96.csv", f"context_counts_{g}_96_exome.csv"]


def wes_to_wgs(counts, species, context_dir=DEFAULT_CONTEXT_DIR):
    """Return (counts, normalized, summary) with WES samples moved onto the WGS basis."""
    counts = counts.astype(float)
    if species not in WES_TO_WGS_BUILDS:
        return counts, counts.div(counts.sum(axis=0), axis=1), f"{species}: unchanged (no correction)"

    contexts = list(counts.index)
    genome = WES_TO_WGS_BUILDS[species]
    wgs = trinuc_rel_freq(context_dir, genome, exome=False)
    wes = trinuc_rel_freq(context_dir, genome, exome=True)
    ratio = per_context(wgs / wes, contexts)

    corrected = counts.copy()
    n = {"WES": 0, "WGS": 0, None: 0}
    for col in counts.columns:
        tech = tech_of(col)
        n[tech] += 1
        if tech == "WES":
            total = counts[col].sum()
            rescaled = counts[col] * ratio.reindex(contexts).values
            corrected[col] = rescaled / rescaled.sum() * total   # reshape only; keep total
    normalized = corrected.div(corrected.sum(axis=0), axis=1)
    summary = (f"{species}: {n['WES']} WES corrected to {genome} WGS basis, "
               f"{n['WGS']} WGS unchanged, {n[None]} unknown technology unchanged")
    return corrected, normalized, summary


def own_opportunity(counts, species, context_dir=DEFAULT_CONTEXT_DIR):
    """Return (counts, normalized, summary); counts unchanged, profiles opportunity-normalized."""
    counts = counts.astype(float)
    if species not in OWN_OPPORTUNITY_BUILDS:
        return counts, counts.div(counts.sum(axis=0), axis=1), \
            f"{species}: plain proportions (no opportunity table)"

    contexts = list(counts.index)
    genome = OWN_OPPORTUNITY_BUILDS[species]
    opp = {"WGS": per_context(trinuc_rel_freq(context_dir, genome, exome=False), contexts),
           "WES": per_context(trinuc_rel_freq(context_dir, genome, exome=True), contexts)}
    cols = {}
    n = {"WES": 0, "WGS": 0, None: 0}
    for col in counts.columns:
        tech = tech_of(col)
        n[tech] += 1
        if tech in opp:
            rate = counts[col] / opp[tech].reindex(contexts).values
            cols[col] = rate / rate.sum()
        else:
            cols[col] = counts[col] / counts[col].sum()
    normalized = pd.DataFrame(cols, index=contexts)[list(counts.columns)]
    summary = (f"{species}: {n['WGS']} WGS ({genome} genome opportunity), "
               f"{n['WES']} WES ({genome} exome opportunity), "
               f"{n[None]} unknown technology (plain proportions)")
    return counts, normalized, summary


SPECIES_PATTERNS = {
    'mouse':    re.compile(r'mouse', re.I),
    'rat':      re.compile(r'(?<![a-z])rat(?![a-z])', re.I),
    'chicken':  re.compile(r'chicken', re.I),
    'human':    re.compile(r'human', re.I),
    'celegans': re.compile(r'c[._]?elegans?', re.I),
}   # same rules as perform_clustering.SPECIES_PATTERNS


def normalize_directory(method, input_dir, output_dir, context_dir=DEFAULT_CONTEXT_DIR):
    """
    Apply a method to every species count file in input_dir (a preprocessed
    folder such as data/input_cleaned/SBS) and write the counts and
    normalized_<name>.tsv profiles to output_dir, in the layout the pipeline
    reads. Normalized files already in input_dir are replaced.
    """
    os.makedirs(output_dir, exist_ok=True)
    files = sorted(f for f in os.listdir(input_dir)
                   if f.endswith((".txt", ".tsv")) and "normaliz" not in f.lower())
    if not files:
        raise FileNotFoundError(f"No count files in {input_dir}")
    for name in files:
        species = next((k for k, pat in SPECIES_PATTERNS.items() if pat.search(name)), None)
        if species is None:
            print(f"  Skipping {name}: species not recognised")
            continue
        counts = pd.read_csv(os.path.join(input_dir, name), sep="\t", index_col=0)
        counts, normalized, summary = normalize(method, counts, species, context_dir)
        counts.to_csv(os.path.join(output_dir, name), sep="\t")
        normalized.to_csv(os.path.join(output_dir, f"normalized_{os.path.splitext(name)[0]}.tsv"),
                          sep="\t")
        print(f"  {summary}")


def normalize(method, counts, species, context_dir=DEFAULT_CONTEXT_DIR):
    if method == "wes-to-wgs":
        return wes_to_wgs(counts, species, context_dir)
    if method == "own-opportunity":
        return own_opportunity(counts, species, context_dir)
    raise ValueError(f"Unknown normalization {method!r}; choose from {METHODS}.")
