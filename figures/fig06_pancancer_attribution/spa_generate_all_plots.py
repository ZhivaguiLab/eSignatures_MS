#!/usr/bin/env python3
"""
Shared pan-cancer analysis library: TMB, smoking, and Fisher-prevalence plotting
for the sample-level SPA decomposition output. Imported by the run2_* scripts.

Input: the SPE de-novo hybrid, penalized run2 decomposition (per-cancer folders,
pooled). Each run writes analysis/ + analysis/fisher/.

Per-cancer-type results are rendered as volcanoes (log2 OR versus -log10 q).

Conventions: cohorts PCAWG + TCGA + Mutographs are all whole-genome, so MB = 2800
for every cohort. Each eSS is labelled by its COSMIC match at cosine >= 0.90, else
"novel".
"""
import os
import re
import sys
import glob
import numpy as np
import pandas as pd
from scipy.stats import ranksums, fisher_exact, pointbiserialr, mannwhitneyu
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import warnings
warnings.simplefilter("ignore")

BASE = os.environ.get("ESS_BASE", "../data/pancan_eSS_assignment")
sys.path.insert(0, os.path.join(BASE, "SPE_denovo_decomposition_eSS", "Smoking_Association"))
from ess_labels import load_labels

META_FILE = os.path.join(BASE, "HD Signatures WGS Datasets v2.0 - final_sample_summary_v2.tsv")
COMP = os.path.join(BASE, "SPE_denovo_decomposition_eSS", "Smoking_Association",
                    "eSS_COSMIC_composition.csv")
# --- cohort subset (env-overridable, defaults preserve the canonical 3-cohort run) ---
# ESS_COHORTS         comma-separated cohort whitelist (default PCAWG,TCGA,Mutographs)
# ESS_MUTO_KIDNEY_ONLY=1  keep Mutographs samples only where HD grouping == Kidney (= RCC)
# ESS_MUTO_EXCLUDE_HDGROUP  comma-separated Cancer_Type_HD_Grouping values to DROP from
#                     Mutographs only (e.g. "Pancreas"); other cohorts keep those tissues
# ESS_MUTO_EXCLUDE_DETAILS  comma-separated Cancer_Type_Details_Publication values to DROP
#                     from Mutographs only (e.g. "EGCA" to drop esophageal adenocarcinoma,
#                     keeping ESCC squamous)
# ESS_EXCLUDE_SAMPLES_FILE  path to a newline-delimited list of Sample_IDs to DROP entirely
#                     (any cohort); used to remove a named sample series (e.g. the Balkan
#                     ccRCC list) before any statistics
# ESS_ANALYSIS_DIR    output subfolder name each run2 script writes into (default "analysis")
COHORTS = os.environ.get("ESS_COHORTS", "PCAWG,TCGA,Mutographs").split(",")
MUTO_KIDNEY_ONLY = os.environ.get("ESS_MUTO_KIDNEY_ONLY") == "1"
MUTO_EXCLUDE_HDGROUP = [x for x in os.environ.get("ESS_MUTO_EXCLUDE_HDGROUP", "").split(",") if x]
MUTO_EXCLUDE_DETAILS = [x for x in os.environ.get("ESS_MUTO_EXCLUDE_DETAILS", "").split(",") if x]
EXCLUDE_SAMPLES_FILE = os.environ.get("ESS_EXCLUDE_SAMPLES_FILE", "")
EXCLUDE_SAMPLES = set()
if EXCLUDE_SAMPLES_FILE and os.path.exists(EXCLUDE_SAMPLES_FILE):
    with open(EXCLUDE_SAMPLES_FILE) as _fh:
        EXCLUDE_SAMPLES = {ln.strip() for ln in _fh if ln.strip()}
ANALYSIS_DIRNAME = os.environ.get("ESS_ANALYSIS_DIR", "analysis")
# subset is "active" whenever any of the cohort-subset knobs are set; when inactive the
# activity loaders behave exactly as before so the original full-cohort outputs are untouched.
SUBSET_ACTIVE = bool(os.environ.get("ESS_COHORTS") or os.environ.get("ESS_MUTO_KIDNEY_ONLY")
                     or MUTO_EXCLUDE_HDGROUP or MUTO_EXCLUDE_DETAILS or EXCLUDE_SAMPLES)
MB = 2800.0                 # whole-genome for ALL cohorts (see memory)
MIN_N = 20                  # min ever & min never per tissue for stratified stats
Q = 0.05
ACT_MIN = 0.05              # eSS "present"/carrier if relative activity >= 5%
col_sm, col_ns, col_g = "#d1495b", "#3a7ca5", "#bdbdbd"

_full, short_lab, _cos = load_labels(COMP)


def lab(s):
    return short_lab.get(s, "novel") if s.startswith("eSS") else s


# ------------------------------------------------------------------ eSS -> COSMIC (TMB tags)
cosmic = pd.read_csv(COMP)
cosmic["top1"] = cosmic["COSMIC_top3"].str.split(";").str[0].str.strip()
cosmic["top1_sbs"] = cosmic["top1"].str.replace(r"\s*\(.*\)", "", regex=True)
ess_to_cosmic = dict(zip(cosmic["eSS"], cosmic["top1"]))
ess_to_cosmic_sbs = dict(zip(cosmic["eSS"], cosmic["top1_sbs"]))

SBS_ETIO = {
    "SBS5": "Clock-like", "SBS40a": "Clock-like", "SBS40c": "Clock-like",
    # SBS100 = COSMIC "associated with tobacco smoking" (lung, smokers) -> Tobacco,
    # NOT the Unknown bucket its v3.6 batch-mates (SBS95/101/102/111) sit in.
    "SBS4": "Tobacco", "SBS29": "Tobacco", "SBS100": "Tobacco",
    "SBS7a": "UV", "SBS38": "UV",
    "SBS11": "Therapy", "SBS31": "Therapy", "SBS90": "Therapy",
    "SBS22a": "Aristolochic acid", "SBS22c": "Aristolochic acid",
    "SBS24": "Aflatoxin",   # COSMIC: "Aflatoxin exposure"
    "SBS88": "Colibactin", "SBS18": "ROS / oxidative",
    "SBS36": "Defective DNA repair", "SBS85": "AID / immune",
    "SBS95": "Unknown (COSMIC)",
    "SBS101": "Unknown (COSMIC)", "SBS102": "Unknown (COSMIC)",
    "SBS111": "Unknown (COSMIC)",
}
ETIO_COLOR = {
    "Clock-like": "#9E9E9E", "Tobacco": "#8B4513", "UV": "#FFB300",
    "Therapy": "#6A3D9A", "Aristolochic acid": "#E31A1C", "Aflatoxin": "#00838F",
    "Colibactin": "#1F78B4",
    "ROS / oxidative": "#33A02C", "Defective DNA repair": "#FB9A99",
    "AID / immune": "#FF7F00", "Unknown (COSMIC)": "#B2DF8A", "novel": "#FFFFFF",
}
ETIO_ORDER = ["Clock-like", "Tobacco", "UV", "Therapy", "Aristolochic acid",
              "Aflatoxin", "Colibactin", "ROS / oxidative", "Defective DNA repair",
              "AID / immune", "Unknown (COSMIC)", "novel"]
# Curated (COSMIC tag, etiology) overrides. RETIRED 2026-07-23: eSS37 used to be forced to
# ("SBS92-like", "Tobacco"). We no longer make that smoking association. SBS92 is a tobacco
# signature and was never in eSS37's top-3 anyway (SBS22c 0.86; SBS12 0.81; SBS26 0.79), so
# the tag was a hand call, not a cosine result. eSS37 now follows the standard 0.90 rule ->
# "novel"; its real training identity (BCA, bromochloroacetic acid) is carried by the
# contributor / group / umbrella dot-matrix tags, where exposure identity belongs.
ESS_OVERRIDE = {}
# ---- etiology overrides from the AUTHORITATIVE training identity (contributor summary) ----
# Applied where the top-1 COSMIC cosine picks a signature whose COSMIC aetiology is wrong for
# this eSS. Only the ETIOLOGY (colour/legend) is overridden; ess_cosmic_tag still reports the
# honest best cosine match, so the COSMIC dot-matrix and volcanoes are unchanged.
#
# eSS10 and eSS11 are BOTH AFB1-trained, but their top-1 matches differ by a hair and used to
# split them across two unrelated buckets:
#   eSS10  SBS29 (0.90) / SBS24 (0.88) / SBS95 (0.84)  -> SBS29 = COSMIC "tobacco chewing"
#   eSS11  SBS95 (0.94) / SBS29 (0.91) / SBS100 (0.90) -> SBS95 = COSMIC "possible artefact"
# COSMIC's aflatoxin signature is SBS24, which is eSS10's #2 at 0.88 and is absent from
# eSS11's top-3. Since the contributor file (not the cosine) is the identity source for eSS,
# both are labelled Aflatoxin. NOTE eSS11's only carriers are 7 low-TMB lymphoid genomes -
# see the eSS11 caveat in the analysis METHODS before interpreting it as real exposure.
ESS_ETIO_OVERRIDE = {"eSS10": "Aflatoxin", "eSS11": "Aflatoxin"}
# eSS50 / eSS51 extend the 49-eSS catalogue to cover the two exposures the 49 miss
# (SBS7d UV and SBS17b/5-FU; see memory "eSS coverage gaps"). They have no row in the
# contributor summary, so their identity is fixed here, version-aware so each plot stays
# in its own convention: COSMIC plot -> SBS name, contributor/group plots -> exposure.
ESS_EXTRA = {
    "eSS50": {"cosmic": "SBS7d",  "contrib": "UVA",  "etio": "UV"},
    "eSS51": {"cosmic": "SBS17b", "contrib": "5-FU", "etio": "Therapy"},
}


COSMIC_COS_MIN = 0.90   # min cosine to assign an eSS a COSMIC identity/etiology (else novel)


def _cosval(c):
    m = re.search(r"\(([\d.]+)\)", ess_to_cosmic.get(c, ""))
    return float(m.group(1)) if m else 0.0


def ess_cosmic_tag(c):
    if c in ESS_EXTRA:
        return ESS_EXTRA[c]["cosmic"]
    if c in ESS_OVERRIDE:
        return ESS_OVERRIDE[c][0]
    return ess_to_cosmic_sbs.get(c, "") if _cosval(c) >= COSMIC_COS_MIN else "novel"


def ess_etiology(c):
    if c in ESS_ETIO_OVERRIDE:      # curated identity beats the top-1 cosine
        return ESS_ETIO_OVERRIDE[c]
    if c in ESS_EXTRA:
        return ESS_EXTRA[c]["etio"]
    if c in ESS_OVERRIDE:
        return ESS_OVERRIDE[c][1]
    if _cosval(c) >= COSMIC_COS_MIN:
        return SBS_ETIO.get(ess_to_cosmic_sbs.get(c, ""), "Unknown (COSMIC)")
    return "novel"


# ---------------------------------------------------------- eSS -> TOP TRAINING CONTRIBUTOR
# Alternative row tag for the dotmatrix: instead of the COSMIC-cosine SBS label, use the
# dominant mutagen the eSS was trained on (e.g. eSS33 -> "AAI" because aristolochic acid I
# is the majority of its training contributors).  Source = the authoritative Cardiff
# per-eSS contributor summary (see project memory "eSS identity source").
CONTRIB_SUMMARY = os.environ.get(
    "ESS_CONTRIB_SUMMARY",
    os.path.join(BASE, "for_decomposition_cluster_summary_for_pie_charts.txt"))

# compound-token -> short acronym.  Ordered longest / most-specific FIRST so combined and
# co-factor tokens (HBV_AFB1, HPV_4NQO, Arsenic_SS, CX-5461_UVA, N-OH-PhIP, *_fat_diet)
# resolve to their intended label before shorter substrings match.
_CONTRIB_PATTERNS = [
    ("Normal_fat_diet", "NFD"), ("High_fat_diet", "HFD"), ("N-OH-PhIP", "N-OH-PhIP"),
    ("Bromochloroacetic_acid", "BCA"), ("Pyrrolizidine_Alkaloid", "PA"),
    ("Diethylnitrosamine", "DEN"), ("3-Nitrobenzanthrone", "3-NBA"),
    ("6-Nitrochrysene", "6-NC"), ("5-Methylchrysene", "5-MC"),
    ("1,8-Dinitropyrene", "1,8-DNP"), ("Potassium_bromate", "KBrO3"),
    ("4-aminobiphenyl", "4-ABP"), ("Cobalt_metal", "Cobalt"), ("Temozolomide", "TMZ"),
    ("Duocarmycin", "Duocarmycin"), ("Colibactin", "Colibactin"),
    ("Glycidamide", "Glycidamide"), ("Acrylamide", "Acrylamide"), ("Cisplatin", "Cisplatin"),
    ("1.2.3_TCP", "1,2,3-TCP"), ("CX-5461_UVA", "CX-5461+UVA"), ("CX-5461", "CX-5461"),
    ("TMP-UVA200J", "TMP+UVA"), ("TMP-UVA", "TMP+UVA"), ("DBA-DE", "DBA-DE"),
    ("DBPDE", "DBPDE"), ("DBP", "DBP"), ("DBA", "DBA"), ("BPDE", "BPDE"), ("BaP", "BaP"),
    ("AFB1", "AFB1"), ("AAI", "AAI"), ("DMBA", "DMBA"), ("PhIP", "PhIP"), ("4NQO", "4NQO"),
    ("MNNG", "MNNG"), ("MNU", "MNU"), ("NNK", "NNK"), ("ENU", "ENU"), ("EMS", "EMS"),
    ("MMS", "MMS"), ("DMS", "DMS"), ("GammaRay", "γRay"), ("Fe-ions", "Fe-ion"),
    ("Arsenic_SS", "As+SS"), ("UVA", "UVA"), ("UVB", "UVB"), ("UVC", "UVC"), ("SS", "SS"),
    ("Arsenic", "Arsenic"), ("VC", "VC"), ("ATC", "ATC"),
]


def _contrib_compound(token):
    for pat, acr in _CONTRIB_PATTERNS:
        if re.search(re.escape(pat), token, flags=re.IGNORECASE):
            return acr
    return token  # unmatched -> surfaces verbatim so it is noticed


# compound acronym -> chemical class, used to relabel MIXED eSS (no single dominant
# compound) with the shared class of their contributors instead of an arbitrary winner.
CONTRIB_CLASS = {
    "ATC": "Nucleoside analog", "Cobalt": "Metal",
    "γRay": "Ionizing radiation", "Fe-ion": "Ionizing radiation",
    "NFD": "Fat diet", "HFD": "Fat diet",
    "UVA": "UV/solar", "UVB": "UV/solar", "UVC": "UV/solar", "SS": "UV/solar", "As+SS": "UV/solar",
    "4-ABP": "Aromatic amine", "PhIP": "Aromatic amine", "N-OH-PhIP": "Aromatic amine",
    "Acrylamide": "Acrylamide", "Glycidamide": "Acrylamide",
    "AFB1": "Aflatoxin", "4NQO": "UV-mimetic", "PA": "Pyrrolizidine alkaloid",
    "Cisplatin": "Platinum", "KBrO3": "Oxidative",
    "3-NBA": "Nitro-PAH", "6-NC": "Nitro-PAH", "1,8-DNP": "Nitro-PAH",
    "BaP": "PAH", "BPDE": "PAH", "DBA": "PAH", "DBA-DE": "PAH", "DBP": "PAH",
    "DBPDE": "PAH", "5-MC": "PAH", "DMBA": "PAH",
    "1,2,3-TCP": "Haloalkane", "VC": "Haloalkane", "BCA": "Haloacetic acid",
    "MNU": "Alkylating", "MNNG": "Alkylating", "TMZ": "Alkylating", "EMS": "Alkylating",
    "MMS": "Alkylating", "DMS": "Alkylating", "ENU": "Alkylating", "DEN": "Alkylating",
    "NNK": "Alkylating",
    "TMP+UVA": "Photochemo (PUVA)", "CX-5461": "CX-5461", "CX-5461+UVA": "CX-5461",
    "Duocarmycin": "Duocarmycin", "Colibactin": "Colibactin",
}
MIX_FRAC = 0.60  # dominant-compound fraction below this = "mixture" -> use chemical class

# umbrella super-group: broader / short class that can represent SEVERAL chemical classes
# at once, so a mixture spanning >1 class collapses to one covering label (nitro-PAH + PAH
# -> "PAH"), or joins compact acronyms (aromatic amine + acrylamide -> "AA/Dietary").
# Classes not listed map to themselves. Coined acronyms are spelled out in ACRONYM_NOTE.
CLASS_SUPERGROUP = {
    "Aromatic amine": "AA",       # AA = aromatic amines
    "Acrylamide": "Dietary",      # Dietary = dietary / cooked-food mutagen (acrylamide)
    "PAH": "PAH", "Nitro-PAH": "PAH",
}
# short legend notes explaining the coined acronyms, per plot (only list what appears):
# γRay is used everywhere the compound acronym shows; AA/Dietary only in the umbrella tags.
GAMMA_NOTE = "γRay = γ-radiation"
ACRONYM_NOTE = ("AA = aromatic amines\n"
                "Dietary = dietary / cooked-food\n"
                "γRay = γ-radiation")


def _supergroup(cls):
    return CLASS_SUPERGROUP.get(cls, cls)


def _load_contrib_labels(path):
    """eSS -> {'compound': (acr, frac), 'group': (class, frac), 'classes': set}
    from the contributor counts aggregated by compound and by chemical class."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            ess = line.split(",", 1)[0]
            rest = line[len(ess) + 1:]
            comp_agg, cls_agg, total = {}, {}, None
            for key, cnt in re.findall(r'([\w.\-\[\],+]+):(\d+)', rest):
                key = key.strip(",")
                if key == "total":
                    total = int(cnt)
                    continue
                comp = _contrib_compound(key)
                comp_agg[comp] = comp_agg.get(comp, 0) + int(cnt)
                cls = CONTRIB_CLASS.get(comp, comp)
                cls_agg[cls] = cls_agg.get(cls, 0) + int(cnt)
            if not comp_agg:
                continue
            tot = total if total else sum(comp_agg.values())
            tc, nc = max(comp_agg.items(), key=lambda kv: kv[1])
            tg, ng = max(cls_agg.items(), key=lambda kv: kv[1])
            out[ess] = {"compound": (tc, nc / tot), "group": (tg, ng / tot),
                        "classes": set(cls_agg)}
    return out


ESS_CONTRIB = _load_contrib_labels(CONTRIB_SUMMARY)
# back-compat: eSS -> (top_compound, fraction)
ESS_TOP_CONTRIB = {e: v["compound"] for e, v in ESS_CONTRIB.items()}


def ess_contrib_tag(c):
    """Dotmatrix row tag = dominant training-mutagen acronym (falls back to novel)."""
    if c in ESS_EXTRA:
        return ESS_EXTRA[c]["contrib"]
    hit = ESS_CONTRIB.get(c)
    return hit["compound"][0] if hit else "novel"


def ess_group_tag(c):
    """Row tag = compound acronym for a clear-majority eSS, else the DOMINANT chemical
    class of its contributors (for mixtures where no compound reaches MIX_FRAC)."""
    if c in ESS_EXTRA:
        return ESS_EXTRA[c]["contrib"]
    hit = ESS_CONTRIB.get(c)
    if not hit:
        return "novel"
    acr, frac = hit["compound"]
    return acr if frac >= MIX_FRAC else hit["group"][0]


def ess_umbrella_tag(c):
    """Row tag = compound acronym for a clear-majority eSS, else a chemical group broad
    enough to cover ALL of its contributors (not just the dominant class): the shared
    super-group if the contributor classes collapse to one, else their joined labels."""
    if c in ESS_EXTRA:
        return ESS_EXTRA[c]["contrib"]
    hit = ESS_CONTRIB.get(c)
    if not hit:
        return "novel"
    acr, frac = hit["compound"]
    if frac >= MIX_FRAC:
        return acr
    supers = sorted({_supergroup(cl) for cl in hit["classes"]})
    return supers[0] if len(supers) == 1 else "/".join(supers)


def _fmt(x):
    return f"{int(round(x)):,}" if x >= 1 else ("%g" % x)


# ------------------------------------------------------------------ metadata (loaded once)
def classify(v):
    if pd.isna(v):
        return np.nan
    s = str(v).strip().lower()
    if s in {"lifelong non-smoker", "lifelong non-smoker (<100 cigarettes smoked in lifetime)",
             "never", "no", "non-smoker", "0"}:
        return 0
    if s in {"current smoker", "current smoker (includes daily smokers non-daily/occasional smokers)",
             "smoker", "yes", "1", "ever smoker", "ex-smoker"} or "reformed smoker" in s:
        return 1
    return np.nan


meta = pd.read_csv(META_FILE, sep="\t", dtype=str, low_memory=False)
meta = meta[meta["Sample_ID"] != "Sample_ID"]
meta = meta[meta["Cohort"].isin(COHORTS)]
if MUTO_KIDNEY_ONLY:
    # "Mutograph RCC only": every Mutographs Kidney sample is RCC/ccRCC in the
    # publication-detail column, so HD grouping == Kidney is an exact RCC filter.
    meta = meta[~((meta["Cohort"] == "Mutographs") &
                  (meta["Cancer_Type_HD_Grouping"] != "Kidney"))]
if MUTO_EXCLUDE_HDGROUP:
    # Drop named HD groupings from Mutographs only (other cohorts keep those tissues).
    meta = meta[~((meta["Cohort"] == "Mutographs") &
                  (meta["Cancer_Type_HD_Grouping"].isin(MUTO_EXCLUDE_HDGROUP)))]
if MUTO_EXCLUDE_DETAILS:
    # Drop named publication-detail subtypes from Mutographs only (e.g. EGCA = esophageal
    # adenocarcinoma, keeping ESCC squamous).
    meta = meta[~((meta["Cohort"] == "Mutographs") &
                  (meta["Cancer_Type_Details_Publication"].isin(MUTO_EXCLUDE_DETAILS)))]
if EXCLUDE_SAMPLES:
    # Drop a named sample series by Sample_ID (any cohort).
    meta = meta[~meta["Sample_ID"].isin(EXCLUDE_SAMPLES)]
meta["smoker"] = meta["Smoking_Status"].apply(classify)
meta["age"] = pd.to_numeric(meta["Age"], errors="coerce")
meta["female"] = meta["Sex"].map({"Female": 1, "Male": 0})
meta["ctype"] = meta["Cancer_Type_HD_Grouping"]
META = meta[["Sample_ID", "smoker", "age", "female", "ctype"]].rename(
    columns={"Sample_ID": "Samples"})
# allowed-sample whitelist used to subset the activity loaders (see load_pooled/load_single)
ALLOWED_SAMPLES = set(META["Samples"])


# ------------------------------------------------------------------ activities loaders
def _read_act(path):
    a = pd.read_csv(path, sep="\t")
    a = a.rename(columns={a.columns[0]: "Samples"})
    return a


def _subset_rows(a):
    """When a cohort subset is active, drop activity rows for samples outside the
    allowed set. No-op otherwise, so canonical full-cohort runs are byte-identical."""
    if SUBSET_ACTIVE:
        return a[a["Samples"].isin(ALLOWED_SAMPLES)].reset_index(drop=True)
    return a


def load_single(path):
    return _subset_rows(_read_act(path))


def load_pooled(parent, sub="Assignment_Solution/Activities/Assignment_Solution_Activities.txt"):
    """Pool per-cancer-type folders into one matrix (union of columns, missing=0)."""
    frames = []
    for d in sorted(glob.glob(os.path.join(parent, "*"))):
        if not os.path.isdir(d):
            continue
        p = os.path.join(d, sub)
        if os.path.exists(p):
            frames.append(_read_act(p))
    if not frames:
        raise FileNotFoundError(f"no per-cancer Activities under {parent}")
    all_cols = []
    for f in frames:
        for c in f.columns:
            if c not in all_cols:
                all_cols.append(c)
    frames = [f.reindex(columns=all_cols).fillna(0.0) if list(f.columns) != all_cols else f
              for f in frames]
    out = pd.concat(frames, ignore_index=True)
    num = [c for c in out.columns if c != "Samples"]
    out[num] = out[num].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return _subset_rows(out.drop_duplicates(subset="Samples"))


# ================================================================== TMB figures
def beeswarm_tmb(groups, outfile, top_rot=40, bottom_line1=None, bottom_line2=None,
                 bottom_colours=None, bottom_font=12, top_in=2.4, col_in=0.62,
                 bottom_rot=45, head_count="n", head1="", head2="", line_gap=0.15):
    groups = [g for g in groups if len(g[1]) > 0]
    if not groups:
        return
    groups.sort(key=lambda g: np.median(g[1]))
    names = [g[0] for g in groups]
    N = len(names)
    allv = np.concatenate([g[1] for g in groups])
    ymin, ymax = np.floor(allv.min()), np.ceil(allv.max())
    if ymax <= ymin:
        ymax = ymin + 1
    ticks = np.arange(ymin, ymax + 1)
    yr = ymax - ymin
    col_w = 1.0
    pad = 0.12 * col_w
    COL_IN, LEFT_IN, RIGHT_IN, TOP_IN, PANEL_H = col_in, 2.1, 0.4, top_in, 4.6
    sinr = np.sin(np.deg2rad(bottom_rot))
    lbl_v = (7 * bottom_font * 0.6 / 72.0) * sinr
    th_in = bottom_font / 72.0
    cnt_off = 0.28
    l1_off = cnt_off + th_in + line_gap
    l2_off = l1_off + lbl_v + line_gap
    BOT_IN = l2_off + (bottom_font + 1) / 72.0 + 0.18
    panel_w_in = N * COL_IN
    fig_w = LEFT_IN + panel_w_in + RIGHT_IN
    fig_h = TOP_IN + PANEL_H + BOT_IN
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([LEFT_IN / fig_w, BOT_IN / fig_h, panel_w_in / fig_w, PANEL_H / fig_h])
    u_per_in = yr / PANEL_H
    bfont, tfont = bottom_font, 16
    for i in range(N):
        face = "#E8E8E8" if i % 2 == 0 else "#FFFFFF"
        ax.add_patch(Rectangle((i * col_w, ymin), col_w, yr, facecolor=face, edgecolor="none", zorder=0))
    for i, (name, vals) in enumerate(groups):
        vals = np.sort(vals)
        n = len(vals)
        xs = np.linspace(i * col_w + pad, (i + 1) * col_w - pad, n) if n > 1 else np.array([i * col_w + col_w / 2])
        ax.scatter(xs, vals, s=8, color="black", zorder=2, edgecolors="none")
        med = np.median(vals)
        ax.plot([i * col_w + pad, (i + 1) * col_w - pad], [med, med], color="red", lw=3, zorder=3, solid_capstyle="butt")
        cx = i * col_w + col_w / 2
        ax.text(cx, ymax + 0.03 * yr, name, ha="left", va="bottom", rotation=top_rot, fontsize=tfont, rotation_mode="anchor")
        ax.text(cx, ymin - cnt_off * u_per_in, str(n), ha="center", va="top", fontsize=bfont - 1)
        if bottom_line1 is not None:
            col = bottom_colours.get(name, "black") if bottom_colours else "black"
            ax.text(cx, ymin - l1_off * u_per_in, bottom_line1.get(name, ""), ha="right", va="top",
                    rotation=bottom_rot, rotation_mode="anchor", fontsize=bfont, color=col, fontweight="bold")
            ax.text(cx, ymin - l2_off * u_per_in, bottom_line2.get(name, ""), ha="center", va="top",
                    fontsize=bfont + 1, color=col)
    hx = -0.4
    th = bfont / 72.0
    ax.text(hx, ymin - (cnt_off + th / 2) * u_per_in, head_count, ha="right", va="center", fontsize=bfont, clip_on=False)
    ax.text(hx, ymin - (l1_off + lbl_v / 2) * u_per_in, head1, ha="right", va="center", fontsize=bfont, fontweight="bold", clip_on=False)
    ax.text(hx, ymin - (l2_off + th / 2) * u_per_in, head2, ha="right", va="center", fontsize=bfont, clip_on=False)
    for t in ticks:
        ax.axhline(t, color="black", lw=1.0 if t in (ymin, ymax) else 0.6,
                   ls="-" if t in (ymin, ymax) else (0, (2, 3)), zorder=1)
    ax.set_xlim(0, N * col_w)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks([])
    ax.set_yticks(ticks)
    ax.set_yticklabels([_fmt(10 ** t) for t in ticks], fontsize=14)
    ax.tick_params(axis="y", length=0, pad=6)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.set_ylabel("Number of eSS mutations per megabase", fontsize=15, fontweight="bold", labelpad=18)
    fig.savefig(outfile, transparent=True)
    plt.close(fig)


def _sig_labels(groups_sig):
    l1, l2 = {}, {}
    for c, _ in groups_sig:
        cs = _cosval(c)
        if c in ESS_OVERRIDE:
            l1[c], l2[c] = ESS_OVERRIDE[c][0], "-"
        elif cs >= 0.85:
            l1[c], l2[c] = ess_to_cosmic_sbs.get(c, ""), f"{cs:.2f}"
        else:
            l1[c], l2[c] = "novel", "-"
    return l1, l2


def bubble_tmb(df, present, ess_cols, outfile, top_rot=40, col_in=0.95, top_in=2.8,
               min_prev=0.05, label_prev=0.20, size_scale=11.0):
    ess_per_mb = df[ess_cols] / MB
    tot_mb = df[ess_cols].sum(axis=1) / MB
    order = (tot_mb[tot_mb > 0].groupby(df["CancerType"]).median().sort_values().index.tolist())
    N = len(order)
    xs, ys, sizes, cols, labels = [], [], [], [], []
    used_etio = set()
    for i, ct in enumerate(order):
        m = df["CancerType"].values == ct
        for c in ess_cols:
            prev = present[c].values[m].mean()
            if prev < min_prev:
                continue
            carriers = m & present[c].values
            med = np.median(ess_per_mb[c].values[carriers])
            if med <= 0:
                continue
            etio = ess_etiology(c)
            used_etio.add(etio)
            xs.append(i + 0.5)
            ys.append(np.log10(med))
            sizes.append(prev * 100 * size_scale)
            cols.append(ETIO_COLOR[etio])
            labels.append((i + 0.5, np.log10(med), c.replace("eSS", ""), prev))
    if not ys:
        return set()
    ymin, ymax = np.floor(min(ys)), np.ceil(max(ys))
    ticks = np.arange(ymin, ymax + 1)
    LEFT_IN, RIGHT_IN, PANEL_H, BOT_IN = 2.1, 1.8, 5.4, 0.5
    panel_w_in = N * col_in
    fig_w = LEFT_IN + panel_w_in + RIGHT_IN
    fig_h = top_in + PANEL_H + BOT_IN
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([LEFT_IN / fig_w, BOT_IN / fig_h, panel_w_in / fig_w, PANEL_H / fig_h])
    yr = ymax - ymin
    for k in range(N):
        face = "#E8E8E8" if k % 2 == 0 else "#FFFFFF"
        ax.add_patch(Rectangle((k, ymin), 1, yr, facecolor=face, edgecolor="none", zorder=0))
    for t in ticks:
        ax.axhline(t, color="black", lw=1.0 if t in (ymin, ymax) else 0.6,
                   ls="-" if t in (ymin, ymax) else (0, (2, 3)), zorder=1)
    ax.scatter(xs, ys, s=sizes, c=cols, alpha=0.7, edgecolors="black", linewidths=0.5, zorder=3)
    for x, y, l, prev in labels:
        if prev >= label_prev:
            ax.text(x, y, l, ha="center", va="center", fontsize=8, fontweight="bold", zorder=4)
    for i, ct in enumerate(order):
        ax.text(i + 0.5, ymax + 0.03 * yr, ct, ha="left", va="bottom", rotation=top_rot, fontsize=16, rotation_mode="anchor")
    ax.set_xlim(0, N)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks([])
    ax.set_yticks(ticks)
    ax.set_yticklabels([_fmt(10 ** t) for t in ticks], fontsize=14)
    ax.tick_params(axis="y", length=0, pad=6)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.set_ylabel("Median eSS mutations per megabase (carriers)", fontsize=15, fontweight="bold", labelpad=18)
    fracs = [0.10, 0.25, 0.50, 1.00]
    ylegs = np.linspace(ymax - 0.4, ymax - 0.4 - 3 * 0.5, 4)
    for frac, yleg in zip(fracs, ylegs):
        ax.scatter([N + 0.6], [yleg], s=frac * 100 * size_scale, c="#999999", alpha=0.7,
                   edgecolors="black", linewidths=0.5, clip_on=False, zorder=3)
        ax.text(N + 1.1, yleg, f"{int(frac*100)}%", ha="left", va="center", fontsize=12, clip_on=False)
    ax.text(N + 0.6, ymax - 0.1, "prevalence", ha="center", va="bottom", fontsize=12, fontweight="bold", clip_on=False)
    fig.savefig(outfile, transparent=True)
    plt.close(fig)
    return used_etio


def etiology_legend(used_etio, outfile):
    etios = [e for e in ETIO_ORDER if e in used_etio]
    n = len(etios)
    if n == 0:
        return
    fig = plt.figure(figsize=(3.0, 0.3 * n + 0.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n + 0.8)
    ax.axis("off")
    ax.text(0.05, n + 0.35, "COSMIC etiology", ha="left", va="center", fontsize=13, fontweight="bold")
    for k, e in enumerate(etios):
        y = n - 1 - k + 0.4
        col = ETIO_COLOR[e]
        if col == "#FFFFFF":          # 'novel' is drawn black in the dotmatrix -> match here
            col = "#000000"
        ax.scatter([0.1], [y], s=140, marker="o", facecolor=col, edgecolors="black", linewidths=0.6)
        ax.text(0.2, y, e, ha="left", va="center", fontsize=12)
    fig.savefig(outfile, transparent=True)
    plt.close(fig)


def dotmatrix_tmb(df, present, ess_cols, outfile, min_prev=0.0, cell_in=0.28, size_scale=3.4,
                  tag_fn=None, note=None, color_labels=True, legend_horizontal=False):
    # color_labels: colour the row labels by COSMIC etiology (default); False -> plain black.
    tag_fn = tag_fn or ess_cosmic_tag
    ess_per_mb = df[ess_cols] / MB
    order_ct = sorted(pd.unique(df["CancerType"]))
    Nc = len(order_ct)
    masks = {ct: (df["CancerType"].values == ct) for ct in order_ct}
    keep = [c for c in ess_cols if max(present[c].values[masks[ct]].sum() for ct in order_ct) > 0]
    rows = sorted(keep, key=lambda c: int(c.replace("eSS", "")))
    Nr = len(rows)
    xs, ys, sizes, vals = [], [], [], []
    for xi, ct in enumerate(order_ct):
        m = masks[ct]
        for yi, c in enumerate(rows):
            carriers = m & present[c].values
            if carriers.sum() == 0:
                continue
            med = np.median(ess_per_mb[c].values[carriers])
            if med <= 0:
                continue
            prev = present[c].values[m].mean()
            if prev < min_prev:
                continue
            xs.append(xi + 0.5)
            ys.append(Nr - 1 - yi + 0.5)
            sizes.append(prev * 100 * size_scale)
            vals.append(np.log10(med))
    if not vals:
        return
    vmin, vmax = np.floor(min(vals)), np.ceil(max(vals))
    norm = plt.Normalize(vmin, vmax)
    cmap = plt.get_cmap("viridis")
    LEFT_IN, RIGHT_IN, TOP_IN, BOT_IN = 2.2, 2.4, 2.8, 0.4
    panel_w, panel_h = Nc * cell_in, Nr * cell_in
    fig_w = LEFT_IN + panel_w + RIGHT_IN
    fig_h = TOP_IN + panel_h + BOT_IN
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([LEFT_IN / fig_w, BOT_IN / fig_h, panel_w / fig_w, panel_h / fig_h])
    for xi in range(Nc):
        ax.axvline(xi, color="#EEEEEE", lw=0.6, zorder=0)
    for yi in range(Nr):
        ax.axhline(yi, color="#EEEEEE", lw=0.6, zorder=0)
    sc = ax.scatter(xs, ys, s=sizes, c=vals, cmap=cmap, norm=norm, edgecolors="black", linewidths=0.4, zorder=3)
    ax.set_xlim(0, Nc)
    ax.set_ylim(0, Nr)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    for xi, ct in enumerate(order_ct):
        ax.text(xi + 0.5, Nr + 0.2, ct, ha="left", va="bottom", rotation=40, rotation_mode="anchor", fontsize=11)
    for yi, c in enumerate(rows):
        y = Nr - 1 - yi + 0.5
        tag = tag_fn(c)
        ecol = ETIO_COLOR[ess_etiology(c)] if color_labels else "#000000"
        if ecol == "#FFFFFF":
            ecol = "#000000"
        ax.text(-0.25, y, f"{c} ({tag})", ha="right", va="center", fontsize=9, fontweight="bold", color=ecol)
    lx = (LEFT_IN + panel_w + 0.35) / fig_w
    panel_top = BOT_IN + panel_h
    LG_H = 1.5
    if legend_horizontal:
        # horizontal colorbar (ticks + label below the bar), prevalence key underneath
        CB_BAR_H, CB_W, GAP = 0.16, RIGHT_IN - 0.7, 0.80
        cb_top = panel_top - 0.15
        cbar_ax = fig.add_axes([lx, (cb_top - CB_BAR_H) / fig_h, CB_W / fig_w, CB_BAR_H / fig_h])
        cb = fig.colorbar(sc, cax=cbar_ax, orientation="horizontal")
        cb.ax.tick_params(labelsize=7.5)
        lg_top = (cb_top - CB_BAR_H) - GAP
    else:
        CB_H, GAP = 1.3, 0.55
        cb_y = panel_top - CB_H
        cbar_ax = fig.add_axes([lx, cb_y / fig_h, 0.14 / fig_w, CB_H / fig_h])
        cb = fig.colorbar(sc, cax=cbar_ax)
        cb.ax.tick_params(labelsize=8)
        lg_top = cb_y - GAP
    cb.set_ticks(np.arange(vmin, vmax + 1))
    cb.set_ticklabels([_fmt(10 ** t) for t in np.arange(vmin, vmax + 1)])
    cb.outline.set_linewidth(0.4)
    cb.set_label("median mut/Mb", fontsize=9, fontweight="bold")
    lax = fig.add_axes([lx, (lg_top - LG_H) / fig_h, (RIGHT_IN - 0.5) / fig_w, LG_H / fig_h])
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    lax.text(0.0, 0.98, "prevalence", ha="left", va="top", fontsize=9, fontweight="bold")
    if legend_horizontal:
        # same bubbles/size, just laid out left-to-right with the % under each
        for frac, px in zip([0.10, 0.25, 0.50, 1.00], [0.13, 0.38, 0.63, 0.88]):
            lax.scatter([px], [0.60], s=frac * 100 * size_scale, c="#777777",
                        edgecolors="black", linewidths=0.4)
            lax.text(px, 0.14, f"{int(frac*100)}%", ha="center", va="center", fontsize=8)
    else:
        for frac, py in zip([0.10, 0.25, 0.50, 1.00], [0.74, 0.53, 0.32, 0.11]):
            lax.scatter([0.16], [py], s=frac * 100 * size_scale, c="#777777", edgecolors="black", linewidths=0.4)
            lax.text(0.45, py, f"{int(frac*100)}%", ha="left", va="center", fontsize=8)
    if note:
        note_h = 1.4
        note_top = (lg_top - LG_H) - 0.4          # inches from fig bottom, just under legend
        nax = fig.add_axes([lx, (note_top - note_h) / fig_h,
                            (RIGHT_IN - 0.35) / fig_w, note_h / fig_h])
        nax.set_xlim(0, 1)
        nax.set_ylim(0, 1)
        nax.axis("off")
        nax.text(0.0, 1.0, "acronyms", ha="left", va="top", fontsize=9, fontweight="bold")
        nax.text(0.0, 0.80, note, ha="left", va="top", fontsize=7.5, linespacing=1.6)
    fig.savefig(outfile, transparent=True)
    plt.close(fig)


def tmb_block(df, ess_cols, outdir):
    """All TMB figures + underlying CSVs into outdir. df needs CancerType.

    Carrier/prevalence uses relative activity vs the sample's TOTAL mutation
    burden (all signature columns: eSS + endogenous COSMIC + background), NOT the
    eSS-only subtotal, so a 5% carrier means >=5% of the whole tumour's mutations.
    Burden (mut/Mb, dot colour) stays absolute and is unaffected.
    """
    counts = df[ess_cols]
    all_sigs = [c for c in df.columns if c not in ("Samples", "CancerType")]
    tot = df[all_sigs].sum(axis=1).replace(0, np.nan)   # TOTAL TMB denominator
    rel = counts.div(tot, axis=0)
    present = rel >= ACT_MIN
    per_mb = counts.div(MB)
    tot_mb = counts.sum(axis=1) / MB

    # per-sample CSV
    ps = pd.DataFrame({"Samples": df["Samples"].values, "CancerType": df["CancerType"].values,
                       "eSS_total": counts.sum(axis=1).values, "eSS_per_mb": tot_mb.values})
    ps.to_csv(os.path.join(outdir, "TMB_eSS_per_sample.csv"), index=False)
    # per-cancer summary CSV
    rows = []
    for ct, g in df.groupby("CancerType"):
        v = (g[ess_cols].sum(axis=1) / MB)
        v = v[v > 0]
        pf = present[df["CancerType"].values == ct].mean(axis=0)
        rows.append({"CancerType": ct, "n": len(g), "n_carriers": int((v > 0).sum()),
                     "median_eSS_per_mb": float(np.median(v)) if len(v) else 0.0,
                     "most_prevalent_eSS": pf.idxmax(), "prevalence": float(pf.max())})
    pd.DataFrame(rows).sort_values("median_eSS_per_mb", ascending=False).to_csv(
        os.path.join(outdir, "TMB_eSS_per_cancer_summary.csv"), index=False)

    # Plot 1: beeswarm per cancer type
    keep = (tot_mb > 0).values
    y1 = np.log10(tot_mb.values[keep])
    ct1 = df["CancerType"].values[keep]
    groups_ct = [(nm, y1[ct1 == nm]) for nm in pd.unique(ct1)]
    groups_ct = [(nm, v) for nm, v in groups_ct if len(v) > 0]
    prevalent, prevalent_pct = {}, {}
    for nm in df["CancerType"].unique():
        pf = present[df["CancerType"].values == nm].mean(axis=0)
        sig = pf.idxmax()
        prevalent[nm] = sig
        prevalent_pct[nm] = f"{pf[sig] * 100:.0f}%"
    beeswarm_tmb(groups_ct, os.path.join(outdir, "TMB_eSS_pan-cancer.pdf"), top_rot=40,
                 bottom_line1=prevalent, bottom_line2=prevalent_pct, bottom_font=12,
                 top_in=2.8, col_in=0.7, head_count="n", head1="eSS", head2="prevalence")
    # Plot 2: per eSS signature
    groups_sig = []
    for c in ess_cols:
        mask = (counts[c].values > 0) & present[c].values
        if mask.sum() > 0:
            groups_sig.append((c, np.log10(per_mb[c].values[mask])))
    if groups_sig:
        l1, l2 = _sig_labels(groups_sig)
        beeswarm_tmb(groups_sig, os.path.join(outdir, "TMB_eSS_per-signature.pdf"), top_rot=90,
                     bottom_line1=l1, bottom_line2=l2, bottom_font=15, top_in=1.8, col_in=0.7,
                     head_count="n", head1="COSMIC", head2="cosine sim", line_gap=0.32)
    # Plot 3 + legend
    used = bubble_tmb(df, present, ess_cols, os.path.join(outdir, "TMB_eSS_pancancer_bubble.pdf"))
    etiology_legend(used, os.path.join(outdir, "TMB_eSS_etiology_legend.pdf"))
    # Plot 4 (COSMIC-cosine row tags) + alternative version tagged by dominant training mutagen
    dotmatrix_tmb(df, present, ess_cols, os.path.join(outdir, "TMB_eSS_dotmatrix.pdf"))
    dotmatrix_tmb(df, present, ess_cols, os.path.join(outdir, "TMB_eSS_dotmatrix_contributor.pdf"),
                  tag_fn=ess_contrib_tag)
    dotmatrix_tmb(df, present, ess_cols, os.path.join(outdir, "TMB_eSS_dotmatrix_group.pdf"),
                  tag_fn=ess_group_tag)
    dotmatrix_tmb(df, present, ess_cols, os.path.join(outdir, "TMB_eSS_dotmatrix_umbrella.pdf"),
                  tag_fn=ess_umbrella_tag, legend_horizontal=True)


# ================================================================== smoking volcanoes
def _ess_cosmic_short(s):
    """COSMIC identity for an eSS at the 0.90-cosine rule (short_lab), but honouring the
    curated overrides for the beyond-49 / special signatures (eSS50=SBS7d, eSS51=SBS17b).
    eSS with cosine < 0.90 stay 'novel' -- including eSS37, whose retired 'SBS92-like'
    tobacco call is no longer made."""
    if s in ESS_EXTRA:
        return ESS_EXTRA[s]["cosmic"]
    if s in ESS_OVERRIDE:
        return ESS_OVERRIDE[s][0]
    return lab(s)


def _volcano_label(s):
    """Point label. eSS keeps its COSMIC identity ("eSS18 | SBS100"); a COSMIC signature
    is shown once ("SBS4") since a second COSMIC tag would be redundant."""
    return f"{s} | {_ess_cosmic_short(s)}" if s.startswith("eSS") else s


def _volcano_label_etio(s):
    """Point label using the eSS etiology (dominant training mutagen / chemical group)
    instead of the COSMIC match, e.g. "eSS18 | 4NQO", "eSS9 | AA/Dietary"."""
    return f"{s} | {ess_umbrella_tag(s)}" if s.startswith("eSS") else s


def volcano(d, xcol, title, outfile, xlab, xclip=None, label_fn=None, label_side=None,
            bold_novel=True, label_dx=None, transparent=False):
    # label_side: optional {signature: "right"|"left"|"center"|"topright"} to override the
    # automatic label placement for specific points.
    # label_dx: optional {signature: points} extra horizontal nudge for specific labels.
    # bold_novel: bold the label of eSS with no COSMIC match (off for etiology-labelled plots).
    label_fn = label_fn or _volcano_label
    v = d.dropna(subset=[xcol, "q"]).copy()
    if v.empty:
        return
    v["y"] = -np.log10(v["q"].clip(lower=1e-300))
    thr = 0.05 if xcol == "r" else 0.2
    v["assoc"] = (v["q"] < Q) & (v[xcol].abs() >= thr)
    v["x"] = v[xcol].clip(-xclip, xclip) if xclip else v[xcol]
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    ns = v[~v.assoc]
    ax.scatter(ns.x, ns.y, s=12, c=col_g, alpha=0.6, edgecolor="none")
    for cond, c in [((v.assoc) & (v[xcol] > 0), col_sm), ((v.assoc) & (v[xcol] < 0), col_ns)]:
        dd = v[cond]
        ax.scatter(dd.x, dd.y, s=40, c=c, edgecolor="black", linewidth=0.4, zorder=3)
    ax.axhline(-np.log10(Q), ls="--", lw=0.7, color="0.5")
    ax.axvline(0, lw=0.5, color="0.3")
    # tight x-limits based on the UNCLIPPED points, so quasi-separation artefacts parked
    # at +/-xclip (huge ORs, non-significant) don't stretch the axis into empty space
    ref = v[v[xcol].abs() < xclip] if xclip else v
    if ref.empty:
        ref = v
    lo, hi = ref.x.min(), ref.x.max()
    xr = max(hi - lo, 1e-6)
    ax.set_xlim(lo - 0.05 * xr - 0.2, hi + 0.05 * xr + 0.2)
    ax.set_ylim(-0.02 * v.y.max(), v.y.max() * 1.12)
    # labels sit a few points ABOVE their marker (clear of the bubble, not off to the right
    # inflating the margin); points near x=0 anchor to their outer side so the label does not
    # cross the x=0 axis; adjustText only declutters overlaps vertically.
    from matplotlib.transforms import offset_copy
    sides = label_side or {}
    nudges = label_dx or {}
    assoc = v[v.assoc]
    texts = []
    for _, r in assoc.iterrows():
        side = sides.get(r.signature)
        if side == "topright":         # snug at the upper-right corner of the marker
            ha, va, dx, dy = "left", "bottom", 3, 3
        elif side == "right":          # beside the marker, to its right
            ha, va, dx, dy = "left", "center", 6, 0
        elif side == "left":
            ha, va, dx, dy = "right", "center", -6, 0
        elif side in ("center", "top"):
            ha, va, dx, dy = "center", "bottom", 0, 4
        elif abs(r.x) < 1.2:           # auto: keep near-zero labels off the x=0 axis
            ha, va, dx, dy = ("left", "bottom", 3, 4) if r.x >= 0 else ("right", "bottom", -3, 4)
        else:                          # auto: centred above the marker
            ha, va, dx, dy = "center", "bottom", 0, 4
        dx += nudges.get(r.signature, 0)
        off = offset_copy(ax.transData, fig=fig, x=dx, y=dy, units="points")
        is_novel = r.signature.startswith("eSS") and _ess_cosmic_short(r.signature) == "novel"
        texts.append(ax.text(r.x, r.y, label_fn(r.signature), fontsize=6, zorder=5,
                             ha=ha, va=va, transform=off,
                             fontweight="bold" if (bold_novel and is_novel) else "normal"))
    if texts:
        try:
            from adjustText import adjust_text   # no leader lines (arrowprops omitted)
            adjust_text(texts, ax=ax, expand=(1.2, 1.4), only_move={"text": "y"})
        except Exception:
            pass  # labels still render at their points if adjustText is unavailable
    ax.set_xlabel(xlab, fontsize=9)
    ax.set_ylabel(r"$-\log_{10}$(q-value)", fontsize=9)
    ax.set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(outfile, bbox_inches="tight", transparent=transparent)
    plt.close()


def refdr(d):
    e = d[d.signature.str.startswith("eSS")].copy()
    pv = e["p"] if "p" in e else e["q"]
    ok = pv.notna()
    e["q"] = np.nan
    e.loc[ok, "q"] = multipletests(pv[ok], method="fdr_bh")[1]
    return e


def smoking_block(df, sigs, outdir, tag):
    ess_sigs = [s for s in sigs if s.startswith("eSS")]
    # pooled univariate
    tot = df[sigs].sum(axis=1).replace(0, np.nan)
    rel = df[sigs].div(tot, axis=0).fillna(0.0)
    smk = df.smoker.values
    uni = []
    for s in sigs:
        a, b = rel.loc[df.smoker == 1, s].values, rel.loc[df.smoker == 0, s].values
        try:
            p = mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            p = np.nan
        r = pointbiserialr(smk, rel[s].values)[0] if rel[s].std() > 0 else 0.0
        uni.append({"signature": s, "r": r, "p": p})
    uni = pd.DataFrame(uni)
    uni["q"] = multipletests(uni.p.fillna(1), method="fdr_bh")[1]
    uni.sort_values("q").to_csv(os.path.join(outdir, "volcano_pooled_univariate.csv"), index=False)
    # cancer-type-adjusted
    det = (df[sigs] > 0).astype(int)
    adj = []
    base = df[["smoker", "age", "female", "ctype"]].copy()
    for s in sigs:
        dd = base.copy()
        dd["y"] = det[s].values
        dd = dd.dropna(subset=["age", "female"])
        keep = [t for t, gg in dd.groupby("ctype") if gg.y.nunique() == 2]
        dd = dd[dd.ctype.isin(keep)]
        if dd.y.nunique() < 2 or dd.smoker.nunique() < 2 or len(dd) < 50:
            adj.append({"signature": s, "log2_OR": np.nan, "p": np.nan})
            continue
        try:
            m = smf.logit("y ~ smoker + age + female + C(ctype)", data=dd).fit(disp=0, maxiter=200)
            adj.append({"signature": s, "log2_OR": m.params["smoker"] / np.log(2), "p": m.pvalues["smoker"]})
        except Exception:
            adj.append({"signature": s, "log2_OR": np.nan, "p": np.nan})
    adj = pd.DataFrame(adj)
    ok = adj.p.notna()
    adj["q"] = np.nan
    adj.loc[ok, "q"] = multipletests(adj.loc[ok, "p"], method="fdr_bh")[1]
    adj.sort_values("q").to_csv(os.path.join(outdir, "volcano_cancertype_adjusted.csv"), index=False)
    uni_e, adj_e = refdr(uni), refdr(adj)
    volcano(uni, "r", f"Smoking association, ALL signatures ({tag})\nPOOLED / univariate — tissue-CONFOUNDED overview",
            os.path.join(outdir, "volcano_pooled_univariate.pdf"), "point-biserial r (<- non-smoker | smoker ->)")
    volcano(adj, "log2_OR", f"Smoking association, cancer-type-ADJUSTED ({tag})\npresence ~ smoker + age + sex + C(cancer type)",
            os.path.join(outdir, "volcano_cancertype_adjusted.pdf"),
            "log2 OR ever vs never, adjusted (<- non-smoker | smoker ->)", xclip=6.0)
    volcano(uni_e, "r", f"Smoking association, eSS ONLY ({tag})\nPOOLED / univariate — tissue-CONFOUNDED overview",
            os.path.join(outdir, "volcano_pooled_univariate_eSS.pdf"), "point-biserial r (<- non-smoker | smoker ->)")
    volcano(adj_e, "log2_OR", f"Smoking association, eSS ONLY, cancer-type-ADJUSTED ({tag})\npresence ~ smoker + age + sex + C(cancer type)",
            os.path.join(outdir, "volcano_cancertype_adjusted_eSS.pdf"),
            "log2 OR ever vs never, adjusted (<- non-smoker | smoker ->)", xclip=6.0)
    uni_e.sort_values("q").to_csv(os.path.join(outdir, "volcano_pooled_univariate_eSS.csv"), index=False)
    adj_e.sort_values("q").to_csv(os.path.join(outdir, "volcano_cancertype_adjusted_eSS.csv"), index=False)

    # eSS-TMB per cancer beeswarm (smoking-script style)
    df2 = df.copy()
    df2["eSS_perMb"] = df2[ess_sigs].sum(axis=1) / MB
    d1 = df2[df2.eSS_perMb > 0]
    groups = sorted([(t, np.log10(g.eSS_perMb.values)) for t, g in d1.groupby("ctype")],
                    key=lambda x: np.median(x[1]))
    if groups:
        fig, ax = plt.subplots(figsize=(0.5 * len(groups) + 2, 5))
        for i, (t, vals) in enumerate(groups):
            vals = np.sort(vals)
            xs = np.linspace(i - 0.35, i + 0.35, len(vals)) if len(vals) > 1 else [i]
            ax.scatter(xs, vals, s=5, c="black", alpha=0.5, edgecolor="none")
            ax.plot([i - 0.35, i + 0.35], [np.median(vals)] * 2, color="red", lw=2)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([t for t, _ in groups], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("log10 eSS mutations / Mb", fontsize=10)
        ax.set_title(f"eSS TMB per cancer type — {tag}", fontsize=11)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "eSS_TMB_per_cancer.pdf"), bbox_inches="tight")
        plt.close()


# ================================================================== Fisher (fisher/ folder, per-cancer VOLCANO)
def per_tissue_stats(df, tissues, sig_list):
    brows, prows = [], []
    for t in tissues:
        g = df[df.ctype == t]
        ev, nv = g[g.smoker == 1], g[g.smoker == 0]
        for s in sig_list:
            a, b = ev[s].values / MB, nv[s].values / MB
            if (a > 0).sum() + (b > 0).sum() >= 5:
                brows.append({"signature": s, "tissue": t, "wilcoxon_p": ranksums(a, b).pvalue,
                              "fold": (a.mean() + 1e-9) / (b.mean() + 1e-9)})
            ce, cn = int((ev[s] > 0).sum()), int((nv[s] > 0).sum())
            if ce + cn >= 5:
                orr, pf = fisher_exact([[ce, len(ev) - ce], [cn, len(nv) - cn]])
                prows.append({"signature": s, "tissue": t, "fisher_p": pf, "OR": orr,
                              "prev_ever": ce / len(ev), "prev_never": cn / len(nv)})
    bb, pp = pd.DataFrame(brows), pd.DataFrame(prows)
    if not bb.empty:
        bb["q"] = multipletests(bb.wilcoxon_p.fillna(1), method="fdr_bh")[1]
    if not pp.empty:
        pp["q"] = multipletests(pp.fisher_p.fillna(1), method="fdr_bh")[1]
    return bb, pp


def _volcano_panel(ax, sub, xcol, title):
    """One per-cancer volcano panel: x=log2(effect), y=-log10(q within tissue)."""
    sub = sub.dropna(subset=[xcol, "q"]).copy()
    if sub.empty:
        ax.axis("off")
        ax.set_title(title, fontsize=9)
        return
    x = np.log2(sub[xcol].clip(lower=1e-2, upper=100.0).values)
    x = np.clip(x, -6, 6)
    y = -np.log10(sub["q"].clip(lower=1e-300).values)
    sig = sub["q"].values < Q
    is_ess = sub["signature"].str.startswith("eSS").values
    ax.scatter(x[~sig], y[~sig], s=16, c=col_g, alpha=0.6, edgecolor="none")
    up = sig & (x > 0)
    dn = sig & (x < 0)
    ax.scatter(x[up], y[up], s=42, c=col_sm, edgecolor="black", linewidth=0.4, zorder=3)
    ax.scatter(x[dn], y[dn], s=42, c=col_ns, edgecolor="black", linewidth=0.4, zorder=3)
    for k, (xi, yi, s, e, sg) in enumerate(zip(x, y, sub["signature"].values, is_ess, sig)):
        if sg:
            ax.annotate(f"{s} {lab(s)}" if e else s, (xi, yi), fontsize=5,
                        ha="left" if xi >= 0 else "right", xytext=(4 if xi >= 0 else -4, 3 if k % 2 == 0 else -6),
                        textcoords="offset points", fontweight="bold" if (e and lab(s) == "novel") else "normal")
    ax.axhline(-np.log10(Q), ls="--", lw=0.7, color="0.5")
    ax.axvline(0, lw=0.5, color="0.3")
    ax.set_xlim(-6.5, 6.5)
    ax.set_title(title, fontsize=9)


def fisher_block(df, sigs, outdir, tag):
    """Fisher prevalence + Wilcoxon burden -> fisher/ folder; per-cancer VOLCANOES."""
    fdir = os.path.join(outdir, "fisher")
    pcdir = os.path.join(fdir, "per_cancer_volcano")
    os.makedirs(pcdir, exist_ok=True)
    ess_sigs = [s for s in sigs if s.startswith("eSS")]
    tissues = sorted([t for t, g in df.groupby("ctype")
                      if (g.smoker == 1).sum() >= MIN_N and (g.smoker == 0).sum() >= MIN_N])
    bdf, pdf_ = per_tissue_stats(df, tissues, sigs)
    bdf_e, pdf_e = per_tissue_stats(df, tissues, ess_sigs)
    for d_, name in [(pdf_, "prevalence_fisher.csv"), (bdf, "burden_wilcoxon.csv"),
                     (pdf_e, "prevalence_fisher_eSS.csv"), (bdf_e, "burden_wilcoxon_eSS.csv")]:
        (d_.sort_values("q") if not d_.empty else d_).to_csv(os.path.join(fdir, name), index=False)

    # --- per-cancer VOLCANO (one PDF per tissue) + a combined faceted overview
    def all_tissue_volcanoes(pv, xcol, kind):
        if pv.empty:
            return
        ts = sorted(pv.tissue.unique())
        # individual PDFs
        for t in ts:
            fig, ax = plt.subplots(figsize=(5.2, 4.4))
            _volcano_panel(ax, pv[pv.tissue == t], xcol, f"{t} — {kind} ever vs never ({tag})")
            ax.set_xlabel(f"log2 {xcol} (<- never | ever ->)", fontsize=9)
            ax.set_ylabel("-log10(q within tissue)", fontsize=9)
            fig.tight_layout()
            fig.savefig(os.path.join(pcdir, f"{kind}_{t}.pdf"), bbox_inches="tight")
            plt.close()
        # faceted overview
        n = len(ts)
        ncol = min(4, n)
        nrow = int(np.ceil(n / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.0 * nrow), squeeze=False)
        for i, t in enumerate(ts):
            _volcano_panel(axes[i // ncol][i % ncol], pv[pv.tissue == t], xcol, t)
        for j in range(n, nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        kind_lab = {"prevalence": "prevalence (Fisher OR)", "burden": "burden (Wilcoxon fold)"}[kind]
        fig.suptitle(f"Per-cancer {kind_lab} — ever vs never ({tag})\n"
                     f"x = log2 {xcol}  (<- never | ever ->),  y = -log10(q within tissue)", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        fig.savefig(os.path.join(fdir, f"{kind}_volcano_by_cancer.pdf"), bbox_inches="tight")
        plt.close()

    all_tissue_volcanoes(pdf_, "OR", "prevalence")
    all_tissue_volcanoes(bdf, "fold", "burden")
    return len(tissues), (0 if pdf_.empty else int((pdf_.q < Q).sum())), (0 if bdf.empty else int((bdf.q < Q).sum()))


# ================================================================== driver
def run(name, act_df, outdir):
    ess_cols = [c for c in act_df.columns if c.startswith("eSS")]
    sigs = [c for c in act_df.columns if c != "Samples"]
    # TMB df (needs CancerType from metadata)
    ct = META[["Samples", "ctype"]].rename(columns={"ctype": "CancerType"})
    tdf = act_df.merge(ct, on="Samples", how="inner").dropna(subset=["CancerType"]).reset_index(drop=True)
    # smoking/fisher df (needs smoker/age/female/ctype)
    sdf = act_df.merge(META, on="Samples", how="inner").dropna(subset=["smoker", "ctype"]).reset_index(drop=True)
    # GUARD: these analyses are patient-level. If the "Samples" are signatures
    # (COSMIC / de-novo decomposition), almost none join the clinical metadata.
    if len(tdf) < 50:
        print(f"\n=== {name} ===\n  SKIPPED: only {len(tdf)}/{len(act_df)} rows join clinical metadata "
              f"-> this is a signature-level decomposition, not a patient-level assignment; "
              f"TMB/smoking/Fisher do not apply.")
        return
    os.makedirs(outdir, exist_ok=True)
    sdf["smoker"] = sdf["smoker"].astype(int)
    n_sm, n_ns = int((sdf.smoker == 1).sum()), int((sdf.smoker == 0).sum())
    print(f"\n=== {name} ===")
    print(f"  activities: {len(act_df)} samples x {len(sigs)} sigs ({len(ess_cols)} eSS)")
    print(f"  TMB df (w/ cancer type): {len(tdf)} | smoking df (w/ smoking): {len(sdf)} "
          f"(ever={n_sm} never={n_ns})")
    tmb_block(tdf, ess_cols, outdir)
    smoking_block(sdf, sigs, outdir, name)
    nt, npv, nbd = fisher_block(sdf, sigs, outdir, name)
    print(f"  fisher/: {nt} tissues, {npv} sig prevalence, {nbd} sig burden (FDR<{Q})")
    print(f"  -> {outdir}")


# Only PATIENT-LEVEL SPA assignments are valid inputs for TMB/smoking/Fisher.
# COSMIC_v3.6_Decomposition_eSS51 (rows = 101 COSMIC SBS signatures) and
# SPE_DeNovo_eSS51 (rows = per-cancer de-novo signatures, SBS96A...) are
# SIGNATURE-level decompositions with no patients, so they are excluded here.
RUNS = [
    ("samples_as_input_eSS_COSMIC", ("dir",
     os.path.join(BASE, "Collinearity_curation_eSS_COSMIC", "outputs_samples_as_input_ess5051")),
     os.path.join(BASE, "Collinearity_curation_eSS_COSMIC", "outputs_samples_as_input_ess5051", "analysis")),
    ("unsupervised_hybrid", ("file",
     os.path.join(BASE, "Collinearity_curation_eSS_COSMIC", "unsupervised_SPA_run", "Unsupervised_Activities.txt")),
     os.path.join(BASE, "Collinearity_curation_eSS_COSMIC", "unsupervised_SPA_run", "analysis")),
]

if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    for name, (kind, src), outdir in RUNS:
        if only and only not in name:
            continue
        act = load_single(src) if kind == "file" else load_pooled(src)
        run(name, act, outdir)
    print("\nAll done.")
