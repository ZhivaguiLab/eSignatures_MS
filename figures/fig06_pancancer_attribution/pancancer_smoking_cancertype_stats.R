# ============================================================
# Script: pancancer_smoking_cancertype_stats.R
# Purpose: Reproduce the pan-cancer eSS smoking and cancer-type
#          association statistics reported in MS1 (Figure 6c-d).
# Input:   pooled per-sample eSS activity matrix + clinical metadata
# Output:  results_smoking_association.csv, results_percancer_association.csv
# Author:  M. Zhivagui
# Date:    2026-07-26
# ============================================================
#
# This script re-derives the smoking-association statistics from the
# de novo pan-cancer decomposition. The publication figures were produced
# by the Python pipeline (run2_*.py); this R version reproduces the identical
# models so the numbers can be regenerated in a clean R session with CRAN
# packages only.
#
# Reported result (cancer-type + cohort adjusted logistic regression):
#   eSS18 (SBS100)  log2 OR = 4.21   q = 9.7e-21
#   eSS47 (MMS/DMS) log2 OR = 3.44   q = 5.7e-4
#   eSS17 (SBS4)    log2 OR = 1.67   q = 1.4e-3
#   eSS33 (SBS22a)  log2 OR = 0.87   q = 1.3e-3
# Per-cancer localisation: eSS18 in lung and head/neck; eSS33 in lung;
#   eSS47 in head/neck.

# --- 1. LOAD LIBRARIES ---
library(readr)   # fast, type-safe CSV/TSV reading
library(dplyr)   # data manipulation and grouping
library(tidyr)   # reshaping the activity matrix

# --- 2. SET PATHS ---
# Set DATA_DIR to the folder holding the two input files (see README):
#   pooled_Decompose_Solution_Activities.txt  pooled per-cancer eSS activities (samples x signatures)
#   HD Signatures WGS Datasets v2.0 - final_sample_summary_v2.tsv  clinical metadata
# The pooled activities file is built by concatenating the per-cancer
# Decompose_Solution_Activities.txt tables from the released data archive.
DATA_DIR   <- "inputs"
ACT_FILE   <- file.path(DATA_DIR, "pooled_Decompose_Solution_Activities.txt")
META_FILE  <- file.path(DATA_DIR, "HD Signatures WGS Datasets v2.0 - final_sample_summary_v2.tsv")
OUT_DIR    <- "results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Analysis constants. MB = 2800 because every cohort here is whole-genome.
MB          <- 2800
MIN_N       <- 20      # minimum ever/never per cancer type for the per-cancer model
Q_THR       <- 0.05    # BH FDR significance threshold
DETECT_FRAC <- 0.05    # a sample "carries" an eSS at >= 5% of its total burden

# --- 3. LOAD DATA ---
# Why: the activity file is a matrix whose first column is the sample ID and
# whose remaining columns are per-signature mutation counts.
stopifnot(file.exists(ACT_FILE), file.exists(META_FILE))
act  <- read_tsv(ACT_FILE, show_col_types = FALSE)
colnames(act)[1] <- "Sample_ID"

meta <- read_tsv(META_FILE, show_col_types = FALSE) %>%
  filter(Sample_ID != "Sample_ID") %>%              # drop any repeated header row
  transmute(
    Sample_ID,
    smoker = case_when(                              # classify smoking status
      grepl("ever|current|former|yes", Smoking_Status, ignore.case = TRUE) ~ 1L,
      grepl("never|no",               Smoking_Status, ignore.case = TRUE) ~ 0L,
      TRUE ~ NA_integer_),
    age    = suppressWarnings(as.numeric(Age)),
    female = ifelse(Sex == "Female", 1L, ifelse(Sex == "Male", 0L, NA_integer_)),
    cohort = Cohort,
    ctype  = Cancer_Type_HD_Grouping)

# --- 4. ANALYSIS ---
set.seed(42)

# 4a. Convert raw activities to per-sample fraction of total burden, then to a
#     binary carrier call at the 5%-of-total-burden threshold.
sig_cols <- setdiff(colnames(act), "Sample_ID")
act_frac <- act %>%
  mutate(total = rowSums(across(all_of(sig_cols)))) %>%
  filter(total > 0)
carrier <- act_frac %>%
  mutate(across(all_of(sig_cols), ~ as.integer(.x / total >= DETECT_FRAC))) %>%
  select(Sample_ID, all_of(sig_cols))

# Restrict to eSS columns (drop endogenous COSMIC / SBS96N background columns).
ess_cols <- grep("^eSS", sig_cols, value = TRUE)

# 4b. Join carriers to clinical metadata; keep smoking-annotated samples.
dat <- carrier %>%
  inner_join(meta, by = "Sample_ID") %>%
  filter(!is.na(smoker), !is.na(ctype))
message(sprintf("Smoking-annotated samples: %d (ever=%d, never=%d)",
                nrow(dat), sum(dat$smoker == 1), sum(dat$smoker == 0)))

# 4c. Pan-cancer, cancer-type + cohort-adjusted logistic regression per eSS.
#     Why: cancer type and cohort confound smoking; adjusting isolates the
#     within-tissue smoking effect. log2 OR is the smoker coefficient / ln(2).
fit_one <- function(sig, d) {
  d$y <- d[[sig]]
  # require variation and both smoking groups among carriers
  if (length(unique(d$y)) < 2) return(NULL)
  df <- d %>% filter(!is.na(age), !is.na(female))
  m <- tryCatch(
    glm(y ~ smoker + age + female + factor(ctype) + factor(cohort),
        data = df, family = binomial()),
    error = function(e) NULL)
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients
  if (!"smoker" %in% rownames(co)) return(NULL)
  # Mann-Whitney on the continuous fraction as an unadjusted effect measure
  frac <- act_frac %>% select(Sample_ID, all_of(sig)) %>%
    inner_join(meta %>% select(Sample_ID, smoker), by = "Sample_ID") %>%
    filter(!is.na(smoker))
  mw <- suppressWarnings(wilcox.test(frac[[sig]] ~ frac$smoker))
  data.frame(
    eSS      = sig,
    log2_OR  = co["smoker", "Estimate"] / log(2),
    p_logit  = co["smoker", "Pr(>|z|)"],
    p_MWU    = mw$p.value,
    n_carrier = sum(df$y == 1))
}

res <- bind_rows(lapply(ess_cols, fit_one, d = dat))
res$q_logit <- p.adjust(res$p_logit, method = "BH")   # Benjamini-Hochberg FDR
res <- res %>% arrange(q_logit)
write_csv(res, file.path(OUT_DIR, "results_smoking_association.csv"))

# 4d. Per-cancer-type model: detected ~ smoker + age + sex + cohort, within each
#     tissue with >= MIN_N ever and >= MIN_N never smokers. Localises the signal.
per_cancer <- function(sig, d) {
  out <- list()
  for (ct in unique(d$ctype)) {
    dc <- d %>% filter(ctype == ct, !is.na(age), !is.na(female))
    if (sum(dc$smoker == 1) < MIN_N || sum(dc$smoker == 0) < MIN_N) next
    dc$y <- dc[[sig]]
    if (length(unique(dc$y)) < 2) next
    m <- tryCatch(glm(y ~ smoker + age + female + factor(cohort),
                      data = dc, family = binomial()), error = function(e) NULL)
    if (is.null(m)) next
    co <- summary(m)$coefficients
    if (!"smoker" %in% rownames(co)) next
    out[[length(out)+1]] <- data.frame(
      eSS = sig, cancer_type = ct,
      log2_OR = co["smoker","Estimate"]/log(2),
      p = co["smoker","Pr(>|z|)"])
  }
  bind_rows(out)
}
pc <- bind_rows(lapply(ess_cols, per_cancer, d = dat))
pc$q <- p.adjust(pc$p, method = "BH")
pc <- pc %>% arrange(q)
write_csv(pc, file.path(OUT_DIR, "results_percancer_association.csv"))

# --- 5. SAVE OUTPUT ---
message("Wrote:")
message("  ", file.path(OUT_DIR, "results_smoking_association.csv"))
message("  ", file.path(OUT_DIR, "results_percancer_association.csv"))
message("Top smoker-enriched eSS (q < 0.05):")
print(res %>% filter(q_logit < Q_THR) %>%
        select(eSS, log2_OR, q_logit), row.names = FALSE)
