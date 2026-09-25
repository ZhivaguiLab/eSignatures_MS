# `results_opportunity_normalized` vs `results_default`

**Method:** per-sample opportunity normalization
([pipeline/normalize_by_own_opportunity.py](../pipeline/normalize_by_own_opportunity.py)),
branch `per-sample-opportunity-normalization`. Every human, mouse, and rat
sample — WGS *and* WES — is normalized by its own species+build+sequencing
-technology trinucleotide opportunity table, independently:

```
rate[c]    = raw_count[c] / opportunity[c]   (that sample's own context-count table)
profile[c] = rate[c] / sum(rate)              (renormalize to sum 1)
```

No common reference genome is chosen; each sample only ever uses its own
opportunity table (human WGS → GRCh38-WGS; human WES → GRCh38-exome; mouse
WGS → mm10-WGS; mouse WES → mm10-exome; rat WGS → rn7-WGS). Chicken and
celegans have no confirmed opportunity table yet, so their input files are
left completely unchanged (still plain `raw / raw.sum()`).

This is a different method from the liftover/translation approach on
`wes-to-wgs-normalization` ([normalize_wes_to_wgs.py](../pipeline/normalize_wes_to_wgs.py)),
not a refinement of it — see that branch's
[COMPARISON_SUMMARY.md](../results_wes_to_wgs_normalized/COMPARISON_SUMMARY.md)
for the liftover version's results. Raw counts are unchanged in both
approaches; only the normalized profile used for clustering/consensus
averaging differs.

## Scope of what changed

This method touches every human/mouse/rat sample — 570 of ~902 total in the
main atlas — versus the liftover approach's 152 (WES-only). That's the main
reason the effect below is substantially larger.

## Headline numbers

| | Default | Per-sample opportunity | Δ |
|---|---|---|---|
| Main clusters | 49 | 52 | +3 |
| Main-cluster samples | 508 | 503 | −5 |
| Small-cluster samples | 32 | 30 | −2 |
| Singletons | 130 | 137 | +7 |
| COSMIC-matched clusters (≥0.85) | 26 | **14** | **−12** |
| De novo clusters | 23 | 38 | +15 |

## Cluster membership

- 22 samples left main-cluster status, 17 newly entered — split roughly
  Mouse (13/12) and Human (9/5). **Zero rat, chicken, or celegans samples
  affected**: chicken/celegans untouched by design, and rat's 5 WGS samples,
  despite being corrected, weren't perturbed enough to flip any cluster
  boundary.
- Among the 486 samples that stayed "main" in both runs, **241 regrouped
  with different cluster partners** (194 mouse, 47 human, 0 other species) —
  a much larger reshuffling than the liftover approach's 4 affected
  groupings.

## The COSMIC-match collapse (26 → 14): why this is expected, not a bug

Under `equal_replicate` averaging, a cluster's consensus profile is the
average of its members' *already-normalized* profiles. In the liftover
approach, every sample was normalized onto the same target (WGS) basis
first, so a cluster's consensus stays on one consistent basis, comparable
to COSMIC. In this method, each sample is normalized by *its own*
opportunity table — so a cluster mixing human-WGS, human-WES, and
mouse-WGS samples now averages together profiles from **three different
opportunity bases at once**, and the resulting consensus isn't on any
single well-defined basis anymore, COSMIC's included. Since COSMIC
signatures are built on a specific reference (standard human genome-wide),
a consensus that's an average across mismatched bases would plausibly
drift away from any COSMIC signature's shape — independent of whether the
underlying biology is being represented more or less correctly.

This tracks with using a genuinely different opportunity table per species
*and* per genome build/technology: clusters with mixed species/genome/tech
composition are the ones expected to be most affected, since they're the
ones averaging across incompatible bases. Clusters that are internally
homogeneous (single species, single tech) shouldn't show this effect, since
all their members share one opportunity table already.

**Open follow-up:** confirm this mechanism directly by checking whether the
clusters that lost their COSMIC match are disproportionately the
mixed-composition ones, versus single-species/single-tech clusters losing
matches too (which would point to something other than basis-mixing).
