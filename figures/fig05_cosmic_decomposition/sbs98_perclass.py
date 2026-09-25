#!/usr/bin/env python3
"""SBS98 C>G rescue: is the match to eSS1's C>G block above chance and stable?

SBS98's whole-profile reconstruction is 0.845 (below 0.90); the manuscript rescues
it by noting its C>G block resembles eSS1's C>G. This tests that specific claim.

Permutation p : permute eSS1's 16 C>G channels; how often does a permuted eSS1 C>G
                match SBS98's C>G as well as the real one?  (specificity of the match)
Bootstrap CI  : multinomial resample of eSS1's C>G counts (total eSS1 = 217,761
                mutations) -> 95% CI on cosine(SBS98 C>G, eSS1 C>G).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ESS1_TOTAL_MUTATIONS = 217761
B = 10000
SEED = 20260902

ess = pd.read_csv(os.path.join(HERE, "input_eSS_reference.txt"), sep="\t").set_index("MutationType")
cosmic = pd.read_csv(os.path.join(HERE, "input_COSMIC_Human_SBS-96_GRCh38_v3.6.tsv"),
                     sep="\t").set_index("Type").loc[ess.index]
order = list(ess.index)
CG = np.array([i for i, mt in enumerate(order) if mt[2:5] == "C>G"])


def cosine(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


ess1 = (ess["eSS1"] / ess["eSS1"].sum()).to_numpy(float)
sbs98 = (cosmic["SBS98"] / cosmic["SBS98"].sum()).to_numpy(float)
ess1_cg = ess1[CG]
sbs98_cg = sbs98[CG]
observed = cosine(sbs98_cg, ess1_cg)

# specificity: eSS1 vs every eSS on the C>G block
sims = sorted(((c, cosine(sbs98_cg, (ess[c] / ess[c].sum()).to_numpy(float)[CG]))
               for c in ess.columns), key=lambda t: -t[1])

# permutation null: shuffle eSS1's C>G channels
rng = np.random.default_rng(SEED)
perm_null = np.array([cosine(sbs98_cg, ess1_cg[rng.permutation(len(CG))]) for _ in range(B)])
perm_p = (1 + int(np.sum(perm_null >= observed))) / (B + 1)

# multinomial bootstrap: eSS1 C>G counts (proportion x total eSS1 mutations)
cg_counts = np.round(ess1_cg * ESS1_TOTAL_MUTATIONS).astype(int)
n_cg = int(cg_counts.sum())
p_cg = cg_counts / cg_counts.sum()
boot = np.empty(B)
for i in range(B):
    resampled = rng.multinomial(n_cg, p_cg)
    boot[i] = cosine(sbs98_cg, resampled.astype(float))
ci = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))

print(f"SBS98 C>G vs eSS1 C>G: observed cosine = {observed:.3f}")
print(f"  next-best eSS on C>G: {sims[1][0]} = {sims[1][1]:.3f}  (eSS1 is #{[c for c,_ in sims].index('eSS1')+1})")
print(f"  permutation p (eSS1 C>G channels shuffled, B={B}) = {perm_p:.4f}")
print(f"  bootstrap 95% CI (multinomial, n_C>G={n_cg}) = [{ci[0]:.3f}, {ci[1]:.3f}]")

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
a1.hist(perm_null, bins=40, color="#8c8cc0", alpha=0.8)
a1.axvline(observed, color="#b3202c", lw=2.4)
a1.set_title(f"Permutation null (eSS1 C>G shuffled)\nobserved {observed:.2f}, p={perm_p:.4f}", fontsize=10)
a1.set_xlabel("cosine(SBS98 C>G, permuted eSS1 C>G)", fontsize=9)
a1.set_ylabel(f"permutations (n={B})", fontsize=9)
a1.set_xlim(0, 1)
a2.hist(boot, bins=40, color="#7fb069", alpha=0.8)
a2.axvline(ci[0], color="#333", ls="--", lw=1.2)
a2.axvline(ci[1], color="#333", ls="--", lw=1.2)
a2.axvline(observed, color="#b3202c", lw=2.4)
a2.set_title(f"Bootstrap of eSS1 C>G\n95% CI [{ci[0]:.2f}, {ci[1]:.2f}]", fontsize=10)
a2.set_xlabel("cosine(SBS98 C>G, resampled eSS1 C>G)", fontsize=9)
for ax in (a1, a2):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
fig.suptitle("SBS98 C>G rescue vs eSS1", fontsize=12)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, f"fig_SBS98_perclass_rescue.{ext}"), dpi=300, bbox_inches="tight")

pd.DataFrame([dict(comparison="SBS98_CG_vs_eSS1_CG", observed_cosine=round(observed, 4),
                   next_best_eSS=sims[1][0], next_best_cosine=round(sims[1][1], 4),
                   permutation_p=round(perm_p, 4), boot_ci_low=round(ci[0], 4),
                   boot_ci_high=round(ci[1], 4))]).to_csv(
    os.path.join(HERE, "results_SBS98_perclass_rescue.csv"), index=False)
