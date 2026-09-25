import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

SUMMARY = '/tscc/lustre/restricted/alexandrov-ddn/users/z8jiang/signature_stress_test/Maria_Z_eSignature/null_SPA_summary.csv'
OUT = '/tscc/lustre/restricted/alexandrov-ddn/users/z8jiang/signature_stress_test/Maria_Z_eSignature'

# LABEL_ALL = False labels label_these only. True labels all 101 targets, which
LABEL_ALL = True

d = pd.read_csv(SUMMARY)

PASS = 'tab:blue'
FAIL = 'silver'
MID = 'tab:orange'

if LABEL_ALL:
    FIGSIZE, LABEL_FONT = (10.5, 10.5), 6.5
else:
    FIGSIZE, LABEL_FONT = (8.2, 8.2), 7.5

fig, ax = plt.subplots(figsize=FIGSIZE)

# ----------------------------------------------------------------- panel A
both = d['beats_chance_fdr5'] & (d['observed_cosine'] >= 0.90)
chance_only = d['beats_chance_fdr5'] & (d['observed_cosine'] < 0.90)
neither = ~d['beats_chance_fdr5']

ax.plot([0.0, 1.05], [0.0, 1.05], color='k', lw=1.0, zorder=1)
ax.axhline(0.90, color='k', ls=':', lw=0.9, zorder=1)
ax.axvline(0.90, color='k', ls=':', lw=0.9, zorder=1)
ax.fill_between([0.0, 1.05], [0.0, 1.05], 1.05, color='whitesmoke', zorder=0)

for mask, colour, lab in ((neither, FAIL, 'Does not exceed random chance reconst cos'),
                          (chance_only, MID, 'Exceed chance, cos < 0.90'),
                          (both, PASS, 'Exceed chance, cos ≥ 0.90')):
    ax.scatter(d.loc[mask, 'null_p95'], d.loc[mask, 'observed_cosine'],
               s=34, c=colour, edgecolor='white', linewidth=0.5,
               label=lab, zorder=3)

ax.set_xlabel('Permutation reconstruction cosine ceiling\n(95th percentile of 1000 shuffled eSS libraries)')
ax.set_ylabel('Observed reconstruction cosine\n(using real eSS library)')
# ax.set_title('A   every target carries its own chance ceiling', loc='left',
#              fontsize=11, fontweight='bold')

# The limits are set before adjust_text runs, because the solver keeps labels
# inside the axes and reads the final limits to do that. Both axes share the
# same 0 to 1.05 range and the same 0.1 tick spacing, so 0.9 is readable on
# each one and the diagonal is a true 45 degree line.
ax.set_xlim(0.1, 1.05)
ax.set_ylim(0.1, 1.05)
ax.set_xticks(np.arange(0.1, 1.01, 0.1))
ax.set_yticks(np.arange(0.1, 1.01, 0.1))
# ax.set_aspect('equal')

# The legend moves to the upper left, which holds no data points. In the lower
# left it lands on SBS17b and SBS60, in the lower right on SBS48.
ax.legend(loc='upper left', fontsize=8, frameon=False)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)

label_these = ['SBS88', 'SBS90', 'SBS98', 'SBS47', 'SBS16', 'SBS51',
               'SBS110', 'SBS32', 'SBS46', 'SBS48', 'SBS28', 'SBS17b']
to_label = d if LABEL_ALL else d[d['cosmic_signature'].isin(label_these)]
texts = [ax.text(r['null_p95'], r['observed_cosine'], r['cosmic_signature'],
                 fontsize=LABEL_FONT, zorder=5)
         for _, r in to_label.iterrows()]

# adjustText is the Python port of ggrepel. Tested against version 1.4.0. The
# 0.8 releases use different argument names, so this call needs rewriting if an
# older build loads.
adjust_text(texts, ax=ax,
            x=to_label['null_p95'].values,
            y=to_label['observed_cosine'].values,
            expand=(1.12, 1.30),
            force_text=(0.6, 0.9),
            force_static=(0.25, 0.4),
            force_explode=(0.35, 0.7),
            explode_radius=18,
            max_move=(40, 40),
            iter_lim=800,
            min_arrow_len=4,
            arrowprops=dict(arrowstyle='-', color='0.55', lw=0.45,
                            shrinkA=1, shrinkB=2))

fig.tight_layout()
fig.savefig(f'{OUT}/panels.png', dpi=600)
print('wrote panels.png')

# ------------------------------------------------- second figure, the forest
sub = pd.concat([d.nlargest(18, 'margin_over_null_p95'),
                 d.nsmallest(18, 'margin_over_null_p95')])
sub = sub.sort_values('margin_over_null_p95')
y = np.arange(len(sub))

fig2, ax = plt.subplots(figsize=(7.4, 9.4))
ax.hlines(y, sub['null_median'], sub['null_max'], color='#d9d9d9', lw=5,
          zorder=1, label='chance range, median to max')
ax.hlines(y, sub['null_median'], sub['null_p95'], color='#9e9e9e', lw=5,
          zorder=2, label='chance range, median to 95th percentile')
cols = np.where(sub['beats_chance_fdr5'] & (sub['observed_cosine'] >= 0.90), PASS,
                np.where(sub['beats_chance_fdr5'], MID, FAIL))
ax.scatter(sub['observed_cosine'], y, s=52, c=cols, edgecolor='k',
           linewidth=0.6, zorder=4, label='observed, real library')
ax.axvline(0.90, color='k', ls=':', lw=0.9)
ax.axhline(len(sub) / 2 - 0.5, color='k', lw=0.7)

ax.set_yticks(y)
ax.set_yticklabels(sub['cosmic_signature'], fontsize=8)
ax.set_xlabel('best achievable reconstruction cosine')
ax.set_title('18 largest and 18 smallest margins over the chance ceiling',
             fontsize=11, fontweight='bold', loc='left')
ax.set_xlim(0.05, 1.03)
ax.legend(loc='upper left', fontsize=7.5, frameon=False)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
fig2.tight_layout()
fig2.savefig(f'{OUT}/forest.png', dpi=600)
print('wrote forest.png')