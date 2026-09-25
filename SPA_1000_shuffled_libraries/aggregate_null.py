import os
import glob
import numpy as np
import pandas as pd

BASE_DIR   = '/tscc/lustre/restricted/alexandrov-ddn/users/z8jiang/signature_stress_test/Maria_Z_eSignature'
RESULT_DIR = os.path.join(BASE_DIR, 'null_cosines')
# Note on the input files:
# We are now parsing the De_Novo_map_to_COSMIC_SBS96.csv files instead of the 
# sample stats files. The sample refitting applies forward stagewise sparsity 
# penalties which ruins the cosine similarities when treating reference signatures 
# as samples (total mutations = 1).
# 
# The mapping file gives us the pure non-negative least squares (NNLS) 
# mathematical decomposition without the biological sample penalties. By running 
# with new_signature_thresh_hold=0.0, everything gets forced to map and we 
# get an unbiased look at the cosines against the null libraries.


# Absolute cosine floors reported side by side. Beating chance just means the fit 
# is better than a shuffled library, but a floor tells us if the fit is actually 
# any good in an absolute sense. Both are needed.
ADEQUACY_FLOORS = [0.85, 0.90, 0.95]



def benjamini_hochberg(pvals):
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted

def load_run(path):
    """
    Reads the collated CSV files written by null_library_fit.py 
    (observed.csv and null_*.csv).
    """
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    
    # Check for the columns output by our null_library_fit.py script
    req_cols = ['cosmic_signature', 'cosine_similarity']
    for col in req_cols:
        if col not in df.columns:
            raise SystemExit(f"Error: {path} is missing the '{col}' column.")
            
    # Clean up the prefix just in case we are reading older files where 
    # the middleman script didn't strip it yet.
    df['cosmic_signature'] = df['cosmic_signature'].astype(str).str.replace('Signature 96-', '').str.strip()
    
    dup = df['cosmic_signature'].duplicated()
    if dup.any():
        print(f"{os.path.basename(path)}: dropped {int(dup.sum())} duplicate rows")
        df = df[~dup]
        
    return df.set_index('cosmic_signature')

def build(files):
    """
    Stacks the observed cosines and the replicates-by-targets null matrix.
    Target order is dictated by the observed file.
    """
    obs_path = os.path.join(RESULT_DIR, 'observed.csv')
    if not os.path.exists(obs_path):
        raise SystemExit(f"Missing {obs_path}. Did the observed run finish?")

    obs = load_run(obs_path)
    target_names = obs.index.tolist()
    observed = obs['cosine_similarity'].to_numpy(dtype=float)

    rows = []
    used_tags = []
    short = []
    
    for f in files:
        tag = os.path.basename(f)[:-4]
        rep = load_run(f)
        col = rep['cosine_similarity'].reindex(target_names)
        
        # If a target is missing from the replicate, track it and skip
        if col.isna().any():
            short.append((tag, int(col.isna().sum())))
            continue
            
        rows.append(col.to_numpy(dtype=float))
        used_tags.append(tag)

    if short:
        print(f"Dropped {len(short)} replicate files due to missing targets. First few: {short[:5]}")
        
    if not rows:
        raise SystemExit("No complete replicate files found. Bailing out.")

    null = np.vstack(rows)
    print(f"Real library: median cosine {np.median(observed):.3f}, max {observed.max():.3f}")
    print(f"Null matrix: {null.shape[0]} replicates by {null.shape[1]} targets\n")
    
    return observed, null, target_names, used_tags

def main():
    # Grab all the null library CSVs
    files = sorted(glob.glob(os.path.join(RESULT_DIR, 'null_*.csv')))
    if not files:
        raise SystemExit(f"No null_*.csv files found in {RESULT_DIR}")

    obs_df = load_run(os.path.join(RESULT_DIR, 'observed.csv'))
    print(f"Processing {len(obs_df)} COSMIC targets across {len(files)} shuffled libraries\n")

    observed, null, target_names, tags = build(files)
    n_reps = null.shape[0]

    # Sanity check: are the shuffles actually moving?
    spread = null.std(axis=0)
    frozen = [target_names[j] for j in np.where(spread < 1e-9)[0]]
    
    print("Checking variance across shuffles:")
    print(f"Targets identical across all {n_reps} replicates: {len(frozen)} / {len(target_names)}")
    if frozen:
        print("  " + ", ".join(frozen[:25]))
    print(f"Median spread of a target across replicates: {np.median(spread):.4f}\n")

    # Global null stats 
    pooled = null.ravel()
    print("Global background stats:")
    print(f"Total cosines from shuffled libraries: {pooled.size}")
    print(f"Median: {np.median(pooled):.3f}, 95th: {np.percentile(pooled, 95):.3f}, 99th: {np.percentile(pooled, 99):.3f}, Max: {pooled.max():.3f}")
    
    for t in (0.80, 0.85, 0.90, 0.95):
        print(f"Fraction reaching {t:.2f}: {np.mean(pooled >= t):.4f}")
    print(f"Real library median: {np.median(observed):.3f}, max: {observed.max():.3f}\n")

    # Per-target p-values 
    rows = []
    for j, name in enumerate(target_names):
        draws = null[:, j]
        o = observed[j]
        k = int((draws >= o).sum())
        
        rows.append({
            'cosmic_signature': name,
            'observed_cosine': round(float(o), 4),
            'null_n': n_reps,
            'null_median': round(float(np.median(draws)), 4),
            'null_p95': round(float(np.percentile(draws, 95)), 4),
            'null_max': round(float(draws.max()), 4),
            'margin_over_null_p95': round(float(o - np.percentile(draws, 95)), 4),
            'null_times_reached_observed': k,
            'p_value': (1 + k) / (n_reps + 1),
        })
        
    summary = pd.DataFrame(rows)
    summary['fdr'] = benjamini_hochberg(summary['p_value'].values)

    n_tested = len(summary)
    p_floor = 1.0 / (n_reps + 1)
    n_tied = int((summary['p_value'] <= p_floor + 1e-12).sum())
    
    print(" Resolution ceiling ")
    print(f"Smallest possible p-value: {p_floor:.5f} (based on {n_reps} replicates)")
    print(f"{n_tied} out of {n_tested} targets hit this floor.")
    if n_tied:
        print(f"Smallest reachable FDR given the ties: {p_floor * n_tested / n_tied:.5f}")
    print()

    # Final classification 
    summary['beats_chance_fdr5'] = summary['fdr'] <= 0.05
    for f in ADEQUACY_FLOORS:
        summary[f'cosine_at_least_{f}'] = summary['observed_cosine'] >= f
        summary[f'explained_at_{f}'] = summary['beats_chance_fdr5'] & summary[f'cosine_at_least_{f}']

    print(f" How many targets qualify? (out of {n_tested}) ")
    print(f"Cosine >= 0.90 alone (ignoring chance): {int((summary['observed_cosine'] >= 0.90).sum())}")
    print(f"Beats chance alone (FDR <= 0.05): {int(summary['beats_chance_fdr5'].sum())}")
    for f in ADEQUACY_FLOORS:
        print(f"Beats chance AND cosine >= {f:.2f}: {int(summary[f'explained_at_{f}'].sum())}")
    print()

    print(" Top 15 targets (ranked by margin over null 95th percentile) ")
    strong = summary.sort_values('margin_over_null_p95', ascending=False)
    print(strong[['cosmic_signature', 'observed_cosine', 'null_median', 'null_p95', 'null_max', 'margin_over_null_p95', 'p_value', 'fdr']].head(15).to_string(index=False))
    print()

    print(" Targets that fail to beat their own chance ceiling ")
    weak = summary[summary['margin_over_null_p95'] <= 0].sort_values('margin_over_null_p95')
    if len(weak):
        print(weak[['cosmic_signature', 'observed_cosine', 'null_median', 'null_p95', 'null_max', 'margin_over_null_p95', 'p_value']].to_string(index=False))
        print("Note: No etiological claim survives here, regardless of the raw cosine.")
    else:
        print("None!")
    print()

    lost = summary[(summary['observed_cosine'] >= 0.90) & ~summary[f'explained_at_{ADEQUACY_FLOORS[1]}']]
    if len(lost):
        print(f" Targets the old 0.90 rule called 'explained' that no longer qualify ({len(lost)}) ")
        print("  " + ", ".join(lost['cosmic_signature'].tolist()))
        print()

    out = os.path.join(BASE_DIR, 'null_SPA_summary.csv')
    summary.sort_values('margin_over_null_p95', ascending=False).to_csv(out, index=False)
    print(f"Saved results to {out}")

if __name__ == '__main__':
    main()