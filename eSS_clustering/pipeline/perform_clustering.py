#!/usr/bin/env python3
"""
Hierarchical clustering pipeline — supports SBS, DBS, and ID mutation types.

For SBS: loads separate raw-count and pre-normalised files.
For DBS/ID: loads raw counts only and normalises internally (column sums to 1).

Output is written to <output_dir>/<MUTATION_TYPE>/ so separate runs do not
overwrite each other.

FIXED VERSION: Exact dendrogram-based label matching + corrected data directory handling
"""

import argparse
import glob
import os
import re
import sys
import shutil

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import seaborn as sns
import yaml
from numpy import dot
from numpy.linalg import norm
from pdf2image import convert_from_path
from scipy.spatial.distance import pdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.naming import remove_replicate_suffix
from utils.mutation_type import get_config, call_plot_function


# ==============================================================================
# 0. SPECIES DETECTION
# ==============================================================================

SPECIES_PATTERNS = {
    'mouse':    (re.compile(r'mouse',          re.I), 'Mouse_'),
    'rat':      (re.compile(r'(?<![a-z])rat(?![a-z])', re.I), 'Rat_'),
    'chicken':  (re.compile(r'chicken',        re.I), 'Chicken_'),
    'human':    (re.compile(r'human',          re.I), 'Human_'),
    'celegans': (re.compile(r'c[._]?elegans?', re.I), 'celegans_'),
}


def detect_species(filename):
    """Return (species_key, prefix) or (None, None) if unrecognised."""
    for key, (pattern, prefix) in SPECIES_PATTERNS.items():
        if pattern.search(filename):
            return key, prefix
    return None, None


# ==============================================================================
# 1. SAMPLE MAPPING
# ==============================================================================

def load_sample_mapping(mapping_file_path):
    if not os.path.exists(mapping_file_path):
        print(f"Warning: Sample mapping file not found at {mapping_file_path}")
        return {}
    mapping_df = pd.read_csv(mapping_file_path, sep='\t')
    mapping_dict = dict(zip(mapping_df['sample_name'], mapping_df['standardized_name']))
    print(f"Loaded {len(mapping_dict)} sample mappings from {mapping_file_path}")
    return mapping_dict


def standardize_sample_name(original_name, species_prefix, mapping_dict):
    if species_prefix.lower().startswith('celegans'):
        return species_prefix + original_name
    name_for_lookup = remove_replicate_suffix(original_name)
    if name_for_lookup in mapping_dict:
        return mapping_dict[name_for_lookup]
    if original_name in mapping_dict:
        return mapping_dict[original_name]
    return species_prefix + remove_replicate_suffix(original_name)


def apply_sample_mapping(dataframe, species_prefix, mapping_dict):
    """Return (mapped_df, {original: standardized}) for one species."""
    new_columns = []
    original_to_standardized = {}
    unmapped = []

    for col in dataframe.columns:
        new_name = standardize_sample_name(col, species_prefix, mapping_dict)
        new_columns.append(new_name)
        original_to_standardized[col] = new_name
        if col not in mapping_dict and remove_replicate_suffix(col) not in mapping_dict:
            unmapped.append(col)

    if unmapped:
        print(f"{species_prefix} — {len(unmapped)} unmapped samples: "
              f"{unmapped[:5]}{'...' if len(unmapped) > 5 else ''}")
    else:
        print(f"{species_prefix} — all samples mapped.")

    # Handle duplicate standardised names
    if len(new_columns) != len(set(new_columns)):
        counts = {}
        deduped = []
        for i, name in enumerate(new_columns):
            orig = dataframe.columns[i]
            if name in counts:
                counts[name] += 1
                final = f"{name}_{counts[name]}"
            else:
                counts[name] = 0
                final = name
            deduped.append(final)
            original_to_standardized[orig] = final
        new_columns = deduped

    result = dataframe.copy()
    result.columns = new_columns
    return result, original_to_standardized


def create_reverse_mapping(all_mappings):
    reverse = {}
    for species, mapping in all_mappings.items():
        for original, standardized in mapping.items():
            reverse[standardized] = original
    return reverse


def convert_to_original_names(standardized_names, all_mappings):
    if all_mappings is None:
        return standardized_names
    reverse_map = {}
    for species, mapping in all_mappings.items():
        _, prefix = SPECIES_PATTERNS.get(species, (None, species + '_'))
        for original, standardized in mapping.items():
            reverse_map[standardized] = (prefix, original)

    result = []
    for name in standardized_names:
        if name in reverse_map:
            prefix, original = reverse_map[name]
            result.append(f"{prefix}{original}")
        else:
            result.append(name)
    return result


def save_sample_mappings(mappings_dict, output_file):
    rows = [
        {'species': sp, 'original_name': orig, 'standardized_name': std}
        for sp, mapping in mappings_dict.items()
        for orig, std in mapping.items()
    ]
    pd.DataFrame(rows).to_csv(output_file, sep='\t', index=False)
    print(f"Sample mappings saved to {output_file}")


# ==============================================================================
# 2. DATA LOADING
# ==============================================================================

def _load_species_files(data_dir, cfg, mapping_dict):
    """
    Discover all .txt/.tsv files in data_dir, detect species from filename,
    and load raw counts.

    Returns (counts_df, all_mappings).
    """
    found = sorted(
        f for ext in ('*.txt', '*.tsv')
        for f in glob.glob(os.path.join(data_dir, ext))
    )

    if not found:
        raise FileNotFoundError(f"No .txt/.tsv files found in {data_dir}")

    frames, all_mappings = [], {}

    for path in found:
        basename = os.path.basename(path)

        # Skip normalised files — those are handled separately for SBS
        if 'normaliz' in basename.lower():
            continue

        species, prefix = detect_species(basename)
        if species is None:
            print(f"  Skipping unrecognised file: {basename}")
            continue

        df = pd.read_csv(path, sep='\t', index_col=0)
        print(f"  {species} ({basename}): {df.shape[1]} samples loaded")

        if species == 'mouse' and cfg.mouse_exclude_pattern:
            df = df.loc[:, ~df.columns.str.contains(cfg.mouse_exclude_pattern)]
            print(f"  Mouse after exclusion filter: {df.shape[1]} samples")

        if species == 'celegans':
            original_cols = df.columns.tolist()
            df.columns = [prefix + c for c in original_cols]
            all_mappings[species] = {orig: prefix + orig for orig in original_cols}
            frames.append(df)
        else:
            mapped_df, mapping = apply_sample_mapping(df, prefix, mapping_dict)
            all_mappings[species] = mapping
            frames.append(mapped_df)

    if not frames:
        raise ValueError(f"No recognisable species files found in {data_dir}")

    counts_df = pd.concat(frames, axis=1)
    print(f"Combined counts shape: {counts_df.shape}")
    return counts_df, all_mappings


def load_data(data_dir, cfg, mapping_file_path):
    """
    Load and return (counts_df, normalized_df, all_mappings).

    For SBS, normalized_df is loaded from pre-computed files (detected by
    'normaliz' in filename + species detection).
    For DBS/ID, normalized_df is derived from counts_df (column-normalised).
    """
    print(f"\nLoading {cfg.name} data from: {data_dir}")
    mapping_dict = load_sample_mapping(mapping_file_path)

    # --- raw counts (all mutation types) ---
    counts_df, all_mappings = _load_species_files(data_dir, cfg, mapping_dict)

    # --- normalised profiles ---
    if cfg.has_prenormalized_files:
        print("Loading pre-normalised files (SBS)...")
        found = sorted(
            f for ext in ('*.txt', '*.tsv')
            for f in glob.glob(os.path.join(data_dir, ext))
        )

        norm_frames = []
        for path in found:
            basename = os.path.basename(path)

            if 'normaliz' not in basename.lower():
                continue

            species, prefix = detect_species(basename)
            if species is None:
                print(f"  Skipping unrecognised normalised file: {basename}")
                continue

            df = pd.read_csv(path, sep='\t', index_col=0)
            print(f"  {species} normalised ({basename}): {df.shape[1]} samples loaded")

            if species == 'mouse' and cfg.mouse_exclude_pattern:
                df = df.loc[:, ~df.columns.str.contains(cfg.mouse_exclude_pattern)]

            if species == 'celegans':
                df.columns = [prefix + c for c in df.columns]
            else:
                df, _ = apply_sample_mapping(df, prefix, mapping_dict)

            norm_frames.append(df)

        if not norm_frames:
            raise ValueError(
                f"No normalised files (containing 'normaliz' in filename) "
                f"found in {data_dir} for SBS."
            )

        normalized_df = pd.concat(norm_frames, axis=1)
        print(f"Combined normalised shape: {normalized_df.shape}")

    else:
        print(f"Normalising {cfg.name} raw counts internally (column sums → 1)...")
        col_sums = counts_df.sum(axis=0)
        zero_cols = col_sums[col_sums == 0].index.tolist()
        if zero_cols:
            print(f"  Warning: {len(zero_cols)} all-zero columns — removing.")
            counts_df = counts_df.drop(columns=zero_cols)
            col_sums  = col_sums.drop(index=zero_cols)
        normalized_df = counts_df.div(col_sums, axis=1)
        print("Normalisation complete.")

    return counts_df, normalized_df, all_mappings


def filter_zero_columns(df, output_dir):
    all_zero_mask = (df == 0).all(axis=0)
    nonzero_df  = df.loc[:, ~all_zero_mask]
    allzero_df  = df.loc[:,  all_zero_mask]
    path = os.path.join(output_dir, "all_zero_samples_removed.txt")
    allzero_df.to_csv(path, sep='\t')
    print(f"Kept {nonzero_df.shape[1]} non-zero samples; "
          f"removed {allzero_df.shape[1]} (saved to {path})")
    return nonzero_df


# ==============================================================================
# 3. CLUSTERING
# ==============================================================================

def cluster_signatures_with_custom_thresholds(X, output_dir, output_name,
                                              num_colors=100,
                                              cos_dist_threshold=0.15,
                                              custom_thresholds=None,
                                              reverse_mapping=None):
    """
    Hierarchical clustering with optional custom per-cluster thresholds.

    Returns
    -------
    (main_clusters, small_clusters, filtered_ordering, all_ordering,
     filtered_colors, all_colors, true_singletons)
    """
    print(f"\nClustering with main threshold: {cos_dist_threshold:.4f}")
    if custom_thresholds:
        print(f"Custom thresholds: {custom_thresholds}")

    # Define special samples that should always be colored (even if cluster size ≤ 2)
    special_patterns = [
        "MEF_AID",
        "MCF10_cisplatin",
        # Add more patterns here as needed
    ]

    np.random.seed(42)
    D = pdist(X.T, 'cosine')
    Z = sch.linkage(D, 'average')

    palette   = sns.color_palette("husl", num_colors)
    hex_colors = [mcolors.rgb2hex(c) for c in palette]
    sch.set_link_color_palette(hex_colors)

    def llf(idx):
        std_name = X.columns[idx]
        if not reverse_mapping or std_name.startswith('celegans_'):
            return std_name
        if std_name in reverse_mapping:
            orig = reverse_mapping[std_name]
            for prefix in ['Mouse_', 'Rat_', 'Chicken_', 'Human_']:
                if std_name.startswith(prefix):
                    return prefix + orig
        return std_name

    # ── Initial clustering ───────────────────────────────────────────────────
    fig_tmp, ax_tmp = plt.subplots(figsize=(10, 6))
    R_orig = sch.dendrogram(Z, color_threshold=cos_dist_threshold,
                            above_threshold_color='gray',
                            leaf_label_func=llf, ax=ax_tmp)
    plt.close(fig_tmp)

    initial_clusters = {}
    initial_colors   = {}
    for color, leaf_idx in zip(R_orig["leaves_color_list"], R_orig['leaves']):
        name = X.columns[leaf_idx]
        initial_clusters.setdefault(color, []).append(name)
        initial_colors[name] = color

    true_singletons = initial_clusters.pop('gray', [])
    print(f"True singletons (gray): {len(true_singletons)}")

    # ── Apply custom thresholds ──────────────────────────────────────────────
    final_clusters = {}
    final_colors   = {}
    available      = list(hex_colors)
    used           = set(initial_colors.values())
    unused         = [c for c in available if c not in used and c != 'gray']
    next_idx       = 0

    if custom_thresholds:
        for color, samples in initial_clusters.items():
            ct = next(
                (t for pat, t in custom_thresholds.items()
                 if any(pat in s for s in samples)),
                None
            )
            if ct is not None and ct < cos_dist_threshold and len(samples) > 1:
                print(f"  Applying custom threshold {ct} to cluster {color}")
                idxs  = [X.columns.get_loc(s) for s in samples]
                sub_D = pdist(X.iloc[:, idxs].T, 'cosine')
                sub_Z = sch.linkage(sub_D, 'average')
                labels = sch.fcluster(sub_Z, ct, criterion='distance')
                subclusters = {}
                for i, lbl in enumerate(labels):
                    subclusters.setdefault(lbl, []).append(samples[i])
                for k, (_, sub_samps) in enumerate(
                        sorted(subclusters.items(),
                               key=lambda x: len(x[1]), reverse=True)):
                    sub_color = color if k == 0 else (
                        unused[next_idx] if next_idx < len(unused)
                        else available[next_idx % len(available)]
                    )
                    if k > 0:
                        next_idx += 1
                    final_clusters[sub_color] = sub_samps
                    for s in sub_samps:
                        final_colors[s] = sub_color
            else:
                final_clusters[color] = samples
                for s in samples:
                    final_colors[s] = color
    else:
        final_clusters = dict(initial_clusters)
        final_colors   = dict(initial_colors)

    # ── Final dendrogram plot ────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(60, 6), dpi=400)
    R_final = sch.dendrogram(Z, color_threshold=cos_dist_threshold,
                             above_threshold_color='gray',
                             leaf_label_func=llf, ax=ax)

    # ── Manual label recoloring (fixed version using exact mapping) ─────────
    # Build exact mapping from dendrogram label text to standardized sample names
    label_to_std = {}
    for leaf_idx in R_final['leaves']:
        std_name = X.columns[leaf_idx]  # True standardized name from data
        label_text = llf(leaf_idx)       # What's displayed in dendrogram
        label_to_std[label_text] = std_name

    # Color labels based on cluster membership
    for label in ax.get_xmajorticklabels():
        label_text = label.get_text()
        
        # Exact lookup - no fuzzy matching
        std_name = label_to_std.get(label_text)
        
        if std_name and std_name in final_colors:
            # Find which cluster this sample belongs to
            cluster_members = next(
                (v for v in final_clusters.values() if std_name in v), 
                None
            )
            
            if cluster_members:
                # Check if cluster should be colored
                is_large = len(cluster_members) > 2
                
                # Check if cluster contains any special samples
                is_special = any(
                    any(pattern in sample for pattern in special_patterns)
                    for sample in cluster_members
                )
                
                # Color if large OR special
                if is_large or is_special:
                    label.set_color(final_colors[std_name])
        
        # Set consistent font size
        label.set_fontsize(4)

    ax.axhline(y=cos_dist_threshold, color='black',  linestyle='-',  linewidth=2,
               label=f'Cosine similarity: {1 - cos_dist_threshold:.2f}')
    ax.axhline(y=0.15,               color='dimgray', linestyle='--', linewidth=2,
               label='Cosine similarity: 0.85')
    ax.legend(fontsize=12)
    ax.set_title('Agglomerative Hierarchical Clustering of eSignature Models',
                 fontsize=18)
    ax.set_xlabel('Signatures', fontsize=16)
    ax.set_ylabel('Cosine Distance',  fontsize=16)
    ax.tick_params(axis='y', labelsize=14)
    plt.tight_layout()

    for ext in ['pdf', 'png', 'svg']:
        path = os.path.join(output_dir, f"{output_name}.{ext}")
        plt.savefig(path, dpi=400, format=ext)
        print(f"Dendrogram saved: {path}")
    plt.close()

    # ── Filter into main / small ─────────────────────────────────────────────
    main_clusters  = {}
    small_clusters = {}
    for color, samples in final_clusters.items():
        is_large   = len(samples) > 2
        is_special = any(any(sp in s for sp in special_patterns) for s in samples)
        if is_large or is_special:
            main_clusters[color] = samples
        else:
            small_clusters[color] = samples

    filtered_colors = {s: c for s, c in final_colors.items()
                       if c in main_clusters}

    filtered_ordering = [s for s in R_final['ivl'] if s in filtered_colors]
    all_ordering      = R_final['ivl']

    print(f"Main clusters: {len(main_clusters)}  |  "
          f"Small clusters: {len(small_clusters)}  |  "
          f"Singletons: {len(true_singletons)}")

    return (main_clusters, small_clusters, filtered_ordering, all_ordering,
            filtered_colors, final_colors, true_singletons)


# ==============================================================================
# 4. CLUSTER SUMMARISATION
# ==============================================================================

def summarize_and_average_clusters(clusters, normalized_df, counts_df,
                                   output_dir, cluster_type, all_mappings,
                                   cluster_prefix, averaging_method="equal_replicate"):
    """
    Save cluster summary TSV and compute per-cluster consensus profiles.

    cluster_prefix : str  e.g. "eSS" / "eDS" / "eIS"
    averaging_method : str
        "equal_replicate" (default) Each sample is normalized to sum to 1
                           first, then the consensus profile is the
                           unweighted mean of those per-sample normalized
                           profiles across the cluster. Every sample/
                           replicate contributes equally regardless of its
                           mutation burden -- avoids one high-burden
                           replicate dominating the consensus.
        "pooled"           Sum raw counts per channel across a cluster's
                           samples, then normalize by the cluster's total
                           mutation count. Samples with a higher mutation
                           burden contribute proportionally more to the
                           consensus profile.

    Note: this only affects the reported consensus profile for each
    cluster (used for plotting / COSMIC comparison) -- it has no effect on
    cluster membership, which is decided upstream during clustering.
    """
    if averaging_method not in ("pooled", "equal_replicate"):
        raise ValueError(
            f"Unknown averaging_method {averaging_method!r}; "
            f"must be 'pooled' or 'equal_replicate'."
        )
    print(f"\n{'='*60}")
    print(f"Summarising {cluster_type.upper()} clusters")
    print(f"{'='*60}")

    # Summary TSV (sample → cluster assignment)
    rows = [(sample, cid + 1, color)
            for cid, (color, samples) in enumerate(clusters.items())
            for sample in samples]
    summary_df = pd.DataFrame(rows, columns=["Sample", "Cluster", "Color"])
    summary_path = os.path.join(output_dir, f"{cluster_type}_clusters_summary.tsv")
    summary_df.to_csv(summary_path, index=False, sep='\t')
    print(f"Cluster summary → {summary_path}")

    # Weighted average profiles
    columns     = {}
    mut_stats   = []

    for num, (color, sigs) in enumerate(clusters.items()):
        cluster_id = f"{cluster_prefix}{num + 1}"

        counts_sub = counts_df.filter(sigs)
        norm_sub   = normalized_df.filter(sigs)

        # Per-sample mutation counts
        per_sample = counts_sub.sum(axis=0)
        cluster_total = int(counts_sub.values.sum())

        print(f"\n  {cluster_id}: {len(sigs)} samples, "
              f"{cluster_total:,} total mutations")
        if all_mappings:
            display = convert_to_original_names(sigs, all_mappings)
        else:
            display = sigs
        for dname, count in zip(display, per_sample):
            print(f"    {dname:50s} {int(count):>10,}")

        if averaging_method == "equal_replicate":
            # Each sample already sums to 1 in normalized_df; take the
            # unweighted mean across samples so no single (high-burden)
            # replicate dominates the consensus profile.
            columns[cluster_id] = norm_sub.mean(axis=1).tolist()
        else:
            # Pooled-counts consensus profile: sum raw counts across samples
            # per channel, then normalize by the cluster's total mutation
            # count. Equivalent to sum_s(counts[c,s]) / sum_s(sample_total[s])
            # for every c -- higher-burden samples dominate proportionally.
            pooled_counts = counts_sub.sum(axis=1)          # per-channel sum across samples
            columns[cluster_id] = (pooled_counts / cluster_total).tolist()

        mut_stats.append({
            'Cluster_ID': cluster_id,
            'N_Samples':  len(sigs),
            'Total_Mutations': cluster_total,
            'Samples': '|'.join(
                convert_to_original_names(sigs, all_mappings)
                if all_mappings else sigs
            ),
        })

    # Mutation stats TSV
    stats_path = os.path.join(output_dir,
                              f"{cluster_type}_clusters_mutation_counts.tsv")
    pd.DataFrame(mut_stats).to_csv(stats_path, sep='\t', index=False)
    print(f"\nMutation stats → {stats_path}")

    # Profiles TSV
    final_df   = pd.DataFrame(columns, index=normalized_df.index)
    col_sums   = final_df.sum(axis=0)
    if not col_sums.between(0.99, 1.01).all():
        final_df = final_df.div(col_sums, axis=1)

    matrix_path = os.path.join(
        output_dir,
        f"{cluster_type}_clusters_NORMALIZED_weighted_avg_profiles.tsv"
    )
    final_df.to_csv(matrix_path, sep='\t')
    print(f"Profiles matrix → {matrix_path}")
    return matrix_path


# ==============================================================================
# 5. SINGLETONS
# ==============================================================================

def save_true_singleton_samples(true_singletons, normalized_df, output_dir):
    if not true_singletons:
        print("No singleton samples.")
        return None
    print(f"Saving {len(true_singletons)} singletons...")
    singleton_data = normalized_df[true_singletons]
    path = os.path.join(output_dir, "singleton_small_clusters_ordered.tsv")
    singleton_data.to_csv(path, sep='\t')
    names_path = os.path.join(output_dir,
                              "singleton_small_clusters_sample_names_ordered.txt")
    with open(names_path, 'w') as f:
        f.writelines(f"{s}\n" for s in true_singletons)
    print(f"Singleton data → {path}")
    return path


# ==============================================================================
# 6. PLOTTING AND PDF CONVERSION
# ==============================================================================

def plot_signatures(matrix_path, output_dir, project_name, cfg):
    """Call the appropriate sigProfilerPlotting function for this mutation type."""
    plot_subdir = os.path.join(output_dir, f"{project_name}_plots")
    os.makedirs(plot_subdir, exist_ok=True)
    print(f"Plotting {cfg.name} signatures for {project_name}...")
    call_plot_function(cfg, matrix_path, plot_subdir, project_name)
    print(f"Plots saved in {plot_subdir}")
    return plot_subdir


def convert_pdf_to_pngs(pdf_path, output_dir, prefix="cluster"):
    if not os.path.exists(pdf_path):
        print(f"Warning: PDF not found — {pdf_path}")
        return
    print(f"Converting {os.path.basename(pdf_path)} to PNGs...")
    images = convert_from_path(pdf_path, dpi=300)
    for i, img in enumerate(images):
        img.save(os.path.join(output_dir, f'{prefix}-{i + 1:02d}.png'), 'PNG')
    print(f"Saved {len(images)} PNGs to {output_dir}")


# ==============================================================================
# 6B. INTEGRATED PREPROCESSING
# ==============================================================================

def load_preprocessing_config(config_path):
    """Load preprocessing configuration from YAML file."""
    if not os.path.exists(config_path):
        return None
    
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        print(f"Warning: Failed to load preprocessing config: {e}")
        return None


def should_exclude_sample(sample_name, exclusion_patterns, case_sensitive=False):
    """Check if sample name matches any exclusion pattern."""
    if not case_sensitive:
        sample_name = sample_name.lower()
        exclusion_patterns = [p.lower() for p in exclusion_patterns]
    
    for pattern in exclusion_patterns:
        if pattern in sample_name:
            return True, pattern
    return False, None


def preprocess_file(input_path, output_path, exclusion_patterns, case_sensitive=False):
    """Remove samples matching exclusion patterns from a single file."""
    df = pd.read_csv(input_path, sep='\t', index_col=0)
    original_count = df.shape[1]
    
    # Identify samples to remove
    removed_samples = []
    kept_columns = []
    
    for col in df.columns:
        exclude, matched_pattern = should_exclude_sample(col, exclusion_patterns, case_sensitive)
        if exclude:
            removed_samples.append((col, matched_pattern))
        else:
            kept_columns.append(col)
    
    # Filter and save
    df_filtered = df[kept_columns]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_filtered.to_csv(output_path, sep='\t')
    
    return {
        'original_count': original_count,
        'removed_count': len(removed_samples),
        'kept_count': len(kept_columns),
        'removed_samples': removed_samples,
    }


def run_preprocessing(input_dir, output_dir, exclusion_patterns, case_sensitive=False):
    """Preprocess all files in input directory."""
    print("\n" + "="*70)
    print("PREPROCESSING")
    print("="*70)
    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}")
    print(f"Exclusion patterns: {exclusion_patterns}")
    print("="*70)
    
    # Find all data files
    files = []
    for ext in ['*.txt', '*.tsv']:
        files.extend(glob.glob(os.path.join(input_dir, ext)))
    
    if not files:
        print(f"\n⚠️  No .txt or .tsv files found in {input_dir}")
        return {}
    
    total_removed = 0
    total_original = 0
    
    # Process each file
    for input_path in sorted(files):
        filename = os.path.basename(input_path)
        output_path = os.path.join(output_dir, filename)
        
        print(f"\n  Processing: {filename}")
        stats = preprocess_file(input_path, output_path, exclusion_patterns, case_sensitive)
        
        total_original += stats['original_count']
        total_removed += stats['removed_count']
        
        print(f"    Original: {stats['original_count']} samples")
        print(f"    Removed:  {stats['removed_count']} samples")
        print(f"    Kept:     {stats['kept_count']} samples")
        
        if stats['removed_samples']:
            print(f"    Matched patterns (first 3):")
            for sample, pattern in stats['removed_samples'][:3]:
                print(f"      - {sample} ('{pattern}')")
            if len(stats['removed_samples']) > 3:
                print(f"      ... and {len(stats['removed_samples']) - 3} more")
    
    print(f"\n{'='*70}")
    print(f"PREPROCESSING COMPLETE")
    print(f"  Total samples: {total_original}")
    print(f"  Removed: {total_removed}")
    print(f"  Kept: {total_original - total_removed}")
    print(f"  Cleaned data saved to: {output_dir}")
    print('='*70 + "\n")
    
    return True


def determine_data_directory(base_input_dir, mutation_type, preprocessing_config_path, 
                             force_preprocess=False, skip_preprocessing=False):
    """
    Determine which data directory to use based on preprocessing config.
    
    Returns
    -------
    str : Path to data directory to use (either original or cleaned)
    """
    original_dir = os.path.join(base_input_dir, mutation_type)
    cleaned_dir = os.path.join(base_input_dir + '_cleaned', mutation_type)
    
    print("\n" + "="*70)
    print("DATA DIRECTORY SELECTION")
    print("="*70)
    
    # Check if user wants to skip preprocessing
    if skip_preprocessing:
        print("✓ Preprocessing skipped (--skip_preprocessing)")
        print(f"✓ Using original data: {original_dir}")
        return original_dir
    
    # Load preprocessing config
    config = load_preprocessing_config(preprocessing_config_path)
    
    if config is None:
        print(f"✗ No preprocessing config found: {preprocessing_config_path}")
        print(f"✓ Using original data: {original_dir}")
        return original_dir
    
    print(f"✓ Found preprocessing config: {preprocessing_config_path}")
    
    # Check if there are exclusion patterns for this mutation type
    exclusion_patterns = config.get(mutation_type, {}).get('exclude', [])
    
    if not exclusion_patterns:
        print(f"✗ No exclusion patterns defined for {mutation_type}")
        print(f"✓ Using original data: {original_dir}")
        return original_dir
    
    print(f"✓ Exclusion patterns for {mutation_type}: {exclusion_patterns}")
    
    # Check if cleaned directory already exists (and we're not forcing reprocessing)
    cleaned_exists = os.path.exists(cleaned_dir)
    has_cleaned_files = False
    
    if cleaned_exists:
        has_cleaned_files = any(
            glob.glob(os.path.join(cleaned_dir, ext))
            for ext in ['*.txt', '*.tsv']
        )
    
    if cleaned_exists and has_cleaned_files and not force_preprocess:
        print(f"✓ Using cached cleaned data: {cleaned_dir}")
        print("  (Use --force_preprocess to regenerate)")
        return cleaned_dir
    
    # Need to run preprocessing
    if force_preprocess and cleaned_exists:
        print(f"⚠️  Forcing reprocessing (cleaning existing directory)")
        shutil.rmtree(cleaned_dir)
    
    print(f"\nRunning preprocessing...")
    case_sensitive = config.get('case_sensitive', False)
    
    success = run_preprocessing(
        original_dir,
        cleaned_dir,
        exclusion_patterns,
        case_sensitive
    )
    
    if success:
        print(f"✓ Using cleaned data: {cleaned_dir}")
        return cleaned_dir
    else:
        print(f"⚠️  Preprocessing failed, falling back to original data")
        return original_dir


def validate_data_dir_argument(data_dir_arg, mutation_type):
    """
    Detect if user passed mutation-type-specific directory instead of base.
    
    Returns the corrected base directory.
    """
    # Check if the path ends with the mutation type
    if data_dir_arg.rstrip('/').endswith(f"/{mutation_type}") or \
       data_dir_arg.rstrip('/').endswith(f"\\{mutation_type}"):
        # User passed the full path, extract the base
        corrected = data_dir_arg.rstrip('/').rstrip('\\')
        if corrected.endswith(mutation_type):
            corrected = os.path.dirname(corrected)
        
        print(f"\n⚠️  WARNING: --data_dir should be the BASE directory")
        print(f"  You passed: {data_dir_arg}")
        print(f"  Auto-correcting to: {corrected}")
        print(f"  (The /{mutation_type} directory will be appended automatically)\n")
        return corrected
    
    return data_dir_arg


# ==============================================================================
# 7. MAIN
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="eSignatures hierarchical clustering pipeline with integrated preprocessing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--output_dir",       required=True,
                        help="Root output directory. Results go to "
                             "<output_dir>/<MUTATION_TYPE>/")
    parser.add_argument("--mutation_type",    default="SBS",
                        choices=["SBS", "DBS", "ID"],
                        help="Mutation type to process.")
    parser.add_argument("--data_dir",         default="data/input",
                        help="Base input directory (e.g., data/input). "
                             "Specific mutation type directory will be appended automatically.")
    parser.add_argument("--cosine_similarity", type=float, default=0.9,
                        help="Main cosine similarity threshold for clustering.")
    parser.add_argument("--custom_thresholds", default="",
                        help="Custom thresholds: 'pattern:value,pattern:value'.")
    parser.add_argument("--mapping_file",     default="sample_mapping.tsv",
                        help="Path to sample name mapping TSV.")
    parser.add_argument("--use_original_names", action="store_true",
                        help="Show original (pre-mapping) names in dendrogram.")
    parser.add_argument("--averaging_method", default="equal_replicate",
                        choices=["equal_replicate", "pooled"],
                        help="Consensus-profile method for main/small clusters. "
                             "'equal_replicate' (default): normalize each sample "
                             "first, then take the unweighted mean across samples "
                             "(every replicate contributes equally, regardless of "
                             "mutation burden). 'pooled': sum raw counts across "
                             "samples, then normalize once (higher-burden samples "
                             "dominate). Does not affect cluster membership.")

    # Preprocessing options
    parser.add_argument("--preprocessing_config", default="config/preprocessing.yaml",
                        help="Path to preprocessing configuration file.")
    parser.add_argument("--force_preprocess", action="store_true",
                        help="Force reprocessing even if cleaned data exists.")
    parser.add_argument("--skip_preprocessing", action="store_true",
                        help="Skip preprocessing even if config exists.")
    
    args = parser.parse_args()

    cfg = get_config(args.mutation_type)

    # Validate and auto-correct data_dir if needed
    args.data_dir = validate_data_dir_argument(args.data_dir, cfg.name)

    # Parse custom thresholds
    custom_thresholds = {}
    if args.custom_thresholds:
        for item in args.custom_thresholds.split(','):
            if ':' in item:
                pat, val = item.split(':', 1)
                custom_thresholds[pat.strip()] = float(val.strip())

    # ── Determine data directory (with integrated preprocessing) ─────────────
    data_dir = determine_data_directory(
        base_input_dir=args.data_dir,
        mutation_type=cfg.name,
        preprocessing_config_path=args.preprocessing_config,
        force_preprocess=args.force_preprocess,
        skip_preprocessing=args.skip_preprocessing
    )

    # ── Output directory namespaced by mutation type ─────────────────────────
    typed_output_dir = os.path.join(args.output_dir, cfg.name)
    interactive_dir  = os.path.join(typed_output_dir, "interactive_clustering")
    os.makedirs(typed_output_dir,  exist_ok=True)
    os.makedirs(interactive_dir,   exist_ok=True)

    print("=" * 60)
    print(f"eSignatures Clustering — {cfg.name}")
    print(f"Output: {typed_output_dir}")
    print(f"Cosine similarity ≥ {args.cosine_similarity}  "
          f"(distance < {1 - args.cosine_similarity:.4f})")
    print(f"Custom thresholds: {custom_thresholds}")
    print(f"Averaging method: {args.averaging_method}")
    print("=" * 60)

    # ── Load data ─────────────────────────────────────────────────────────────
    counts_df, normalized_df, all_mappings = load_data(
        data_dir, cfg, args.mapping_file
    )

    counts_df     = filter_zero_columns(counts_df, typed_output_dir)
    normalized_df = normalized_df[counts_df.columns]

    # Save applied mappings for downstream scripts
    save_sample_mappings(
        all_mappings,
        os.path.join(typed_output_dir, "applied_sample_name_mappings.tsv")
    )

    # ── Clustering ───────────────────────────────────────────────────────────
    reverse_mapping = (create_reverse_mapping(all_mappings)
                       if args.use_original_names else None)
    cos_dist = 1 - args.cosine_similarity

    (main_clusters, small_clusters, filtered_ordering, all_ordering,
     filtered_colors, all_colors, true_singletons) = (
        cluster_signatures_with_custom_thresholds(
            X=normalized_df,
            output_dir=interactive_dir,
            output_name="all_esignature_models_dendrogram",
            cos_dist_threshold=cos_dist,
            custom_thresholds=custom_thresholds,
            reverse_mapping=reverse_mapping,
        )
    )

    # ── Sub-directories ───────────────────────────────────────────────────────
    main_dir     = os.path.join(typed_output_dir, "main_clusters")
    small_dir    = os.path.join(typed_output_dir, "small_clusters")
    singletons_dir = os.path.join(typed_output_dir, "singletons")
    for d in [main_dir, small_dir, singletons_dir]:
        os.makedirs(d, exist_ok=True)

    # ── Main clusters ─────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("MAIN CLUSTERS")
    print("=" * 60)
    main_matrix = summarize_and_average_clusters(
        main_clusters, normalized_df, counts_df, main_dir,
        "main", all_mappings, cfg.cluster_prefix,
        averaging_method=args.averaging_method
    )
    main_project  = f"main_clusters_hierarchical"
    main_plot_dir = plot_signatures(main_matrix, main_dir, main_project, cfg)
    main_pdf      = os.path.join(main_plot_dir,
                                 f"{cfg.pdf_prefix}_{main_project}.pdf")
    convert_pdf_to_pngs(main_pdf, main_plot_dir, "main_cluster")

    # ── Small clusters ────────────────────────────────────────────────────────
    if small_clusters:
        print("\n" + "=" * 60)
        print("SMALL CLUSTERS")
        print("=" * 60)
        small_matrix = summarize_and_average_clusters(
            small_clusters, normalized_df, counts_df, small_dir,
            "small", all_mappings, cfg.cluster_prefix,
            averaging_method=args.averaging_method
        )
        small_project  = "small_clusters_hierarchical"
        small_plot_dir = plot_signatures(small_matrix, small_dir,
                                         small_project, cfg)
        small_pdf      = os.path.join(small_plot_dir,
                                      f"{cfg.pdf_prefix}_{small_project}.pdf")
        convert_pdf_to_pngs(small_pdf, small_plot_dir, "small_cluster")
    else:
        print("No small clusters.")

    # ── Singletons ───────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("SINGLETONS")
    print("=" * 60)
    singleton_matrix = save_true_singleton_samples(
        true_singletons, normalized_df, singletons_dir
    )
    if singleton_matrix:
        singleton_project  = "singleton_samples"
        singleton_plot_dir = plot_signatures(singleton_matrix, singletons_dir,
                                             singleton_project, cfg)
        singleton_pdf      = os.path.join(
            singleton_plot_dir,
            f"{cfg.pdf_prefix}_{singleton_project}.pdf"
        )
        convert_pdf_to_pngs(singleton_pdf, singleton_plot_dir, "singleton")

    # ── Done ─────────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print(f"COMPLETE — {cfg.name}")
    print(f"Output: {typed_output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()