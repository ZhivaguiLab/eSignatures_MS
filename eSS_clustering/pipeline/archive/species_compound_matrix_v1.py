#!/usr/bin/env python3
"""
Species-Model vs Compound Matrix Heatmap Generator

Creates a matrix heatmap where:
- Columns = Compound exposures (replicate suffixes removed)
- Rows    = Species-model combinations (grouped by species)
- Cells   = Pie charts showing clustering outcomes
             (main clusters, small clusters, singletons)

Outputs three PDF files when --other_species_plot is used:
  <output_path>_human.pdf
  <output_path>_mouse.pdf
  <output_path>_other_species.pdf
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from collections import defaultdict

# ---------------------------------------------------------------------------
# Allow running from the pipeline root without installing as a package
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.naming import (
    load_acronym_mapping,
    parse_sample_into_species_model_compound,
)

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
CLUSTER_OUTCOME_COLORS = {
    'main_cluster':  '#61BEEA',   # blue
    'small_cluster': '#F5CB62',   # amber
    'singleton':     '#F97070',   # red
    'not_in_cluster':'#D6D6D6',   # grey (reserved)
}


# ===========================================================================
# DATA LOADING
# ===========================================================================

def _load_summary_samples(cluster_type_dir, cluster_type):
    """Return sample list from <type>_clusters_summary.tsv, or []."""
    path = os.path.join(cluster_type_dir, f"{cluster_type}_clusters_summary.tsv")
    if not os.path.exists(path):
        print(f"No {cluster_type} cluster summary at: {path}")
        return []
    try:
        df = pd.read_csv(path, sep='\t')
        if 'Sample' not in df.columns:
            print(f"Warning: 'Sample' column missing in {path}")
            return []
        samples = df['Sample'].tolist()
        print(f"  Loaded {len(samples)} {cluster_type}-cluster samples")
        return samples
    except Exception as e:
        print(f"Error loading {cluster_type} cluster samples: {e}")
        return []


def _load_singleton_samples(singletons_dir):
    """Return sample list from singleton_small_clusters_ordered.tsv, or []."""
    path = os.path.join(singletons_dir, "singleton_small_clusters_ordered.tsv")
    if not os.path.exists(path):
        print(f"No singleton file at: {path}")
        return []
    try:
        df = pd.read_csv(path, sep='\t', index_col=0)
        samples = df.columns.tolist()
        print(f"  Loaded {len(samples)} singleton samples")
        return samples
    except Exception as e:
        print(f"Error loading singleton samples: {e}")
        return []


def collect_species_compound_data(main_output_dir, acronym_mapping=None):
    """
    Walk the three cluster output directories and build a summary dict.

    Returns
    -------
    dict
        {(species_model, compound): {
            'total': int, 'clustered': int,
            'small_clusters': int, 'singletons': int
        }}
    """
    print("Collecting species-model vs compound clustering data...")

    data = defaultdict(lambda: {
        'total': 0, 'clustered': 0, 'small_clusters': 0, 'singletons': 0
    })

    # Prefer capitalised compound names when the same compound appears under
    # different casings across files.
    _case_pref = {}   # lowercase compound → preferred display form

    def _normalize(compound):
        key = compound.lower()
        if key not in _case_pref:
            _case_pref[key] = compound
        elif compound != compound.lower():   # has capitals → prefer it
            _case_pref[key] = compound
        return _case_pref[key]

    def _process(sample, source):
        result = parse_sample_into_species_model_compound(
            sample, acronym_mapping, debug=False
        )
        if result is None:
            return
        species_model, compound = result
        compound = _normalize(compound)
        key = (species_model, compound)
        data[key]['total'] += 1
        data[key][source]  += 1

    # Main clusters
    main_dir = os.path.join(main_output_dir, "main_clusters")
    if os.path.exists(main_dir):
        summary = os.path.join(main_dir, "main_clusters_summary.tsv")
        if os.path.exists(summary):
            df = pd.read_csv(summary, sep='\t')
            print(f"Processing {len(df)} main-cluster entries...")
            for sample in df['Sample']:
                _process(sample, 'clustered')

    # Small clusters
    small_dir = os.path.join(main_output_dir, "small_clusters")
    if os.path.exists(small_dir):
        for sample in _load_summary_samples(small_dir, "small"):
            _process(sample, 'small_clusters')

    # Singletons
    singletons_dir = os.path.join(main_output_dir, "singletons")
    if os.path.exists(singletons_dir):
        for sample in _load_singleton_samples(singletons_dir):
            _process(sample, 'singletons')

    result = dict(data)
    print(f"Found {len(result)} unique (species_model, compound) combinations")
    return result


# ===========================================================================
# DATA ORGANISATION
# ===========================================================================

def organize_matrix_data(species_compound_data):
    """
    Sort species-models and compounds ready for plotting.

    Species order: Human → Mouse → Rat → Chicken → C. elegans → other.
    Within a species, models are sorted alphabetically.
    Compounds are sorted alphabetically.

    Returns
    -------
    tuple
        (sorted_species_models, sorted_compounds, species_compound_data)
    """
    all_species_models = set()
    all_compounds      = set()

    for (sm, compound) in species_compound_data:
        all_species_models.add(sm)
        all_compounds.add(compound)

    species_order = {'human': 1, 'mouse': 2, 'rat': 3, 'chicken': 4}

    def _sort_key(sm):
        if sm == "C.elegans":
            return (5, "elegans")
        parts   = sm.split('_')
        species = parts[0].lower()
        model   = '_'.join(parts[1:]) if len(parts) > 1 else ""
        return (species_order.get(species, 999), model.lower())

    sorted_species_models = sorted(all_species_models, key=_sort_key)
    sorted_compounds      = sorted(all_compounds)

    print(f"Matrix: {len(sorted_species_models)} species-models × "
          f"{len(sorted_compounds)} compounds")
    return sorted_species_models, sorted_compounds, species_compound_data


# ===========================================================================
# VISUALISATION HELPERS
# ===========================================================================

def _pie_patch(ax, cx, cy, radius, proportions, colors, counts):
    """Draw a pie chart at (cx, cy) using Wedge patches."""
    if not proportions or sum(proportions) == 0:
        return

    single_color = sum(1 for p in proportions if p > 0) == 1
    angle = 90.0   # start at 12 o'clock

    for prop, color, count in zip(proportions, colors, counts):
        if prop <= 0:
            continue
        sweep = prop * 360
        ax.add_patch(patches.Wedge(
            center=(cx, cy), r=radius,
            theta1=angle, theta2=angle + sweep,
            facecolor=color, edgecolor='white', linewidth=0.2
        ))
        if count > 0:
            if single_color:
                lx, ly = cx, cy
            else:
                mid = np.radians(angle + sweep / 2)
                lx = cx + radius * 0.6 * np.cos(mid)
                ly = cy + radius * 0.6 * np.sin(mid)
            ax.text(lx, ly, str(count),
                    ha='center', va='center',
                    fontsize=30, fontweight='bold', color='black')
        angle += sweep


def _abbreviate(name, max_len=25):
    """Break long names across two lines for axis labels."""
    parts = name.split('_')
    if len(name) <= 15:
        return name
    if len(parts) > 1 and len(parts[0]) >= 8:
        return f"{parts[0]}\n{parts[1]}"
    if name == "Nickel_sulfate_hexahydrate":
        return f"{'_'.join(parts[:2])}\n{parts[-1]}"
    if name in ("m_Intestinal_organoids", "h_Intestinal_organoids"):
        return f"{'_'.join(parts[:2])}\n{parts[-1]}"
    return f"{name[:max_len]}\n{name[max_len:]}"


def _legend_elements(marker_size=12):
    return [
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=CLUSTER_OUTCOME_COLORS['main_cluster'],
                   markersize=marker_size, label='Main Cluster'),
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=CLUSTER_OUTCOME_COLORS['small_cluster'],
                   markersize=marker_size, label='Small Cluster (n=2)'),
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=CLUSTER_OUTCOME_COLORS['singleton'],
                   markersize=marker_size, label='Singleton'),
    ]


# ===========================================================================
# PLOT BUILDERS
# ===========================================================================

def _draw_grid(ax, n_cols, n_rows, bold_boundaries=None):
    """Draw grid lines; bold_boundaries is a set of row indices for species dividers."""
    for i in range(n_cols + 1):
        ax.axvline(i, color='lightgray', alpha=0.7, linewidth=1)
    for i in range(n_rows + 1):
        if bold_boundaries and i in bold_boundaries:
            ax.axhline(i, color='black', alpha=0.9, linewidth=3.0)
        else:
            ax.axhline(i, color='lightgray', alpha=0.7, linewidth=1)


def _fill_cells(ax, species_models_subset, group_compounds, matrix_data, pie_radius):
    """Render all pie-chart cells for one group."""
    n_rows = len(species_models_subset)
    cells  = 0
    for row_idx, sm in enumerate(species_models_subset):
        for col_idx, compound in enumerate(group_compounds):
            key = (sm, compound)
            if key not in matrix_data:
                continue
            d = matrix_data[key]
            if d['total'] == 0:
                continue
            cells += 1
            cx = col_idx + 0.5
            cy = n_rows - 1 - row_idx + 0.5

            props, cols, cnts = [], [], []
            for cat, color_key in [('clustered',     'main_cluster'),
                                    ('small_clusters','small_cluster'),
                                    ('singletons',    'singleton')]:
                if d[cat] > 0:
                    props.append(d[cat] / d['total'])
                    cols.append(CLUSTER_OUTCOME_COLORS[color_key])
                    cnts.append(d[cat])

            _pie_patch(ax, cx, cy, pie_radius, props, cols, cnts)
    return cells


def _make_group_plot(group_name, species_models_subset, group_compounds,
                     matrix_data, output_path):
    """Render and save the plot for one species group."""
    cell_size  = 2.5
    n_cols     = len(group_compounds)
    n_rows     = len(species_models_subset)
    fig_width  = max(10, n_cols * cell_size + 4)
    fig_height = max(8,  n_rows * cell_size + 4)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    cells = _fill_cells(ax, species_models_subset, group_compounds,
                        matrix_data, pie_radius=0.45)
    print(f"  {group_name}: drew {cells} cells")

    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, n_rows)
    ax.set_aspect('equal', adjustable='box')

    # Column labels (compound names, rotated)
    for col_idx, compound in enumerate(group_compounds):
        ax.text(col_idx + 0.5, n_rows + 0.1,
                _abbreviate(compound, 25),
                ha='center', va='bottom', fontsize=50, rotation=90)

    ax.set_xticks([])
    ax.set_xticklabels([])

    # Row labels (model names)
    ax.set_yticks([i + 0.5 for i in range(n_rows)])
    model_labels = []
    for sm in reversed(species_models_subset):
        if sm == "C.elegans":
            model_labels.append("C. elegans")
        else:
            parts   = sm.split('_', 1)
            species = parts[0]
            model   = parts[1] if len(parts) > 1 else sm
            label   = f"{species} {model}" if group_name == "Other" else model
            model_labels.append(label)

    ax.set_yticklabels([_abbreviate(l, 20) for l in model_labels], fontsize=50)
    for tick in ax.get_yticklabels():
        tick.set_multialignment('center')

    _draw_grid(ax, n_cols, n_rows)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, top=False)

    ax.legend(handles=_legend_elements(marker_size=30),
              loc='upper center', bbox_to_anchor=(0.5, -0.025),
              ncol=3, frameon=True, fontsize=40)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.2)
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path}")


# ===========================================================================
# PUBLIC API
# ===========================================================================

def create_species_compound_matrix(species_models, compounds, matrix_data,
                                    output_path):
    """
    Create three separate PDF plots (Human / Mouse / Other Species).

    Output files are derived from output_path:
      <base>_human.<ext>
      <base>_mouse.<ext>
      <base>_other_species.<ext>

    Parameters
    ----------
    species_models : list[str]
    compounds      : list[str]
    matrix_data    : dict  (from collect_species_compound_data)
    output_path    : str   Base path; extension determines format.
    """
    print("Creating 3 separate plots: Human, Mouse, Other Species...")

    groups = {'Human': [], 'Mouse': [], 'Other': []}
    for sm in species_models:
        if sm == "C.elegans":
            groups['Other'].append(sm)
        else:
            species = sm.split('_')[0]
            if species == 'Human':
                groups['Human'].append(sm)
            elif species == 'Mouse':
                groups['Mouse'].append(sm)
            else:
                groups['Other'].append(sm)

    base, ext = os.path.splitext(output_path)
    ext = ext or '.pdf'

    suffix_map = {'Human': 'human', 'Mouse': 'mouse', 'Other': 'other_species'}

    for group_name, subset in groups.items():
        if not subset:
            print(f"  No models for {group_name}, skipping.")
            continue

        # Only include compounds that have data for this group
        group_compounds = sorted(
            c for c in compounds
            if any(matrix_data.get((sm, c), {}).get('total', 0) > 0
                   for sm in subset)
        )
        if not group_compounds:
            print(f"  No compound data for {group_name}, skipping.")
            continue

        out = f"{base}_{suffix_map[group_name]}{ext}"
        print(f"\n--- {group_name}: {len(subset)} models, "
              f"{len(group_compounds)} compounds ---")
        _make_group_plot(group_name, subset, group_compounds, matrix_data, out)

    print("\nAll 3 plots created.")


# ===========================================================================
# CLI
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate species-model vs compound clustering matrix.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--main_output_dir", required=True,
                        help="Directory containing main_clusters/, small_clusters/, singletons/.")
    parser.add_argument("--output_path", required=True,
                        help="Base output path (e.g. results/matrix.pdf). "
                             "Three files will be created with species suffixes.")
    parser.add_argument("--abbreviation_file",
                        help="TSV with 'compound' and 'acronym' columns (optional).")
    parser.add_argument("--exclude_species", nargs='+', default=[],
                        help="Species to exclude (currently only applies to combined "
                             "plot mode, not used in default 3-plot mode).")
    args = parser.parse_args()

    acronym_mapping = (load_acronym_mapping(args.abbreviation_file)
                       if args.abbreviation_file else {})

    print("=" * 60)
    print("SPECIES-MODEL vs COMPOUND MATRIX")
    print("=" * 60)
    print(f"Input dir:   {args.main_output_dir}")
    print(f"Output base: {args.output_path}")
    print(f"Acronym mappings: {len(acronym_mapping)}")

    species_compound_data = collect_species_compound_data(
        args.main_output_dir, acronym_mapping)

    if not species_compound_data:
        print("No data found — exiting.")
        return

    species_models, compounds, matrix_data = organize_matrix_data(
        species_compound_data)

    create_species_compound_matrix(
        species_models, compounds, matrix_data, args.output_path)

    print("\n" + "=" * 60)
    print("MATRIX GENERATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
