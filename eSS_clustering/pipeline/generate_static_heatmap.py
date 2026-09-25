#!/usr/bin/env python

import argparse
import sys
import os
import math
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.mutation_type import get_config


def generate_static_heatmap(esignature_path, cosmic_path, output_path, cfg):
    """
    Generates three static heatmaps:
      1. All COSMIC signatures
      2. Artifacts removed
      3. Artifacts kept but low-similarity COSMIC signatures filtered out

    cfg : MutationTypeConfig  — controls artifact list and axis labels
    """

    def _create_and_save_plot(sim_df, file_path, plot_title):
        print(f"\n--- Generating Plot: {plot_title} ---")

        annot_mask = sim_df.applymap(lambda x: f"{x:.2f}" if x >= 0.8 else "")

        colors     = ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725']
        boundaries = [0, 0.8, 0.845, 0.895, 0.945, 1.0]
        cmap       = ListedColormap(colors)
        norm       = BoundaryNorm(boundaries, cmap.N, clip=True)

        fig, ax = plt.subplots(
            figsize=(max(12, sim_df.shape[1] * 0.4),
                     max(8,  sim_df.shape[0] * 0.4))
        )

        sns.heatmap(ax=ax, data=sim_df, annot=annot_mask, fmt="",
                    cmap=cmap, norm=norm, linewidths=0,
                    annot_kws={"fontsize": 8}, cbar=False)

        label_fontsize = 18
        tick_fontsize  = 14

        cbar_ax   = fig.add_axes([0.95, 0.11, 0.02, 0.25])
        cbar_ticks = [0, 0.8, 0.845, 0.895, 0.945, 1.0]
        cbar = fig.colorbar(ax.collections[0], cax=cbar_ax,
                            orientation='vertical', ticks=cbar_ticks)
        cbar.ax.set_yticklabels(['0.00', '0.80', '0.85', '0.90', '0.95', '1.00'])
        cbar.ax.set_title("Cosine\nSimilarity", fontsize=label_fontsize, pad=20)
        cbar.ax.tick_params(labelsize=tick_fontsize)

        # Colour annotation text for readability
        valid_texts = [t for t in ax.texts if t.get_text()]
        for text_obj in valid_texts:
            x, y  = text_obj.get_position()
            val   = sim_df.iloc[int(y), int(x)]
            rgba  = cmap(norm(val))
            brightness = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
            text_obj.set_color("black" if brightness > 0.6 else "white")

        matched_esigs  = sim_df.index[sim_df.max(axis=1) >= 0.845]
        de_novo_esigs  = sim_df.index[sim_df.max(axis=1) <  0.845]
        print(f"eSignatures with matches >=0.85: {list(matched_esigs)}")
        print(f"De novo eSignatures:             {list(de_novo_esigs)}")

        # Dynamic axis labels based on mutation type
        x_label = f"COSMIC v3.5 {cfg.name} Signatures"
        ax.set_xlabel(x_label, fontsize=label_fontsize,
                      fontweight='bold', labelpad=20)
        ax.set_ylabel("eSignature clusters", fontsize=label_fontsize,
                      fontweight='bold', labelpad=20)
        ax.tick_params(axis='x', labelsize=tick_fontsize, rotation=90)
        ax.tick_params(axis='y', labelsize=tick_fontsize)

        for label in ax.get_yticklabels():
            if label.get_text() in matched_esigs:
                label.set_fontweight('bold')

        fig.text(-0.06, 0.09, r"$\mathbf{bold}$: indicates",
                 fontsize=12, ha='left', va='center')
        fig.text(-0.04, 0.08, r"similarity ≥0.85",
                 fontsize=12, ha='left', va='center')

        print(f"--- Saving plot to {file_path} ---")
        plt.savefig(file_path, dpi=400, bbox_inches='tight')
        plt.close()

    # --- Load data ---
    print("--- Loading Data ---")
    try:
        a_df = pd.read_csv(esignature_path, index_col=0, sep='\t')
        b_df = pd.read_csv(cosmic_path,     index_col=0, sep='\t')
    except FileNotFoundError as e:
        print(f"❌ ERROR: Input file not found.\nDetails: {e}")
        return

    a_t = a_df.T
    b_t = b_df.T

    base, ext = os.path.splitext(output_path)

    # --- Version 1: all signatures ---
    sim_df_all = pd.DataFrame(
        cosine_similarity(a_t, b_t),
        index=a_t.index, columns=b_t.index
    )
    _create_and_save_plot(sim_df_all, output_path,
                          f"Cosine Similarity ({cfg.name} — All Signatures)")

    # --- Version 2: artifacts removed ---
    # Use the artifact list from the config; skip sigs not in this COSMIC file
    artifact_signatures   = cfg.artifact_signatures
    artifacts_in_df       = [s for s in artifact_signatures if s in b_df.columns]
    b_df_no_artifacts     = b_df.drop(columns=artifacts_in_df)
    b_t_no_artifacts      = b_df_no_artifacts.T

    if artifacts_in_df:
        print(f"\nRemoving {len(artifacts_in_df)} artifact signatures: {artifacts_in_df}")
    else:
        print(f"\nNo artifact signatures found in COSMIC file for {cfg.name}.")

    sim_df_no_art = pd.DataFrame(
        cosine_similarity(a_t, b_t_no_artifacts),
        index=a_t.index, columns=b_t_no_artifacts.index
    )
    output_no_art = f"{base}_no_artifacts{ext}"
    _create_and_save_plot(sim_df_no_art, output_no_art,
                          f"Cosine Similarity ({cfg.name} — Artifacts Removed)")

    # --- Version 3: keep artifacts, remove low-similarity COSMIC sigs ---
    high_sim_cols = sim_df_all.columns[sim_df_all.max(axis=0) > 0.8]
    sim_df_high   = sim_df_all[high_sim_cols]

    print(f"\n--- Filtering Summary ---")
    print(f"Original COSMIC signatures:               {len(b_df.columns)}")
    print(f"After removing artifacts:                 {len(b_df_no_artifacts.columns)}")
    print(f"Keeping artifacts, removing low-sim sigs: {len(sim_df_high.columns)}")
    print(f"Removed {len(b_df.columns) - len(sim_df_high.columns)} low-similarity signatures")

    output_high = f"{base}_high_similarity_only{ext}"
    _create_and_save_plot(sim_df_high, output_high,
                          f"Cosine Similarity ({cfg.name} — High Similarity Only)")

    print("\n✅ Done! All three heatmaps generated.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate three static cosine similarity heatmaps.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--esignature_profiles", required=True,
                        help="Path to the eSignature profiles TSV.")
    parser.add_argument("--cosmic_profiles", required=True,
                        help="Path to the COSMIC profiles file.")
    parser.add_argument("--output_file", required=True,
                        help="Base output path (e.g. heatmap.pdf). Two additional "
                             "files are created with '_no_artifacts' and "
                             "'_high_similarity_only' suffixes.")
    parser.add_argument("--mutation_type", default="SBS",
                        choices=["SBS", "DBS", "ID"],
                        help="Mutation type — controls artifact list and axis labels.")
    args = parser.parse_args()

    cfg = get_config(args.mutation_type)
    generate_static_heatmap(
        args.esignature_profiles,
        args.cosmic_profiles,
        args.output_file,
        cfg,
    )


if __name__ == "__main__":
    main()
