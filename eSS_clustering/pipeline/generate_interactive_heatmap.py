import pandas as pd
import json
import os
import sys
import argparse
from sklearn.metrics.pairwise import cosine_similarity

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.mutation_type import get_config


def generate_interactive_heatmap(esignature_path, cosmic_path, image_dir,
                                  output_path, threshold, cfg,
                                  output_high_path=None, output_low_path=None):
    """
    Generates an interactive HTML heatmap and optionally saves filtered profiles.

    cfg : MutationTypeConfig — used for axis label
    """
    print("--- Starting Interactive Heatmap Generation ({}) ---".format(cfg.name))

    # === Step 1: Load data ===
    try:
        print(f"Loading eSignature profiles from: {esignature_path}")
        A_df = pd.read_csv(esignature_path, sep='\t', index_col=0)

        print(f"Loading COSMIC profiles from: {cosmic_path}")
        B_df = pd.read_csv(cosmic_path, sep='\t', index_col=0)
    except FileNotFoundError as e:
        print(f"❌ ERROR: Input file not found.\nDetails: {e}")
        return

    A_t = A_df.T
    B_t = B_df.T

    # === Step 2: Cosine similarity & filter ===
    print("Calculating cosine similarity matrix...")
    sim_matrix = cosine_similarity(A_t, B_t)
    sim_df     = pd.DataFrame(sim_matrix, index=A_t.index, columns=B_t.index)

    print(f"Filtering eSignatures with similarity > {threshold}...")
    rows_to_keep = (sim_df >= threshold).any(axis=1)

    high_sigs = rows_to_keep[rows_to_keep].index
    low_sigs  = rows_to_keep[~rows_to_keep].index

    df_high = A_df[high_sigs]
    df_low  = A_df[low_sigs]

    if output_high_path:
        df_high.to_csv(output_high_path, sep='\t')
        print(f"✅ Saved {len(df_high.columns)} high-similarity profiles → "
              f"{os.path.abspath(output_high_path)}")

    if output_low_path:
        df_low.to_csv(output_low_path, sep='\t')
        print(f"✅ Saved {len(df_low.columns)} low-similarity profiles → "
              f"{os.path.abspath(output_low_path)}")

    filtered_sim_df = sim_df[rows_to_keep]
    print(f"✅ Kept {len(filtered_sim_df)} of {len(sim_df)} eSignatures for the heatmap.")

    # === Step 3: JSON for embedding ===
    sim_json = {
        "index":   list(filtered_sim_df.index),
        "columns": list(filtered_sim_df.columns),
        "data":    filtered_sim_df.values.tolist()
    }
    sim_json_str = json.dumps(sim_json)

    # === Step 4: Relative image paths ===
    output_dir = os.path.dirname(output_path) or '.'
    relative_image_base = os.path.relpath(image_dir, output_dir)
    esignature_img_path = os.path.join(relative_image_base, "cluster_plots").replace("\\", "/")
    cosmic_img_path     = os.path.join(relative_image_base, "cosmic_plots").replace("\\", "/")

    # === Step 5: Axis label from mutation type ===
    x_axis_label = f"COSMIC {cfg.name} Signatures"

    # === Step 6: HTML template ===
    html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Cosine Similarity Explorer</title>
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
  <style>
    html, body {
      height: 100%; margin: 0; padding: 0; display: block; font-family: Arial, sans-serif;
    }
    body {
      display: flex; flex-direction: row-reverse; gap: 40px; padding: 5px;
    }
    #heatmap { width: 65%; transition: all 0.3s ease; }
    #heatmap.shrunk { width: 30%; }
    #details {
      width: 30%; display: flex; flex-direction: column; justify-content: center;
      overflow-y: auto; max-height: 100vh; transition: all 0.3s ease;
      padding-right: 12px; scrollbar-gutter: stable;
    }
    #details.expanded { width: 70%; }
    img { width: 100%; object-fit: contain; border: 1px solid #ccc; margin-bottom: 20px; }
    .badge {
      padding: 3px 10px; border-radius: 8px; font-weight: bold;
      font-size: 0.9rem; margin-left: 10px;
    }
    #zoom-btn {
      margin-bottom: 20px; padding: 6px 12px; font-size: 0.9rem;
      cursor: pointer; align-self: flex-start;
    }
  </style>
</head>
<body>

<div id="details">
  <button id="zoom-btn">🔍 Expand for Comparison</button>
  <div><strong>eSignature:</strong> <span id="esig-label"></span></div>
  <img id="esig-img" src="">
  <div><strong>COSMIC:</strong> <span id="cosmic-label"></span></div>
  <img id="cosmic-img" src="">
</div>

<div id="heatmap"></div>

<script>
  const simData = __SIM_JSON_PLACEHOLDER__;
  const esigImgDir   = "__ESIG_IMG_PATH_PLACEHOLDER__";
  const cosmicImgDir = "__COSMIC_IMG_PATH_PLACEHOLDER__";
  const xAxisLabel   = "__X_AXIS_LABEL_PLACEHOLDER__";

  const xLabels = simData.columns;
  const yLabels = simData.index;
  const zValues = simData.data;

  const hoverText = zValues.map((rowVals, i) =>
    rowVals.map((val, j) =>
      `eSignatures: ${yLabels[i]}<br>COSMIC: ${xLabels[j]}<br>Cosine Similarity: ${val.toFixed(3)}`
    )
  );

  const zmin = 0.8, zmax = 1.0, range = zmax - zmin;
  const p1 = (0.85 - zmin) / range;
  const p2 = (0.90 - zmin) / range;
  const p3 = (0.95 - zmin) / range;

  const discreteColorscale = [
    [0, '#440154'], [0.000001, '#440154'],
    [0.000001, '#3b528b'], [p1, '#3b528b'],
    [p1, '#21918c'], [p2, '#21918c'],
    [p2, '#5ec962'], [p3, '#5ec962'],
    [p3, '#fde725'], [1.0, '#fde725']
  ];

  const heatmap = {
    z: zValues, x: xLabels, y: yLabels,
    type: 'heatmap', text: hoverText, hoverinfo: 'text',
    colorscale: discreteColorscale,
    zmin: zmin, zmax: zmax,
    colorbar: {
      title: { text: 'Cosine<br>Similarity', side: 'top' },
      y: 0.5, yanchor: 'middle', len: 0.9,
      tickvals: [0.80, 0.85, 0.90, 0.95, 1.0],
      ticktext: ['0.80', '0.85', '0.90', '0.95', '1.0']
    }
  };

  const layout = {
    title: 'Cosine Similarity: eSignatures vs COSMIC',
    xaxis: { title: xAxisLabel, tickangle: 90 },
    yaxis: { title: { text: 'eSignature Clusters', standoff: 20 }, automargin: true },
    margin: { l: 150, r: 20, t: 80, b: 120 },
    autosize: true
  };

  Plotly.newPlot('heatmap', [heatmap], layout, {responsive: true});

  document.getElementById('heatmap').on('plotly_hover', function(data) {
    const pt = data.points[0];
    const rowIndex = yLabels.indexOf(pt.y);
    const colIndex = xLabels.indexOf(pt.x);
    const rowShape = { type: 'rect', xref: 'x', yref: 'y', x0: -0.5, x1: xLabels.length - 0.5, y0: rowIndex - 0.5, y1: rowIndex + 0.5, fillcolor: 'rgba(255,255,255,0.2)', line: { width: 0 }, layer: 'above' };
    const colShape = { type: 'rect', xref: 'x', yref: 'y', x0: colIndex - 0.5, x1: colIndex + 0.5, y0: -0.5, y1: yLabels.length - 0.5, fillcolor: 'rgba(255,255,255,0.2)', line: { width: 0 }, layer: 'above' };
    Plotly.relayout('heatmap', { shapes: [rowShape, colShape] });
  });

  document.getElementById('heatmap').on('plotly_unhover', function() {
    Plotly.relayout('heatmap', { shapes: [] });
  });

  document.getElementById('heatmap').on('plotly_click', function(data) {
    const pt = data.points[0];
    const esig = pt.y, cosmic = pt.x, cosVal = pt.z;
    document.getElementById('esig-label').innerHTML =
      `${esig} <span class="badge" style="background-color:${getColor(cosVal)};color:${getTextColor(cosVal)}">${cosVal.toFixed(2)}</span>`;
    document.getElementById('cosmic-label').textContent = cosmic;
    document.getElementById('esig-img').src   = `${esigImgDir}/${esig}.png`;
    document.getElementById('cosmic-img').src = `${cosmicImgDir}/${cosmic}.png`;
  });

  function getColor(cos) {
    if (cos <= 0.80) return '#440154';
    if (cos <= 0.85) return '#3b528b';
    if (cos <= 0.90) return '#21918c';
    if (cos <= 0.95) return '#5ec962';
    return '#fde725';
  }

  function getTextColor(cos) { return cos <= 0.90 ? 'white' : 'black'; }

  const zoomBtn    = document.getElementById("zoom-btn");
  const details    = document.getElementById("details");
  const heatmapDiv = document.getElementById("heatmap");

  zoomBtn.addEventListener("click", () => {
    details.classList.toggle("expanded");
    heatmapDiv.classList.toggle("shrunk");
    zoomBtn.textContent = details.classList.contains("expanded")
      ? "❌ Collapse View" : "🔍 Expand for Comparison";
    setTimeout(function() { Plotly.Plots.resize(heatmapDiv); }, 300);
  });
</script>
</body>
</html>
"""

    final_html = html_template.replace("__SIM_JSON_PLACEHOLDER__",      sim_json_str)
    final_html = final_html.replace("__ESIG_IMG_PATH_PLACEHOLDER__",    esignature_img_path)
    final_html = final_html.replace("__COSMIC_IMG_PATH_PLACEHOLDER__",  cosmic_img_path)
    final_html = final_html.replace("__X_AXIS_LABEL_PLACEHOLDER__",     x_axis_label)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(final_html)

    print(f"\n✅ HTML file created: {os.path.abspath(output_path)}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate an interactive HTML cosine similarity heatmap.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--esignature_profiles", required=True)
    parser.add_argument("--cosmic_profiles",     required=True)
    parser.add_argument("--image_directory",     required=True)
    parser.add_argument("--output_file",         required=True)
    parser.add_argument("--threshold",  type=float, default=0.85)
    parser.add_argument("--output_high_similarity")
    parser.add_argument("--output_low_similarity")
    parser.add_argument("--mutation_type", default="SBS",
                        choices=["SBS", "DBS", "ID"],
                        help="Mutation type — used for the heatmap x-axis label.")
    args = parser.parse_args()

    cfg = get_config(args.mutation_type)

    generate_interactive_heatmap(
        args.esignature_profiles,
        args.cosmic_profiles,
        args.image_directory,
        args.output_file,
        args.threshold,
        cfg,
        args.output_high_similarity,
        args.output_low_similarity,
    )


if __name__ == "__main__":
    main()
