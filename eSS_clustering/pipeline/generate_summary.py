#!/usr/bin/env python3
"""
Mutational Signature Clustering Report Generator

Generates HTML reports for all three cluster types (main, small, singletons).
Imports shared naming utilities from utils/ rather than maintaining local copies.
"""

import argparse
import os
import shutil
import sys
from collections import Counter

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.naming import (
    load_acronym_mapping,
    apply_acronym_mapping,
    extract_sample_prefix,
    apply_celegans_collapsing,
)
from utils.mutation_type import get_config

# ---------------------------------------------------------------------------
# Colour palette (consistent across all visualisations)
# ---------------------------------------------------------------------------
SIGNATURE_COLORS = [
    '#8DD3C7', '#FFD700', '#BEBADA', '#FB8072',
    '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5',
]


# ==============================================================================
# DATA LOADING HELPERS
# ==============================================================================

def load_cluster_data(cluster_type_dir, cluster_type):
    """Load <type>_clusters_summary.tsv → grouped DataFrame or None."""
    path = os.path.join(cluster_type_dir, f"{cluster_type}_clusters_summary.tsv")
    if not os.path.exists(path):
        print(f"No {cluster_type} cluster summary at: {path}")
        return None
    try:
        df = pd.read_csv(path, sep='\t')
        if not {'Cluster', 'Color', 'Sample'}.issubset(df.columns):
            print(f"Warning: missing required columns in {path}")
            return None
        summary = (
            df.groupby(['Cluster', 'Color'])
            .agg(Count=('Sample', 'count'),
                 Samples=('Sample', lambda x: "\n".join(x)))
            .reset_index()
            .set_index('Cluster')
            .sort_index(key=lambda x: pd.to_numeric(x, errors='coerce'))
        )
        print(f"Loaded {len(summary)} {cluster_type} clusters")
        return summary
    except Exception as e:
        print(f"Error loading {cluster_type} clusters: {e}")
        return None


def load_singleton_data(singletons_dir):
    """Load singleton_small_clusters_ordered.tsv → cluster-style DataFrame or None."""
    path = os.path.join(singletons_dir, "singleton_small_clusters_ordered.tsv")
    if not os.path.exists(path):
        print(f"No singleton file at: {path}")
        return None
    try:
        df = pd.read_csv(path, sep='\t', index_col=0)
        rows = [{'Cluster': i, 'Color': '#808080', 'Count': 1, 'Samples': col}
                for i, col in enumerate(df.columns, 1)]
        summary = pd.DataFrame(rows).set_index('Cluster')
        print(f"Loaded {len(summary)} singletons")
        return summary
    except Exception as e:
        print(f"Error loading singletons: {e}")
        return None


def find_plots_directory(cluster_type_dir, cluster_type):
    """Return path to the PNG plots directory for a cluster type, or None."""
    for pattern in [
        f"{cluster_type}_clusters_hierarchical_plots",
        f"singleton_samples_plots",
        f"{cluster_type}_plots",
        "plots",
    ]:
        d = os.path.join(cluster_type_dir, pattern)
        if os.path.exists(d):
            pngs = [f for f in os.listdir(d)
                    if f.endswith('.png') and ('cluster' in f or 'singleton' in f)]
            if pngs:
                print(f"Found {len(pngs)} plots in {d}")
                return d
    # Check directly in cluster dir
    pngs = [f for f in os.listdir(cluster_type_dir)
            if f.endswith('.png') and ('cluster' in f or 'singleton' in f)]
    if pngs:
        return cluster_type_dir
    print(f"No plots found for {cluster_type} in {cluster_type_dir}")
    return None


def load_original_sample_mappings(mappings_file_path):
    """Load applied_sample_name_mappings.tsv → {standardized: [original, ...]}."""
    if not mappings_file_path or not os.path.exists(mappings_file_path):
        print(f"Original sample mappings not found: {mappings_file_path}")
        return {}
    try:
        df = pd.read_csv(mappings_file_path, sep='\t')
        if not {'original_name', 'standardized_name'}.issubset(df.columns):
            print("Warning: expected 'original_name' and 'standardized_name' columns")
            return {}
        result = {}
        for _, row in df.iterrows():
            std  = str(row['standardized_name']).strip()
            orig = str(row['original_name']).strip()
            result.setdefault(std, []).append(orig)
        print(f"Loaded reverse mappings for {len(result)} standardised names")
        return result
    except Exception as e:
        print(f"Error loading original sample mappings: {e}")
        return {}


# ==============================================================================
# SHARED PREFIX PROCESSING
# ==============================================================================

def _collapse_and_map_prefixes(sample_lines, acronym_mapping):
    """
    Extract prefixes, apply C. elegans collapsing and acronym mapping.
    Returns list of collapsed/mapped prefix strings (same length as sample_lines).
    """
    result = []
    for sample in sample_lines:
        prefix = extract_sample_prefix(sample)
        # C. elegans collapsing via shared util
        prefix = apply_celegans_collapsing(prefix)
        # Acronym mapping via shared util
        if acronym_mapping:
            prefix = apply_acronym_mapping(prefix, acronym_mapping)
        result.append(prefix)
    return result


# ==============================================================================
# SVG PIE CHART & LEGEND HELPERS
# ==============================================================================

def _svg_pie(samples, size=30):
    """Return an inline SVG pie chart for use in HTML."""
    import math
    radius, center = size * 0.4, size / 2
    paths, angle = [], 0.0
    for idx, s in enumerate(samples):
        pct   = s['percentage']
        sweep = pct / 100 * 360
        a1    = math.radians(angle - 90)
        a2    = math.radians(angle + sweep - 90)
        x1, y1 = center + radius * math.cos(a1), center + radius * math.sin(a1)
        x2, y2 = center + radius * math.cos(a2), center + radius * math.sin(a2)
        large  = 1 if pct > 50 else 0
        css    = f"color-{idx % 8}"
        if pct == 100:
            paths.append(
                f'<circle cx="{center}" cy="{center}" r="{radius}" '
                f'class="{css}" stroke="white" stroke-width="1"/>'
            )
        else:
            paths.append(
                f'<path d="M {center} {center} L {x1:.2f} {y1:.2f} '
                f'A {radius} {radius} 0 {large} 1 {x2:.2f} {y2:.2f} Z" '
                f'class="{css}" stroke="white" stroke-width="0.5"/>'
            )
        angle += sweep
    return (f'<svg width="{size}" height="{size}" class="pie-chart">'
            + "".join(paths) + '</svg>')


def _legend_html(samples):
    items = []
    for idx, s in enumerate(samples):
        items.append(
            f'<div class="legend-item">'
            f'<div class="legend-color color-{idx % 8}"></div>'
            f'<div class="legend-text">{s["name"]} ({s["percentage"]}%)</div>'
            f'</div>'
        )
    return "".join(items)


def _smart_wrap(text, max_len=20):
    if len(text) <= max_len:
        return text
    mid   = len(text) // 2
    best  = None
    mindist = float('inf')
    for char in ['_', '-', ' ']:
        for pos, c in enumerate(text):
            if c == char and pos > 10:
                d = abs(pos - mid)
                if d < mindist:
                    mindist, best = d, pos
    if best:
        return text[:best] + '<br>' + text[best + 1:]
    return text[:mid] + '<br>' + text[mid:]


def _color_css():
    return "\n".join(
        f'.color-{i} {{ fill: {c}; background-color: {c}; }}'
        for i, c in enumerate(SIGNATURE_COLORS)
    )


# ==============================================================================
# UNIFORM GRID HTML (main cluster summary without dendrogram)
# ==============================================================================

_GRID_CSS = """
body { font-family: Arial, sans-serif; margin: 0; padding: 0; background: white; line-height: 1.2; }
.cluster-container { display: flex; flex-wrap: wrap; gap: 0; width: 100%; }
.cluster-row { display: flex; align-items: center; gap: 1px; padding: 1px 0 1px 0; width: 33.333%; box-sizing: border-box; overflow: hidden; break-inside: avoid; page-break-inside: avoid; }
.signature-plot { width: 220px; height: 55px; flex-shrink: 0; background: #f9f9f9; border-radius: 2px; overflow: hidden; }
.signature-plot img { width: 100%; height: 100%; object-fit: fill; display: block; }
.pie-section { display: flex; flex-direction: column; align-items: center; gap: 2px; flex-shrink: 0; }
.pie-chart { width: 40px; height: 40px; }
.sample-count { font-size: 8px; font-weight: bold; color: #333; text-align: center; }
.legend { font-size: 6px; line-height: 1.3; flex: 1; min-width: 0; overflow: hidden; }
.legend-item { display: flex; align-items: center; margin-bottom: 2px; width: 100%; }
.legend-color { width: 4px; height: 4px; margin-right: 3px; border-radius: 50%; flex-shrink: 0; }
.legend-text { color: #444; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; flex: 1; min-width: 0; }
@media print {
  body { margin: 0; padding: 0; }
  @page { size: A4 landscape; margin: 3mm; }
  * { -webkit-print-color-adjust: exact !important; print-color-adjust: exact !important; }
}
"""


def _build_cluster_sample_data(row, acronym_mapping):
    """Return (sample_data_list, plot_image_path_template_idx) from a cluster row."""
    sample_lines = [l.strip() for l in str(row['Samples']).split('\n') if l.strip()]
    collapsed    = _collapse_and_map_prefixes(sample_lines, acronym_mapping)
    counts       = Counter(collapsed)
    total        = sum(counts.values())
    sample_data  = [
        {'name': _smart_wrap(p), 'percentage': round(c / total * 100, 1)}
        for p, c in counts.items()
    ]
    return sample_data


def create_html_summary_figure_uniform(cluster_summary, plots_dir,
                                       output_dir, acronym_mapping=None,
                                       plot_prefix="main_cluster",
                                       output_filename="cluster_summary.html",
                                       title="Cluster Summary"):
    """Uniform 3-column grid of cluster signature plots (main or small)."""
    print(f"Creating uniform HTML summary for {len(cluster_summary)} clusters "
          f"(plot_prefix={plot_prefix!r})...")
    rows_html = []
    output_html_path = os.path.join(output_dir, output_filename)
    for i, (_, row) in enumerate(cluster_summary.iterrows(), 1):
        sample_data = _build_cluster_sample_data(row, acronym_mapping)
        plot_path   = os.path.join(plots_dir, f"{plot_prefix}-{i:02d}.png")
        # Compute relative path from HTML to image (like generate_html_report does)
        try:
            rel_plot_path = os.path.relpath(plot_path, os.path.dirname(output_html_path))
        except ValueError:
            rel_plot_path = plot_path
        img_html    = (f'<img src="{rel_plot_path}" alt="Cluster {i}">'
                       if os.path.exists(plot_path) else
                       f'<div style="color:#999;font-size:10px">C{i}</div>')
        rows_html.append(f'''
        <div class="cluster-row">
          <div class="signature-plot">{img_html}</div>
          <div class="pie-section">{_svg_pie(sample_data, size=40)}
            <div class="sample-count">n={row["Count"]}</div></div>
          <div class="legend">{_legend_html(sample_data)}</div>
        </div>''')

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<title>{title}</title>
<style>{_GRID_CSS}\n{_color_css()}</style></head><body>
<div class="cluster-container">{"".join(rows_html)}</div>
</body></html>"""

    out = os.path.join(output_dir, output_filename)
    with open(out, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"Uniform summary → {out}")
    return out


# ==============================================================================
# INTERACTIVE REPORT WITH DENDROGRAM
# ==============================================================================

def generate_html_report(df_path, dendrogram_svg_path, plots_dir,
                         output_file_path, acronym_mapping=None):
    """
    Full interactive report: zoomable dendrogram + per-cluster details.
    Returns cluster_summary DataFrame.
    """
    print(f"\nGenerating HTML report: {output_file_path}")
    df = pd.read_csv(df_path, sep='\t')
    cluster_summary = (
        df.groupby(['Cluster', 'Color'])
        .agg(Count=('Sample', 'count'),
             Samples=('Sample', lambda x: "\n".join(x)))
        .reset_index().set_index('Cluster')
        .sort_index(key=lambda x: pd.to_numeric(x, errors='coerce'))
    )

    with open(dendrogram_svg_path, 'r', encoding='utf-8') as f:
        svg_content = f.read().replace('`', '\\`')

    # Build HTML
    parts = [f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Cluster Summary Report</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;padding:0}}
.dendrogram-container{{position:sticky;top:0;z-index:10;background:white;
  border-bottom:2px solid #ccc;overflow:auto;height:350px;padding:20px}}
.zoom-buttons-wrapper{{margin:20px}}
.zoom-buttons button{{margin-right:5px;padding:6px 10px}}
details{{margin:20px;border:1px solid #ccc;border-radius:4px;padding:10px}}
summary{{font-weight:bold;font-size:16px;cursor:pointer}}
.cluster-container{{display:flex;flex-direction:row;gap:20px;margin-top:10px}}
.sample-list{{flex:1;max-height:300px;overflow-y:auto;white-space:pre-wrap;
  font-family:monospace;background:#f9f9f9;padding:10px;border:1px solid #eee}}
.plot{{flex:1;display:flex;justify-content:center;align-items:center;padding:10px}}
.plot img{{max-width:100%;max-height:400px;border:1px solid #ddd}}
#svg-wrapper svg{{display:block}}
</style></head><body>
<div class="zoom-buttons-wrapper">
  <div class="zoom-buttons">
    <button onclick="zoomIn()">+</button>
    <button onclick="zoomOut()">&minus;</button>
    <button onclick="resetZoom()">Reset</button>
  </div>
</div>
<div class="dendrogram-container" id="dendrogram-container">
  <div class="dendrogram-wrapper"><div id="svg-wrapper"></div></div>
</div>
<script>
let zoomLevel=1.0;const ZOOM_FACTOR=1.5;let container;const sampleXMap={{}};
window.addEventListener('load',function(){{
  container=document.getElementById("dendrogram-container");
  const svgRaw=`{svg_content}`;
  const parser=new DOMParser();const svgDoc=parser.parseFromString(svgRaw,"image/svg+xml");
  const svgElement=svgDoc.documentElement;
  document.getElementById("svg-wrapper").appendChild(svgElement);
  setTimeout(()=>{{container.scrollTop=container.scrollHeight*0.1;}},50);
  const allTexts=svgElement.querySelectorAll("g[id^='text_']");
  allTexts.forEach(group=>{{
    const comment=Array.from(group.childNodes).find(n=>n.nodeType===Node.COMMENT_NODE);
    const label=comment?.textContent?.trim();
    const gEl=group.querySelector("g[transform]");
    const transform=gEl?.getAttribute("transform");
    const match=transform?.match(/translate\\(([^\\"]+)\\s/);
    if(label&&match){{const x=parseFloat(match[1]);if(!isNaN(x))sampleXMap[label]=x;}}
  }});
  updateZoom();
  const svg=document.querySelector('#svg-wrapper svg');
  if(svg)svg.addEventListener('click',function(e){{zoomAtCursor(e,!e.shiftKey);}});
}});
function updateZoom(){{
  const svg=document.querySelector('#svg-wrapper svg');
  if(svg){{svg.style.transform=`scale(${{zoomLevel}})`;svg.style.transformOrigin="top left";
    const wrapper=document.getElementById("svg-wrapper");
    const bbox=svg.getBBox();wrapper.style.minWidth=(bbox.width*zoomLevel+2000)+"px";}}
}}
function zoomAtCursor(event,zIn=true){{
  const rect=container.getBoundingClientRect();const old=zoomLevel;
  zoomLevel*=zIn?ZOOM_FACTOR:1/ZOOM_FACTOR;const ratio=zoomLevel/old;
  const ox=event.clientX-rect.left+container.scrollLeft;
  const oy=event.clientY-rect.top+container.scrollTop;
  updateZoom();
  container.scrollLeft=ox*ratio-(event.clientX-rect.left);
  container.scrollTop=oy*ratio-(event.clientY-rect.top);
}}
function zoomIn(){{zoomAtCursor({{clientX:container.clientWidth/2,clientY:container.clientHeight/2}},true);}}
function zoomOut(){{zoomAtCursor({{clientX:container.clientWidth/2,clientY:container.clientHeight/2}},false);}}
function resetZoom(){{zoomLevel=1;updateZoom();
  setTimeout(()=>{{container.scrollLeft=0;container.scrollTop=container.scrollHeight*0.1;}},50);}}
function zoomToCluster(sampleList){{
  const positions=sampleList.map(l=>sampleXMap[l.trim()]).filter(x=>typeof x==="number");
  if(!positions.length)return;
  const mid=positions.reduce((a,b)=>a+b,0)/positions.length;
  zoomLevel=2;updateZoom();
  setTimeout(()=>{{
    container.scrollLeft=mid*zoomLevel*1.3-container.clientWidth/2;
    container.scrollTop=container.scrollHeight*0.4;
  }},50);
}}
</script>"""]

    for cluster_id, row in cluster_summary.iterrows():
        cluster_str = str(cluster_id).zfill(2)
        color  = row['Color']
        count  = row['Count']
        samples_raw = row['Samples']

        # Apply acronym mapping for display
        if acronym_mapping:
            sample_lines  = [l.strip() for l in samples_raw.strip().split('\n') if l.strip()]
            mapped_samples = []
            for s in sample_lines:
                prefix = extract_sample_prefix(s)
                mapped = apply_acronym_mapping(prefix, acronym_mapping)
                mapped_samples.append(s.replace(prefix, mapped) if prefix != mapped else s)
            samples_display = '\n'.join(mapped_samples)
        else:
            samples_display = samples_raw

        js_list = '[' + ','.join(
            f"'{s.strip()}'"
            for s in samples_raw.strip().split('\n') if s.strip()
        ) + ']'
        plot_path = os.path.join(plots_dir, f"main_cluster-{cluster_str}.png")
        try:
            rel_path = os.path.relpath(plot_path, os.path.dirname(output_file_path))
        except ValueError:
            rel_path = plot_path
        plot_tag = (f'<img class="plot-image" src="{rel_path}" '
                    f'alt="Cluster {cluster_str}">'
                    if os.path.exists(plot_path) else "<em>Plot not available</em>")

        parts.append(f"""
<details>
  <summary style="color:{color};cursor:pointer;" onclick="zoomToCluster({js_list})">
    Cluster {cluster_id} (n={count})
  </summary>
  <div class="cluster-container">
    <div class="sample-list">{samples_display}</div>
    <div class="plot">{plot_tag}</div>
  </div>
</details>""")

    parts.append("</body></html>")
    with open(output_file_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))
    print(f"Interactive report → {output_file_path}")
    return cluster_summary


# ==============================================================================
# DENDROGRAM SUMMARY (main clusters only)
# ==============================================================================

def _build_dendrogram_clusters_html(cluster_summary, plots_dir, acronym_mapping,
                                     plot_prefix="cluster"):
    """Render cluster rows for the with-dendrogram HTML."""
    rows = []
    output_dir_for_rel_path = None  # Will be set when we know where HTML is saved
    
    for i, (cluster_id, row) in enumerate(cluster_summary.iterrows(), 1):
        sample_data = _build_cluster_sample_data(row, acronym_mapping)
        plot_path   = os.path.join(plots_dir, f"{plot_prefix}-{i:02d}.png")
        img_html    = (f'<img src="{plot_path}" alt="Cluster {i}">'
                       if os.path.exists(plot_path) else
                       f'<div style="color:#999;font-size:10px">C{i}</div>')
        border = f'border: 1px solid {row["Color"]};'
        rows.append(f'''
        <div class="cluster-row" style="{border}">
          <div class="signature-plot">{img_html}</div>
          <div class="pie-section">{_svg_pie(sample_data)}
            <div class="sample-count">n={row["Count"]}</div></div>
          <div class="legend">{_legend_html(sample_data)}</div>
        </div>''')
    return rows


_DENDRO_CSS = """
body{font-family:Arial,sans-serif;margin:0;padding:0;background:white;line-height:1.2}
.dendrogram-section{width:100%;margin-bottom:10px;text-align:center}
.dendrogram-section svg{width:100%;height:auto;max-height:200px}
.cluster-container{display:flex;flex-wrap:wrap;gap:0;width:100%}
.cluster-row{display:flex;align-items:center;gap:1px;padding:1px 4px 1px 0;width:33.333%;
  min-height:25px;box-sizing:border-box;overflow:hidden;break-inside:avoid;page-break-inside:avoid}
.signature-plot{width:220px;height:55px;flex-shrink:0;background:#f9f9f9;
  border-radius:2px;overflow:hidden}
.signature-plot img{width:100%;height:100%;object-fit:fill;display:block}
.pie-section{display:flex;flex-direction:column;align-items:center;gap:2px;flex-shrink:0}
.pie-chart{width:40px;height:40px}
.sample-count{font-size:8px;font-weight:bold;color:#333;text-align:center}
.legend{font-size:6px;line-height:1.3;flex:1;min-width:0;overflow:hidden}
.legend-item{display:flex;align-items:center;margin-bottom:2px}
.legend-color{width:4px;height:4px;margin-right:3px;border-radius:50%;flex-shrink:0}
.legend-text{color:#444;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;max-width:100%}
@media print{body{margin:0;padding:0}@page{size:A4 landscape;margin:3mm}
  *{-webkit-print-color-adjust:exact !important;print-color-adjust:exact !important}}
"""


def create_html_summary_figure_with_dendrogram(cluster_summary, plots_dir,
                                               output_path,
                                               dendrogram_svg_path=None,
                                               acronym_mapping=None,
                                               plot_prefix="cluster"):
    """HTML summary with embedded dendrogram above cluster rows."""
    dendro_html = ""
    if dendrogram_svg_path and os.path.exists(dendrogram_svg_path):
        with open(dendrogram_svg_path, 'r', encoding='utf-8') as f:
            svg = f.read()
        dendro_html = f'<div class="dendrogram-section">{svg}</div>'

    # Build rows with relative paths
    rows_html = []
    for i, (cluster_id, row) in enumerate(cluster_summary.iterrows(), 1):
        sample_data = _build_cluster_sample_data(row, acronym_mapping)
        plot_path   = os.path.join(plots_dir, f"{plot_prefix}-{i:02d}.png")
        # Compute relative path from HTML to image
        try:
            rel_plot_path = os.path.relpath(plot_path, os.path.dirname(output_path))
        except ValueError:
            rel_plot_path = plot_path
        img_html = (f'<img src="{rel_plot_path}" alt="Cluster {i}">'
                    if os.path.exists(plot_path) else
                    f'<div style="color:#999;font-size:10px">C{i}</div>')
        border = f'border: 1px solid {row["Color"]};'
        rows_html.append(f'''
        <div class="cluster-row" style="{border}">
          <div class="signature-plot">{img_html}</div>
          <div class="pie-section">{_svg_pie(sample_data, size=40)}
            <div class="sample-count">n={row["Count"]}</div></div>
          <div class="legend">{_legend_html(sample_data)}</div>
        </div>''')

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<title>Cluster Summary</title>
<style>{_DENDRO_CSS}\n{_color_css()}</style></head><body>
{dendro_html}
<div class="cluster-container">{"".join(rows_html)}</div>
</body></html>"""
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"Dendrogram summary → {output_path}")
    return output_path


def create_dendrogram_summary(main_output_dir, reports_dir,
                              dendrogram_svg_path=None, acronym_mapping=None):
    """Create main_clusters_with_dendrogram.html only."""
    main_clusters_dir = os.path.join(main_output_dir, "main_clusters")
    if not os.path.exists(main_clusters_dir):
        print("main_clusters directory not found; skipping dendrogram summary.")
        return

    main_summary  = load_cluster_data(main_clusters_dir, "main")
    main_plots_dir = find_plots_directory(main_clusters_dir, "main")
    if main_summary is None:
        print("No main cluster data; skipping dendrogram summary.")
        return

    # Copy plots to a standard cluster-XX.png naming so template can find them
    std_plots_dir = os.path.join(reports_dir, "main_clusters_standard_plots")
    os.makedirs(std_plots_dir, exist_ok=True)

    for i, (cluster_id, _) in enumerate(main_summary.iterrows(), 1):
        for candidate in [f"main_cluster-{cluster_id:02d}.png",
                          f"cluster-{cluster_id:02d}.png"]:
            src = main_plots_dir and os.path.join(main_plots_dir, candidate)
            if src and os.path.exists(src):
                shutil.copy2(src, os.path.join(std_plots_dir, f"cluster-{i:02d}.png"))
                break

    # Renumber cluster IDs 1..N for the template
    renum = main_summary.copy().reset_index()
    renum['Cluster'] = range(1, len(renum) + 1)
    renum = renum.set_index('Cluster')

    create_html_summary_figure_with_dendrogram(
        renum,
        std_plots_dir,
        os.path.join(reports_dir, "main_clusters_with_dendrogram.html"),
        dendrogram_svg_path=dendrogram_svg_path,
        acronym_mapping=acronym_mapping,
        plot_prefix="cluster",
    )


# ==============================================================================
# PIE CHART TSV (counts only)
# ==============================================================================

def create_pie_chart_tsv(cluster_summary, output_dir, cluster_type,
                         acronym_mapping=None, original_mappings=None,
                         cluster_prefix="eSS"):
    """
    Save a compact TSV for pie chart recreation.
    Format: <cluster_id>,<name>:<count>,...,total:<n>
    Only the counts file is produced (no-counts and original-mapping
    files were not needed downstream and have been removed).
    """
    print(f"Creating pie chart TSV for {cluster_type} clusters...")
    rows = []
    for i, (cluster_id, row) in enumerate(cluster_summary.iterrows()):
        sample_lines = [l.strip() for l in str(row['Samples']).split('\n') if l.strip()]
        collapsed    = _collapse_and_map_prefixes(sample_lines, acronym_mapping)
        counts       = Counter(collapsed)

        label   = f"{cluster_prefix}{cluster_id}"
        entries = ",".join(f"{p}:{c}" for p, c in counts.items())
        total   = len(sample_lines)
        if original_mappings:
            orig_names = []
            for std in sample_lines:
                orig_names.extend(original_mappings.get(std, [std]))
            total = len(orig_names)
        rows.append(f"{label},{entries},total:{total}")

    path = os.path.join(output_dir,
                        f"cluster_summary_for_pie_charts_{cluster_type}_counts.txt")
    with open(path, 'w') as f:
        f.writelines(f"{r}\n" for r in rows)
    print(f"Pie chart counts → {path}")
    return path


# ==============================================================================
# ORIGINAL SAMPLE NAMES FILE
# ==============================================================================

def create_original_sample_names_file(cluster_summary, output_dir, cluster_type,
                                      original_mappings=None, cluster_prefix="eSS"):
    """Save <cluster_id>,<original_sample1>,... one line per cluster."""
    print(f"Creating original sample names file for {cluster_type} clusters...")
    lines = []
    for cluster_id, row in cluster_summary.iterrows():
        std_samples = [l.strip() for l in str(row['Samples']).split('\n') if l.strip()]
        if original_mappings:
            originals = []
            for s in std_samples:
                originals.extend(original_mappings.get(s, [s]))
        else:
            originals = std_samples
        label = f"{cluster_prefix}{cluster_id}"
        lines.append(f"{label}," + ",".join(originals))

    path = os.path.join(output_dir,
                        f"cluster_original_sample_names_{cluster_type}.txt")
    with open(path, 'w') as f:
        f.writelines(f"{l}\n" for l in lines)
    print(f"Original sample names → {path}")
    return path


# ==============================================================================
# CLUSTER COMPOSITION SUMMARY TABLE
# ==============================================================================

def create_cluster_summary_table(cluster_summary, output_dir,
                                 acronym_mapping=None):
    """
    Saves cluster_summary_table.tsv and cluster_details_breakdown.tsv.
    (detected_categories.txt removed as it was not needed downstream.)
    """
    print("Analysing cluster composition...")

    def _parse(prefix):
        parts = prefix.split('_')
        if len(parts) < 2:
            return parts[0], parts[0], None
        species = parts[0]
        # Special three-part models
        three_part = {
            'Mouse_Bone_Marrow', 'Mouse_Lymph_Node',
            'Human_Breast_organoids', 'Human_Colon_organoids',
        }
        if '_'.join(parts[:3]) in three_part and len(parts) > 3:
            model    = '_'.join(parts[:3])
            compound = '_'.join(parts[3:]) or None
        elif len(parts) >= 3:
            model    = '_'.join(parts[:2])
            compound = '_'.join(parts[2:])
        else:
            model    = '_'.join(parts[:2])
            compound = None
        return species, model, compound

    total = len(cluster_summary)
    eligible = multi_sp = multi_model = multi_cmp = 0
    details  = []

    for _, row in cluster_summary.iterrows():
        sample_lines = [l.strip() for l in str(row['Samples']).split('\n') if l.strip()]
        collapsed    = _collapse_and_map_prefixes(sample_lines, acronym_mapping)
        unique       = set(collapsed)
        if len(unique) <= 1:
            continue
        eligible += 1
        sp_set, m_set, c_set = set(), set(), set()
        for p in unique:
            sp, m, c = _parse(p)
            sp_set.add(sp); m_set.add(m)
            if c:
                c_set.add(c)
        if len(sp_set) > 1:  multi_sp    += 1
        if len(m_set)  > 1:  multi_model += 1
        if len(c_set)  > 1:  multi_cmp   += 1
        details.append({
            'unique_sample_types': len(unique),
            'species':   sorted(sp_set),
            'models':    sorted(m_set),
            'compounds': sorted(c_set),
        })

    def pct(n, d):
        return round(n / d * 100, 1) if d else 0

    summary_data = {
        'Category': ['Total Clusters', 'Eligible (>1 type)',
                     'Multi-species', 'Multi-model', 'Multi-compound'],
        'Count':    [total, eligible, multi_sp, multi_model, multi_cmp],
        'Pct of All':      [100.0, pct(eligible, total),
                            pct(multi_sp, total), pct(multi_model, total),
                            pct(multi_cmp, total)],
        'Pct of Eligible': ['-', 100.0,
                            pct(multi_sp, eligible), pct(multi_model, eligible),
                            pct(multi_cmp, eligible)],
    }
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, "cluster_summary_table.tsv"),
                      sep='\t', index=False)
    pd.DataFrame(details).to_csv(
        os.path.join(output_dir, "cluster_details_breakdown.tsv"),
        sep='\t', index=False)
    print(summary_df.to_string(index=False))
    return summary_df, pd.DataFrame(details)


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate HTML reports for eSignatures clustering output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_df_path",    required=True,
                        help="Path to main_clusters_summary.tsv.")
    parser.add_argument("--plots_dir",        required=True,
                        help="Directory containing main_cluster-XX.png files.")
    parser.add_argument("--main_output_dir",  required=True,
                        help="Typed output dir (e.g. results/SBS/).")
    parser.add_argument("--main_working_dir", required=True,
                        help="Same as main_output_dir (kept for compatibility).")
    parser.add_argument("--dendrogram_svg",
                        help="Path to all_esignature_models_dendrogram.svg.")
    parser.add_argument("--create_dendrogram_summaries", action="store_true",
                        help="Also generate main_clusters_with_dendrogram.html.")
    parser.add_argument("--abbreviation_file",
                        help="TSV with 'compound' and 'acronym' columns.")
    parser.add_argument("--mutation_type", default="SBS",
                        choices=["SBS", "DBS", "ID"],
                        help="Mutation type (used for cluster ID prefix).")
    args = parser.parse_args()

    cfg             = get_config(args.mutation_type)
    acronym_mapping = load_acronym_mapping(args.abbreviation_file)
    original_mappings = load_original_sample_mappings(
        os.path.join(args.main_output_dir, "applied_sample_name_mappings.tsv")
    )
    reports_dir = os.path.join(args.main_output_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)

    print("=" * 60)
    print(f"GENERATING REPORTS — {cfg.name}")
    print(f"Output: {reports_dir}")
    print("=" * 60)

    # 1. Interactive report (main clusters)
    print("\n--- 1. Interactive HTML report ---")
    cluster_summary = generate_html_report(
        args.input_df_path,
        args.dendrogram_svg,
        args.plots_dir,
        os.path.join(reports_dir, "cluster_report.html"),
        acronym_mapping=acronym_mapping,
    )

    # 2. Uniform grid summary (main clusters)
    print("\n--- 2. Uniform grid summary ---")
    main_clusters_dir = os.path.join(args.main_output_dir, "main_clusters")
    if os.path.exists(main_clusters_dir):
        main_summary = load_cluster_data(main_clusters_dir, "main")
        if main_summary is not None and os.path.exists(args.plots_dir):
            renum = main_summary.copy().reset_index()
            renum['Cluster'] = range(1, len(renum) + 1)
            renum = renum.set_index('Cluster')
            create_html_summary_figure_uniform(
                renum, args.plots_dir, reports_dir, acronym_mapping=acronym_mapping
            )

    # 3. Dendrogram summary (optional)
    if args.create_dendrogram_summaries:
        print("\n--- 3. Dendrogram summary ---")
        create_dendrogram_summary(
            args.main_output_dir, reports_dir,
            args.dendrogram_svg, acronym_mapping=acronym_mapping,
        )
    else:
        print("\n--- 3. Skipping dendrogram summary ---")

    # 4. Original sample names — main
    print("\n--- 4. Original sample names (main) ---")
    create_original_sample_names_file(
        cluster_summary, args.main_output_dir, "main",
        original_mappings=original_mappings,
        cluster_prefix=cfg.cluster_prefix,
    )

    # 5. Pie chart TSV — main
    print("\n--- 5. Pie chart TSV (main) ---")
    create_pie_chart_tsv(
        cluster_summary, args.main_output_dir, "main",
        acronym_mapping=acronym_mapping,
        original_mappings=original_mappings,
        cluster_prefix=cfg.cluster_prefix,
    )

    # 6. Composition summary table
    print("\n--- 6. Cluster composition table ---")
    create_cluster_summary_table(
        cluster_summary, args.main_output_dir, acronym_mapping=acronym_mapping
    )

    # 7. Interactive report (small clusters)
    print("\n--- 7. Small clusters report ---")
    small_summary_path = os.path.join(
        args.main_working_dir, "small_clusters", "small_clusters_summary.tsv"
    )
    small_plots_dir    = os.path.join(
        args.main_working_dir, "small_clusters",
        "small_clusters_hierarchical_plots"
    )
    small_summary = None
    if os.path.exists(small_summary_path):
        small_summary = generate_html_report(
            small_summary_path,
            args.dendrogram_svg,
            small_plots_dir,
            os.path.join(reports_dir, "small_clusters_cluster_report.html"),
            acronym_mapping=acronym_mapping,
        )

    # 7b. Uniform grid summary (small clusters) -- same template as main clusters
    if small_summary is not None and os.path.exists(small_plots_dir):
        print("\n--- 7b. Uniform grid summary (small clusters) ---")
        small_renum = small_summary.copy().reset_index()
        small_renum['Cluster'] = range(1, len(small_renum) + 1)
        small_renum = small_renum.set_index('Cluster')
        create_html_summary_figure_uniform(
            small_renum, small_plots_dir, reports_dir,
            acronym_mapping=acronym_mapping,
            plot_prefix="small_cluster",
            output_filename="small_clusters_summary.html",
            title="Small Clusters Summary",
        )

    if small_summary is not None:
        # 8. Original sample names — small
        print("\n--- 8. Original sample names (small) ---")
        create_original_sample_names_file(
            small_summary, args.main_output_dir, "small",
            original_mappings=original_mappings,
            cluster_prefix=cfg.cluster_prefix,
        )
        # 9. Pie chart TSV — small
        print("\n--- 9. Pie chart TSV (small) ---")
        create_pie_chart_tsv(
            small_summary, args.main_output_dir, "small",
            acronym_mapping=acronym_mapping,
            original_mappings=original_mappings,
            cluster_prefix=cfg.cluster_prefix,
        )

    print("\n" + "=" * 60)
    print("REPORTS COMPLETE")
    print(f"All files in: {reports_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()