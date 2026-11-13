#!/usr/bin/env python3
"""
combined_matrix_figure.py — one figure from the three species × compound
matrix PDFs (Step 6), saved as HTML and as a one-page PDF.

    python pipeline/combined_matrix_figure.py --run_dir results/min307/SBS

Reads, from --run_dir (written by species_compound_matrix.py, Step 6):

    species_compound_matrix_human.pdf
    species_compound_matrix_mouse.pdf
    species_compound_matrix_other_species.pdf

and writes, in the same folder:

    species_compound_matrix_combined.html
    species_compound_matrix_combined.pdf     (one page; needs Google Chrome)

The panels are the Step 6 PDFs themselves, converted to vector SVG
(pdftocairo, part of Poppler), so the figure always shows exactly what Step 6
drew. Each PDF's own legend is cropped off and one shared legend is added.
The three panels are scaled so their grid cells are the same size.

Panel titles: a Human / b Mouse / c Other experimental models, each with
"<models> models · <exposures> exposures · <profiles> profiles", counted from
the same cluster outputs Step 6 uses (species_compound_matrix.py).
"""

import argparse
import base64
import html
import os
import re
import subprocess
import sys
import tempfile

from pdf2image import convert_from_path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import species_compound_matrix as scm  # noqa: E402
import html_to_pdf  # noqa: E402
from utils.naming import load_acronym_mapping  # noqa: E402

REPO_ROOT = os.path.dirname(HERE)

PANELS = [  # (group in species_compound_matrix.py, PDF suffix, title, rule colour)
    ("Human", "human", "Human experimental models", "#7cbf55"),
    ("Mouse", "mouse", "Mouse experimental models", "#6ec9d8"),
    ("Other", "other_species", "Other experimental models", "#f4b63c"),
]

FIGURE_WIDTH_PX = 1100     # widest panel; ~291 mm, the A4-landscape print width
MEASURE_DPI = 40           # resolution used only to find the legend and grid
WHITE = 253                # a row is blank if every pixel is at least this light


# ---------------------------------------------------------------------------
# Counts for the panel titles (same data and grouping as Step 6)
# ---------------------------------------------------------------------------

def panel_counts(run_dir, abbreviation_file, compound_grouping):
    acronyms = load_acronym_mapping(abbreviation_file) if abbreviation_file else {}
    grouping = scm.load_compound_grouping_config(compound_grouping)
    data = scm.collect_species_compound_data(run_dir, acronyms, grouping)
    counts = {g: {"models": set(), "compounds": set(), "profiles": 0} for g, *_ in PANELS}
    for (sm, compound), d in data.items():
        if not d["total"]:
            continue
        sp = sm.split("_")[0]
        group = sp if sp in ("Human", "Mouse") and sm != "C.elegans" else "Other"
        counts[group]["models"].add(sm)
        counts[group]["compounds"].add(compound)
        counts[group]["profiles"] += d["total"]
    return {g: {"models": len(c["models"]), "compounds": len(c["compounds"]),
                "profiles": c["profiles"]} for g, c in counts.items()}


# ---------------------------------------------------------------------------
# Reading the Step 6 PDFs
# ---------------------------------------------------------------------------

def measure(pdf_path):
    """Return (crop_fraction, grid_width_fraction).

    crop_fraction: fraction of the page height to keep, cutting in the blank
    gap between the matrix and the legend box under it.
    grid_width_fraction: width of the grid's bottom line as a fraction of the
    page width (divided by the number of columns, this is the cell size).
    """
    im = convert_from_path(pdf_path, dpi=MEASURE_DPI)[0].convert("L")
    w, h = im.size
    px = im.load()

    def blank(y):
        return all(px[x, y] >= WHITE for x in range(w))

    y = h - 1
    while y > 0 and blank(y):            # white margin under the legend
        y -= 1
    while y > 0 and not blank(y):        # the legend box
        y -= 1
    legend_gap_bottom = y
    while y > 0 and blank(y):            # gap between grid and legend
        y -= 1
    grid_bottom = y
    if grid_bottom <= 0 or legend_gap_bottom - grid_bottom < 2:
        raise RuntimeError(f"Could not find the legend under the matrix in {pdf_path}")

    # Longest run of grey pixels in the grid's bottom rows = grid width.
    best = 0
    for yy in range(max(0, grid_bottom - 2), grid_bottom + 1):
        run = 0
        for x in range(w):
            if px[x, yy] < 240:
                run += 1
                best = max(best, run)
            else:
                run = 0
    cut = (grid_bottom + legend_gap_bottom) / 2 + 0.5
    return cut / h, best / w


def pdf_to_svg(pdf_path, keep_fraction):
    """Vector SVG of the PDF page, cropped to the top keep_fraction."""
    with tempfile.TemporaryDirectory() as tmp:
        svg_path = os.path.join(tmp, "page.svg")
        subprocess.run(["pdftocairo", "-svg", pdf_path, svg_path], check=True)
        with open(svg_path, encoding="utf-8") as f:
            svg = f.read()
    m = re.search(r'<svg[^>]*\bviewBox="([\d.\s-]+)"', svg)
    x0, y0, vw, vh = (float(v) for v in m.group(1).split())
    new_h = vh * keep_fraction
    root = m.group(0)
    new_root = (root.replace(m.group(1), f"{x0} {y0} {vw} {new_h:.3f}")
                    .replace(re.search(r'height="[^"]*"', root).group(0), f'height="{new_h:.3f}pt"'))
    return svg.replace(root, new_root, 1), vw, new_h


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------

CSS = """
@page { size: A4 landscape; margin: 3mm; }
* { box-sizing: border-box; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
body { margin: 0; font-family: Arial, Helvetica, sans-serif; color: #222; background: #fff; }
.fig { width: %(width)dpx; padding: 4px 6px; }
.panel { margin-bottom: 14px; }
.ptitle { display: flex; align-items: baseline; gap: 10px; font-size: 15px; font-weight: bold;
          padding: 0 0 3px 2px; border-bottom: 3px solid var(--accent); margin-bottom: 6px; }
.ptitle .letter { font-size: 20px; }
.ptitle small { font-weight: normal; color: #666; font-size: 11px; }
.row { display: flex; align-items: flex-end; gap: 48px; }
.row img { display: block; }
.legend { display: flex; flex-direction: column; gap: 6px; font-size: 12px; margin-bottom: 14px; }
.legend .li { display: flex; align-items: center; gap: 10px; }
.legend .dot { width: 10px; height: 10px; border-radius: 50%%; display: inline-block; }
"""


def legend_html():
    c = scm.CLUSTER_OUTCOME_COLORS
    items = [("Main Cluster", c["main_cluster"]), ("Small Cluster (n=2)", c["small_cluster"]),
             ("Singleton", c["singleton"])]
    return ('<div class="legend">' + "".join(
        f'<div class="li"><span class="dot" style="background:{col}"></span>{lab}</div>'
        for lab, col in items) + "</div>")


def build(run_dir, abbreviation_file, compound_grouping):
    counts = panel_counts(run_dir, abbreviation_file, compound_grouping)

    panels = []
    for group, suffix, title, colour in PANELS:
        pdf = os.path.join(run_dir, f"species_compound_matrix_{suffix}.pdf")
        if not os.path.exists(pdf) or counts[group]["models"] == 0:
            print(f"  {group}: no matrix, skipped")
            continue
        keep, grid_frac = measure(pdf)
        svg, vw, vh = pdf_to_svg(pdf, keep)
        cell_pt = grid_frac * vw / counts[group]["compounds"]
        panels.append(dict(group=group, title=title, colour=colour, svg=svg,
                           vw=vw, vh=vh, cell_pt=cell_pt))

    # One scale (px per cell) for all panels; the widest panel fills the width.
    cell_px = min(FIGURE_WIDTH_PX * p["cell_pt"] / p["vw"] for p in panels)

    body = []
    for i, p in enumerate(panels):
        scale = cell_px / p["cell_pt"]
        width, height = p["vw"] * scale, p["vh"] * scale
        c = counts[p["group"]]
        data_uri = "data:image/svg+xml;base64," + base64.b64encode(p["svg"].encode()).decode()
        img = (f'<img src="{data_uri}" width="{width:.1f}" height="{height:.1f}" '
               f'alt="{html.escape(p["title"])}">')
        extra = legend_html() if i == len(panels) - 1 else ""
        body.append(
            f'<div class="panel" style="--accent:{p["colour"]}"><div class="ptitle">'
            f'<span class="letter">{"abcdefg"[i]}</span>{html.escape(p["title"])}'
            f'<small>{c["models"]} models · {c["compounds"]} exposures · '
            f'{c["profiles"]} profiles</small></div>'
            f'<div class="row">{img}{extra}</div></div>')
        print(f"  {p['group']}: {c['models']} models, {c['compounds']} exposures, "
              f"{c['profiles']} profiles; drawn at {width:.0f} x {height:.0f} px")

    return ("<!DOCTYPE html><html><head><meta charset='utf-8'>"
            "<title>Species compound matrix</title>"
            f"<style>{CSS % {'width': FIGURE_WIDTH_PX + 12}}</style></head>"
            f"<body><div class='fig'>{''.join(body)}</div></body></html>")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run_dir", required=True, help="<run folder>/<TYPE>, e.g. results/min307/SBS")
    parser.add_argument("--abbreviation_file",
                        default=os.path.join(REPO_ROOT, "config", "abv_table_clusters.txt"))
    parser.add_argument("--compound_grouping",
                        default=os.path.join(REPO_ROOT, "config", "compound_grouping.yaml"))
    parser.add_argument("--no_pdf", action="store_true", help="write the HTML only")
    args = parser.parse_args()

    print("Combined species × compound figure")
    page = build(args.run_dir, args.abbreviation_file, args.compound_grouping)
    out_html = os.path.join(args.run_dir, "species_compound_matrix_combined.html")
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(page)
    print(f"  → {out_html}")

    if args.no_pdf:
        return 0
    chrome = html_to_pdf.find_chrome()
    if chrome is None:
        print("Warning: Google Chrome not found (set CHROME_PATH); skipping the PDF.")
        return 0
    with tempfile.TemporaryDirectory(prefix="chrome_profile_") as profile_dir:
        html_to_pdf.html_to_one_page_pdf(chrome, out_html, profile_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
