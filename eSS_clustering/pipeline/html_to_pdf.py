#!/usr/bin/env python3
"""
html_to_pdf.py — print the cluster summary HTML reports to one-page PDFs.

    python pipeline/html_to_pdf.py --reports_dir results/min307/SBS/reports

For each summary report in --reports_dir (cluster_summary.html,
main_clusters_with_dendrogram.html, small_clusters_summary.html) it writes a
PDF and a 300 DPI PNG with the same name next to it, e.g.
reports/cluster_summary.pdf and reports/cluster_summary.png.

Each PDF is ONE page: the report's A4-landscape width, and as tall as the
report's content, so nothing is split across pages (printing from the browser
gives A4 pages, with the last rows on a second page). The page is printed by
headless Google Chrome, the same engine as printing from the browser, so it
looks the same; the empty space under the last row is trimmed.

Chrome is looked up in $CHROME_PATH, then the usual macOS and Linux
locations. If it is not found, the script prints a warning and exits without
an error, so a run never fails over the PDFs.

The PNG is the whole page as one image (pdftoppm, part of Poppler).

Other HTML files can be given with --html (paths relative to --reports_dir).

    python pipeline/html_to_pdf.py --check

checks the tools before a run (run_pipeline.sh does this first): stops with
an error if Poppler is missing (the pipeline needs it, through pdf2image, from
Step 1 on) and warns if Google Chrome is missing (only the PDFs and PNGs of
Steps 2b and 6b need it).
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

from pdf2image import convert_from_path
from pypdf import PdfReader, PdfWriter

SUMMARY_REPORTS = [
    "cluster_summary.html",
    "main_clusters_with_dendrogram.html",
    "small_clusters_summary.html",
]

MM = 72 / 25.4                 # points per mm
PAGE_WIDTH_MM = 297            # A4 landscape, as in the reports' print CSS
PAGE_MARGIN_MM = 3             # same margin as the reports' print CSS
START_HEIGHT_MM = 600          # first try; doubled until everything fits
MAX_HEIGHT_MM = 5000           # Chrome's page size limit is about 5 m
BOTTOM_MARGIN_PT = PAGE_MARGIN_MM * MM
PNG_DPI = 300
POPPLER_TOOLS = ("pdftoppm", "pdftocairo")

POPPLER_INSTALL = ("Install Poppler: brew install poppler (macOS), "
                   "sudo apt-get install poppler-utils (Linux), "
                   "or conda install -c conda-forge poppler (in the conda environment).")
CHROME_INSTALL = ("Install Google Chrome (https://www.google.com/chrome/) "
                  "or set CHROME_PATH to its executable.")


def find_chrome():
    candidates = [
        os.environ.get("CHROME_PATH"),
        "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
        "/Applications/Chromium.app/Contents/MacOS/Chromium",
        shutil.which("google-chrome"),
        shutil.which("google-chrome-stable"),
        shutil.which("chromium"),
        shutil.which("chromium-browser"),
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None


def print_page(chrome, html_path, pdf_path, height_mm, profile_dir):
    """Print html_path to pdf_path on one page of PAGE_WIDTH_MM x height_mm."""
    # A copy of the report, in the same folder so its relative image paths
    # still work, with a later @page rule that replaces the A4 page size.
    folder = os.path.dirname(html_path)
    with open(html_path, encoding="utf-8") as f:
        html = f.read()
    override = (f"<style>@media print {{ @page {{ size: {PAGE_WIDTH_MM}mm "
                f"{height_mm}mm; margin: {PAGE_MARGIN_MM}mm; }} }}</style>")
    html = html.replace("</head>", override + "</head>", 1) if "</head>" in html \
        else override + html
    fd, tmp_html = tempfile.mkstemp(suffix=".html", prefix=".print_", dir=folder)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(html)
        result = subprocess.run(
            [chrome, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
             f"--user-data-dir={profile_dir}", "--virtual-time-budget=10000",
             f"--print-to-pdf={pdf_path}", "file://" + os.path.abspath(tmp_html)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=300)
    finally:
        os.remove(tmp_html)
    if result.returncode != 0 or not os.path.exists(pdf_path):
        raise RuntimeError(f"Chrome could not print {html_path}:\n{result.stderr[-2000:]}")
    return len(PdfReader(pdf_path).pages)


def trim_bottom(pdf_path):
    """Cut the page off just under the lowest drawn content."""
    image = convert_from_path(pdf_path, dpi=50, first_page=1, last_page=1)[0].convert("L")
    width, height = image.size
    pixels = image.load()
    last_row = 0
    for y in range(height - 1, -1, -1):
        if any(pixels[x, y] < 245 for x in range(0, width, 2)):
            last_row = y
            break
    reader = PdfReader(pdf_path)
    page = reader.pages[0]
    page_height = float(page.mediabox.height)
    content_bottom = page_height * (1 - (last_row + 1) / height)   # PDF y runs upward
    new_bottom = max(0.0, content_bottom - BOTTOM_MARGIN_PT)
    for box in (page.mediabox, page.cropbox):
        box.lower_left = (float(box.left), new_bottom)
    writer = PdfWriter()
    writer.add_page(page)
    with open(pdf_path, "wb") as f:
        writer.write(f)
    return (page_height - new_bottom) / MM


def pdf_to_png(pdf_path, dpi=PNG_DPI):
    """Save the PDF's page as <name>.png at dpi, next to the PDF."""
    if shutil.which("pdftoppm") is None:
        print(f"  Warning: pdftoppm (Poppler) not found; no PNG for "
              f"{os.path.basename(pdf_path)}. {POPPLER_INSTALL}")
        return None
    base = os.path.splitext(pdf_path)[0]
    subprocess.run(["pdftoppm", "-png", "-r", str(dpi), "-singlefile", pdf_path, base],
                   check=True)
    print(f"  {os.path.basename(base)}.png: {dpi} DPI")
    return base + ".png"


def check_tools():
    """Print what the PDF/PNG steps will use; return False if Poppler is missing."""
    ok = True
    missing = [t for t in POPPLER_TOOLS if shutil.which(t) is None]
    if missing:
        print(f"ERROR: Poppler not found ({', '.join(missing)} missing). The pipeline "
              f"needs it from Step 1 on (pdf2image) and for the PNGs. {POPPLER_INSTALL}")
        ok = False
    else:
        print(f"Poppler:       {os.path.dirname(shutil.which('pdftoppm'))}")
    chrome = find_chrome()
    if chrome is None:
        print("WARNING: Google Chrome not found. The run continues, but the one-page "
              "PDFs and PNGs of the cluster summaries (Step 2b) and the combined "
              f"species x compound figure (Step 6b) will be skipped. {CHROME_INSTALL}")
    else:
        print(f"Google Chrome: {chrome}")
    return ok


def html_to_one_page_pdf(chrome, html_path, profile_dir):
    pdf_path = os.path.splitext(html_path)[0] + ".pdf"
    height_mm = START_HEIGHT_MM
    while True:
        pages = print_page(chrome, html_path, pdf_path, height_mm, profile_dir)
        if pages == 1:
            break
        if height_mm >= MAX_HEIGHT_MM:
            raise RuntimeError(f"{html_path} does not fit on one {MAX_HEIGHT_MM} mm page")
        height_mm = min(height_mm * 2, MAX_HEIGHT_MM)
    final_height = trim_bottom(pdf_path)
    print(f"  {os.path.basename(pdf_path)}: 1 page, {PAGE_WIDTH_MM} x {final_height:.0f} mm")
    pdf_to_png(pdf_path)
    return pdf_path


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reports_dir", help="<run folder>/<TYPE>/reports")
    parser.add_argument("--check", action="store_true",
                        help="only check for Poppler and Google Chrome")
    parser.add_argument("--html", nargs="+", default=None,
                        help="HTML files to print (default: the cluster summaries)")
    args = parser.parse_args()
    if args.check:
        return 0 if check_tools() else 1
    if not args.reports_dir:
        parser.error("--reports_dir is required (or use --check)")

    chrome = find_chrome()
    if chrome is None:
        print("WARNING: Google Chrome not found; skipping the PDF and PNG versions of "
              f"the summary reports. {CHROME_INSTALL}")
        return 0

    names = args.html or [n for n in SUMMARY_REPORTS
                          if os.path.exists(os.path.join(args.reports_dir, n))]
    if not names:
        print(f"No summary reports found in {args.reports_dir}")
        return 0

    print(f"Printing {len(names)} report(s) to one-page PDFs and {PNG_DPI} DPI PNGs with {chrome}")
    with tempfile.TemporaryDirectory(prefix="chrome_profile_") as profile_dir:
        for name in names:
            html_to_one_page_pdf(chrome, os.path.join(args.reports_dir, name), profile_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
