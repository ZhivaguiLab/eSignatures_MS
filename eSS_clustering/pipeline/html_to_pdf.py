#!/usr/bin/env python3
"""
html_to_pdf.py — print the cluster summary HTML reports to one-page PDFs.

    python pipeline/html_to_pdf.py --reports_dir results/min307/SBS/reports

For each summary report in --reports_dir (cluster_summary.html,
main_clusters_with_dendrogram.html, small_clusters_summary.html) it writes a
PDF with the same name next to it, e.g. reports/cluster_summary.pdf.

Each PDF is ONE page: the report's A4-landscape width, and as tall as the
report's content, so nothing is split across pages (printing from the browser
gives A4 pages, with the last rows on a second page). The page is printed by
headless Google Chrome, the same engine as printing from the browser, so it
looks the same; the empty space under the last row is trimmed.

Chrome is looked up in $CHROME_PATH, then the usual macOS and Linux
locations. If it is not found, the script prints a warning and exits without
an error, so a run never fails over the PDFs.

Other HTML files can be given with --html (paths relative to --reports_dir).
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
    return pdf_path


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reports_dir", required=True, help="<run folder>/<TYPE>/reports")
    parser.add_argument("--html", nargs="+", default=None,
                        help="HTML files to print (default: the cluster summaries)")
    args = parser.parse_args()

    chrome = find_chrome()
    if chrome is None:
        print("Warning: Google Chrome not found (set CHROME_PATH); "
              "skipping the PDF versions of the summary reports.")
        return 0

    names = args.html or [n for n in SUMMARY_REPORTS
                          if os.path.exists(os.path.join(args.reports_dir, n))]
    if not names:
        print(f"No summary reports found in {args.reports_dir}")
        return 0

    print(f"Printing {len(names)} report(s) to one-page PDFs with {chrome}")
    with tempfile.TemporaryDirectory(prefix="chrome_profile_") as profile_dir:
        for name in names:
            html_to_one_page_pdf(chrome, os.path.join(args.reports_dir, name), profile_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
