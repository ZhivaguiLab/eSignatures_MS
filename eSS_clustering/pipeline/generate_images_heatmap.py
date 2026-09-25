import os
import sys
import argparse
import pandas as pd
import sigProfilerPlotting as sigPlt
from pdf2image import convert_from_path
import shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.mutation_type import get_config, call_plot_function


def process_esignatures(pdf_path, output_dir, cfg):
    """
    Converts each page of a single multi-page PDF into a separate PNG.
    Images are named using the mutation-type cluster prefix (eSS, eDS, eIS).
    """
    print("--- Starting eSignatures PDF to PNG Conversion ---")
    if not os.path.exists(pdf_path):
        print("❌ ERROR: eSignatures PDF not found at '{}'. Skipping.".format(pdf_path))
        return

    try:
        print("Reading PDF: {}".format(pdf_path))
        images = convert_from_path(pdf_path, dpi=300)

        for i, image in enumerate(images):
            output_png_path = os.path.join(output_dir, '{}{}.png'.format(cfg.cluster_prefix, i + 1))
            image.save(output_png_path, 'PNG')

        print("✅ Successfully converted {} pages into PNGs in '{}'".format(len(images), output_dir))

    except Exception as e:
        print("❌ An error occurred during eSignatures processing: {}".format(e))


def process_cosmic(txt_path, output_dir, cfg):
    """
    Generates a PNG plot for each column in a COSMIC data file.
    Uses the appropriate sigProfilerPlotting function for the mutation type.
    """
    print("\n--- Starting COSMIC Data to PNG Conversion ({}) ---".format(cfg.name))
    if not os.path.exists(txt_path):
        print("❌ ERROR: COSMIC data file not found at '{}'. Skipping.".format(txt_path))
        return

    temp_dir     = os.path.join(output_dir, "temp_files")
    temp_tsv_dir = os.path.join(temp_dir, "tsvs")
    temp_pdf_dir = os.path.join(temp_dir, "pdfs")
    os.makedirs(temp_tsv_dir, exist_ok=True)
    os.makedirs(temp_pdf_dir, exist_ok=True)

    print("Reading COSMIC data from: {}".format(txt_path))
    try:
        df = pd.read_csv(txt_path, sep='\t', index_col=0, header=0)
        total_samples = len(df.columns)
        print("Found {} signatures to process.".format(total_samples))

        for i, sample in enumerate(df.columns):
            print("Processing {}/{}: {}...".format(i + 1, total_samples, sample))

            # 1. Write single-column TSV
            sample_df   = df.filter([sample])
            matrix_path = os.path.join(temp_tsv_dir, "{}.tsv".format(sample))
            sample_df.to_csv(matrix_path, index=True, header=True, sep='\t')

            # 2. Plot using the correct function for this mutation type
            call_plot_function(cfg, matrix_path, temp_pdf_dir, sample)

            # 3. Convert the resulting PDF to PNG
            # PDF name follows sigProfilerPlotting convention: <prefix>_<project>.pdf
            pdf_path = os.path.join(temp_pdf_dir,
                                    "{}_{}.pdf".format(cfg.pdf_prefix, sample))
            if os.path.exists(pdf_path):
                images = convert_from_path(pdf_path, dpi=300)
                if images:
                    images[0].save(
                        os.path.join(output_dir, "{}.png".format(sample)), 'PNG'
                    )
            else:
                print("  ⚠️ Warning: PDF for sample {} was not generated.".format(sample))

        shutil.rmtree(temp_dir)
        print("✅ Successfully processed {} samples → '{}'".format(
            total_samples, output_dir))
        print("   Temporary files removed.")

    except Exception as e:
        print("❌ An error occurred during COSMIC processing: {}".format(e))


def main():
    parser = argparse.ArgumentParser(
        description="Convert PDF plots from eSignatures and COSMIC data into PNG images.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--esignature_pdf",
                        help="(Optional) Path to the multi-page eSignature PDF file.")
    parser.add_argument("--cosmic_txt",
                        help="(Optional) Path to the COSMIC data file.")
    parser.add_argument("--output_dir", required=True,
                        help="The main output directory where results will be saved.")
    parser.add_argument("--mutation_type", default="SBS",
                        choices=["SBS", "DBS", "ID"],
                        help="Mutation type — controls which sigProfilerPlotting "
                             "function is called and which PDF name pattern is used.")
    args = parser.parse_args()

    cfg = get_config(args.mutation_type)

    os.makedirs(args.output_dir, exist_ok=True)
    cluster_plots_dir = os.path.join(args.output_dir, "cluster_plots")
    cosmic_plots_dir  = os.path.join(args.output_dir, "cosmic_plots")
    os.makedirs(cluster_plots_dir, exist_ok=True)
    os.makedirs(cosmic_plots_dir,  exist_ok=True)

    if args.esignature_pdf:
        process_esignatures(args.esignature_pdf, cluster_plots_dir, cfg)

    if args.cosmic_txt:
        process_cosmic(args.cosmic_txt, cosmic_plots_dir, cfg)

    if not args.esignature_pdf and not args.cosmic_txt:
        print("No input files provided. "
              "Please specify --esignature_pdf and/or --cosmic_txt.")

    print("\nAll tasks complete.")


if __name__ == "__main__":
    main()