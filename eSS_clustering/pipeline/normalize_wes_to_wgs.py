#!/usr/bin/env python3
"""
normalize_wes_to_wgs.py — standalone wes-to-wgs opportunity normalization.

Method 1 (liftover, manuscript revision item A4): human and mouse WES samples
are rescaled onto their own genome's WGS trinucleotide opportunity (GRCh38,
mm10), keeping each sample's total. WGS samples and rat, chicken and
C. elegans samples are unchanged. See NORMALIZATION_APPROACH.md.

The same normalization runs inside the pipeline with
`bash run_pipeline.sh SBS 0.9 0.85 wes-to-wgs`; both use
pipeline/utils/opportunity_normalization.py. This script writes a normalized
copy of a preprocessed folder, e.g. to inspect or reuse the profiles.

Usage
-----
    python pipeline/normalize_wes_to_wgs.py \\
        --input_dir data/input_cleaned/SBS \\
        --output_dir data/normalized_wes-to-wgs/SBS

The input is the preprocessed data (after the 307 filter and exclusions),
made by any pipeline run. To cluster the output directly:
    python pipeline/perform_clustering.py --mutation_type SBS \\
        --data_dir data/normalized_wes-to-wgs --skip_preprocessing \\
        --output_dir results/wes_to_wgs_standalone
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils.opportunity_normalization import DEFAULT_CONTEXT_DIR, normalize_directory  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input_dir", default="data/input_cleaned/SBS",
                        help="Preprocessed SBS folder (<species>_SBS96.txt count files)")
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--context_dist_dir", default=DEFAULT_CONTEXT_DIR,
                        help="Folder with context_counts_<genome>_96[_exome].csv")
    args = parser.parse_args()
    normalize_directory("wes-to-wgs", args.input_dir, args.output_dir, args.context_dist_dir)
    print(f"\nDone. Output -> {args.output_dir}")


if __name__ == "__main__":
    main()
