#!/usr/bin/env python

import argparse
import sys
import pandas as pd
from . import utility


def fetch_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=argparse.SUPPRESS,
        description="MAGpurify: coverage module: find contigs with outlier coverage profile",
    )
    parser.add_argument("program", help=argparse.SUPPRESS)
    parser.add_argument("fna", type=str, help="""Path to input genome in FASTA format""")
    parser.add_argument(
        "out",
        type=str,
        help="""Output directory to store results and intermediate files""",
    )
    parser.add_argument(
        "--bams",
        nargs="+",
        type=str,
        required=True,
        help="""Paths to input sorted BAM files""",
    )
    parser.add_argument(
        "--max-deviation",
        type=float,
        default=5.0,
        help="""Contigs with coverage greater than [max-deviation * mean coverage] or less than [(1/max-deviation) * mean coverage] will be flagged as outliers""",
    )
    args = vars(parser.parse_args())
    return args


def main():
    args = fetch_args()
    utility.add_tmp_dir(args)
    utility.check_input(args)
    utility.check_dependencies(["coverm"])
    print("\n## Computing contig coverage")
    utility.run_coverm(args["bams"], args["tmp_dir"])
    coverage_df = pd.read_csv(f"{args['tmp_dir']}/coverage.tsv", sep="\t", index_col=0)
    print("\n## Identifying outlier contigs")
    mag_id_list = []
    for id, seq in utility.parse_fasta(args["fna"]):
        mag_id_list.append(id)
    mag_coverage_df = coverage_df.loc[mag_id_list]
    largest_mean_coverage_sample = mag_coverage_df.mean(axis=0).idxmax()
    if mag_coverage_df.shape[1] > 1:
        print(
            f"\n## Sample being used for outlier detection: {largest_mean_coverage_sample.split()[0]}"
        )
    mag_coverage_df = mag_coverage_df.loc[:, largest_mean_coverage_sample]
    if mag_coverage_df.mean() < 1:
        sys.exit(
            "\nError: The average coverage is less than 1 in all the supplied BAM files"
        )
    outliers = ((mag_coverage_df / mag_coverage_df.mean()) >= args["max_deviation"]) | (
        (mag_coverage_df / mag_coverage_df.mean()) <= 1 / args["max_deviation"]
    )
    flagged = mag_coverage_df.loc[outliers].index.tolist()
    out = f"{args['tmp_dir']}/flagged_contigs"
    print(f"   {len(flagged)} flagged contigs: {out}")
    with open(out, "w") as f:
        for contig in flagged:
            f.write(contig + "\n")

