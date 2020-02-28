# -*- coding: utf-8 -*-
#
#   This file is part of the magpurify package, available at:
#   https://github.com/snayfach/MAGpurify
#
#   Magpurify is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <https://www.gnu.org/licenses/>.

import argparse
import sys
import pandas as pd
from magpurify import utilities


def fetch_args(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="coverage")
    parser.add_argument("fna", type=str, help="Path to input genome in FASTA format")
    parser.add_argument(
        "out", type=str, help="Output directory to store results and intermediate files",
    )
    parser.add_argument(
        "bams", nargs="+", type=str, help="Path to input sorted BAM file(s)",
    )
    parser.add_argument(
        "--max-deviation",
        type=float,
        default=5.0,
        help="Contigs with coverage greater than [max-deviation * mean coverage] or less than [(1/max-deviation) * mean coverage] will be flagged as outliers",
    )


def main(args):
    utilities.add_tmp_dir(args)
    utilities.check_input(args)
    utilities.check_dependencies(["coverm"])
    print("\u001b[1m" + "• Computing contig coverage" + "\u001b[0m")
    utilities.run_coverm(args["bams"], args["tmp_dir"])
    coverage_df = pd.read_csv(f"{args['tmp_dir']}/coverage.tsv", sep="\t", index_col=0)
    print("\u001b[1m" + "\n• Identifying outlier contigs" + "\u001b[0m")
    mag_id_list = []
    for id, seq in utilities.parse_fasta(args["fna"]):
        mag_id_list.append(id)
    mag_coverage_df = coverage_df.loc[mag_id_list]
    largest_mean_coverage_sample = mag_coverage_df.mean(axis=0).idxmax()
    if mag_coverage_df.shape[1] > 1:
        print(
            "\u001b[1m"
            + f"\n• Sample being used for outlier detection: {largest_mean_coverage_sample.split()[0]}"
            + "\u001b[0m"
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
    print(f"  {len(flagged)} flagged contigs: {out}")
    with open(out, "w") as f:
        for contig in flagged:
            f.write(contig + "\n")

