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
import os
import sys
from operator import itemgetter
from magpurify import utilities


def fetch_args(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="clean-bin")
    parser.add_argument(
        "fna",
        type=str,
        help="Path to input genome in FASTA format"
    )
    parser.add_argument(
        "out",
        type=str,
        help="Output directory to store results and intermediate files",
    )
    parser.add_argument(
        "out_fna",
        type=str,
        help="Path to the output FASTA file",
    )


def main(args):
    utilities.check_input(args)
    print("\u001b[1m" + "• Reading genome bin" + "\u001b[0m")
    bin = {}
    for id, seq in utilities.parse_fasta(args["fna"]):
        bin[id] = seq
    bin_length = round(sum(len(_) for _ in bin.values()) / 1000, 2)
    print(f"  genome length: {len(bin)} contigs, {bin_length} Kbp")
    print("\u001b[1m" + "\n• Reading flagged contigs" + "\u001b[0m")
    flagged_contigs = []
    programs = [
        "phylo-markers",
        "clade-markers",
        "conspecific",
        "tetra-freq",
        "gc-content",
        "coverage",
        "known-contam",
    ]
    for program in programs:
        path = f"{args['out']}/{program}/flagged_contigs"
        if not os.path.exists(path):
            print(f"  {program}: no output file found")
        else:
            contigs = [_.rstrip() for _ in open(path)]
            bases = round(sum(len(bin[id]) for id in contigs) / 1000, 2)
            flagged_contigs += contigs
            print(f"  {program}: {len(contigs)} contigs, {bases} Kbp")
    flagged_contigs = list(set(flagged_contigs))
    flagged_length = round(sum(len(bin[id]) for id in flagged_contigs) / 1000, 2)
    print("\u001b[1m" + "\n• Removing flagged contigs" + "\u001b[0m")
    clean = bin.copy()
    for id in flagged_contigs:
        del clean[id]
    clean_length = round(sum(len(_) for _ in clean.values()) / 1000, 2)
    print(f"  removed: {len(flagged_contigs)} contigs, {flagged_length} Kbp")
    print(f"  remains: {len(clean)} contigs, {clean_length} Kbp")
    with open(args['out_fna'], "w") as f:
        for id, seq in clean.items():
            f.write(">" + id + "\n" + seq + "\n")
    print(f"  cleaned bin: {args['out_fna']}")
