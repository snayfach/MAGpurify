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
from magpurify import phylo, clade, conspecific, tetra, gc, coverage, contam, clean


def cli():
    parser = argparse.ArgumentParser(
        description="Identify and remove incorrectly binned contigs from metagenome-assembled genomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version="%(prog)s 2.0")
    subparsers = parser.add_subparsers()
    phylo_parser = subparsers.add_parser(
        "phylo-markers",
        help="find taxonomic discordant contigs using a database of phylogenetic marker genes.",
        description="Find taxonomic discordant contigs using a database of phylogenetic marker genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    phylo.fetch_args(phylo_parser)
    clade_parser = subparsers.add_parser(
        "clade-markers",
        help="find taxonomic discordant contigs using a database of clade-specific marker genes.",
        description="Find taxonomic discordant contigs using a database of clade-specific marker genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    clade.fetch_args(clade_parser)
    conspecific_parser = subparsers.add_parser(
        "conspecific",
        help="find contigs that fail to align to closely related genomes.",
        description="Find contigs that fail to align to closely related genomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    conspecific.fetch_args(conspecific_parser)
    tetra_parser = subparsers.add_parser(
        "tetra-freq",
        help="find contigs with outlier tetranucleotide frequency.",
        description="Find contigs with outlier tetranucleotide frequency.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    tetra.fetch_args(tetra_parser)
    gc_parser = subparsers.add_parser(
        "gc-content",
        help="find contigs with outlier GC content.",
        description="Find contigs with outlier GC content.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    gc.fetch_args(gc_parser)
    coverage_parser = subparsers.add_parser(
        "coverage",
        help="find contigs with outlier coverage profile.",
        description="Find contigs with outlier coverage profile.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    coverage.fetch_args(coverage_parser)
    contam_parser = subparsers.add_parser(
        "known-contam",
        help="find contigs that match a database of known contaminants.",
        description="Find contigs that match a database of known contaminants.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    contam.fetch_args(contam_parser)
    clean_parser = subparsers.add_parser(
        "clean-bin",
        help="remove putative contaminant contigs from bin.",
        description="Remove putative contaminant contigs from bin.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    clean.fetch_args(clean_parser)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "phylo-markers":
            phylo_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "clade-markers":
            clade_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "conspecific":
            conspecific_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "tetra-freq":
            tetra_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "gc-content":
            gc_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "coverage":
            coverage_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "known-contam":
            contam_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "clean-bin":
            clean_parser.print_help()
            sys.exit(0)
    args = vars(parser.parse_args())
    args["func"](args)
