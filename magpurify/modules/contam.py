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
from magpurify import utilities


def fetch_args(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="known-contam")
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
        "--threads",
        type=int,
        default=1,
        help="Number of CPUs to use",
    )
    parser.add_argument(
        "--db",
        type=str,
        help="Path to reference database. By default, the IMAGEN_DB environmental variable is used",
    )
    parser.add_argument(
        "--pid",
        type=float,
        default=98,
        help="Minimum %% identity to reference",
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5, help="Maximum evalue"
    )
    parser.add_argument(
        "--qcov",
        type=float,
        default=25,
        help="Minimum percent query coverage",
    )


def run_blastn(query, db, out, threads, qcov=25, pid=98, evalue=1e-5):
    cmd = "blastn "
    cmd += f"-query {query} "
    cmd += f"-db {db} "
    cmd += f"-out {out} "
    cmd += "-outfmt '6 std qlen slen' "
    cmd += "-max_target_seqs 1 "
    cmd += "-max_hsps 1 "
    cmd += f"-qcov_hsp_perc {qcov} "
    cmd += f"-perc_identity {pid} "
    cmd += f"-evalue {evalue} "
    cmd += f"-num_threads {threads} "
    utilities.run_process(cmd)


def main(args):
    utilities.add_tmp_dir(args)
    utilities.check_input(args)
    utilities.check_dependencies(["blastn"])
    utilities.check_database(args)
    utilities.add_tmp_dir(args)
    print("\u001b[1m" + "• Searching database with BLASTN" + "\u001b[0m")
    for target in ["hg38", "phix"]:
        db = f"{args['db']}/known-contam/{target}/{target}"
        out = f"{args['tmp_dir']}/{target}.m8"
        run_blastn(
            args["fna"],
            db,
            out,
            args["threads"],
            args["qcov"],
            args["pid"],
            args["evalue"],
        )
    print("\u001b[1m" + "\n• Identifying contigs with hits to database" + "\u001b[0m")
    flagged = set([])
    for target in ["hg38", "phix"]:
        out = f"{args['tmp_dir']}/{target}.m8"
        for r in utilities.parse_blast(out):
            flagged.add(r["qname"])
    flagged = list(flagged)
    out = f"{args['tmp_dir']}/flagged_contigs"
    print(f"  {len(flagged)} flagged contigs: {out}")
    with open(out, "w") as f:
        for contig in flagged:
            f.write(contig + "\n")
