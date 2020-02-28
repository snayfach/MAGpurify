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
    parser.set_defaults(program="conspecific")
    parser.add_argument("fna", type=str, help="Path to input genome in FASTA format")
    parser.add_argument(
        "out", type=str, help="Output directory to store results and intermediate files",
    )
    parser.add_argument(
        "mash_sketch", type=str, help="Path to Mash sketch of reference genomes",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of CPUs to use",
    )
    parser.add_argument(
        "--mash-dist",
        type=float,
        default=0.05,
        help="Mash distance to reference genomes",
    )
    parser.add_argument(
        "--max-genomes", type=int, default=25, help="Max number of genomes to use",
    )
    parser.add_argument(
        "--min-genomes", type=int, default=1, help="Min number of genomes to use",
    )
    parser.add_argument(
        "--contig-aln",
        type=float,
        default=0.50,
        help="Minimum fraction of contig aligned to reference",
    )
    parser.add_argument(
        "--contig-pid",
        type=float,
        default=95.0,
        help="Minimum percent identity of contig aligned to reference",
    )
    parser.add_argument(
        "--hit-rate", type=float, default=0.00, help="Hit rate for flagging contigs",
    )
    parser.add_argument(
        "--exclude", nargs="+", default="", help="List of references to exclude",
    )


def run_mash(mash_sketch, fna_path, tmp_dir, threads=1):
    out_path = f"{tmp_dir}/mash.dist"
    command = f"mash dist -p {threads} -d 0.25 {fna_path} {mash_sketch} > {out_path}"
    out, err = utilities.run_process(command)
    with open(tmp_dir + "/id_map.tsv", "w") as f:
        for id, rec in enumerate(utilities.parse_mash(out_path)):
            f.write(str(id) + "\t" + rec["target"] + "\n")


def find_conspecific(tmp_dir, max_dist, exclude):
    targets = []
    for rec in utilities.parse_mash(f"{tmp_dir}/mash.dist"):
        if rec["query"] == rec["target"]:
            continue
        elif rec["dist"] > max_dist:
            continue
        elif rec["pvalue"] > 1e-3:
            continue
        elif rec["target"] in exclude:
            continue
        else:
            targets.append([rec["target"], rec["dist"]])
    targets = sorted(targets, key=itemgetter(1))
    return targets


def blastn(query, target, outdir, id):
    out_path = outdir + "/" + id + ".m8"
    if not os.path.exists(out_path):
        cmd = "blastn -outfmt '6 std qlen slen' "
        cmd += "-max_target_seqs 1 -max_hsps 1 "
        cmd += f"-query {query} -subject {target} "
        out, err = utilities.run_process(cmd)
        out = out.decode("utf-8")
        open(out_path, "w").write(out)
    else:
        out = open(out_path).read()
    return out


def id_blast_hits(blast_out, min_aln, min_pid):
    blast_hits = set([])
    for rec in utilities.parse_blast(blast_out, type="string"):
        if rec["qcov"] < min_aln:
            continue
        elif rec["pid"] < min_pid:
            continue
        else:
            blast_hits.add(rec["qname"])
    return blast_hits


def align_contigs(args, genomes):
    target_to_id = {}
    for line in open(args["tmp_dir"] + "/id_map.tsv"):
        id, target = line.split()
        target_to_id[target] = id
    alignments = []
    for genome_path, mash_dist in genomes:
        blast_out = blastn(
            args["fna"], genome_path, args["tmp_dir"], target_to_id[target]
        )
        alignments.append(blast_out)
    return alignments


def find_contig_targets(args, genomes, alignments):
    contigs = dict(
        [
            (id, {"hits": 0, "len": len(seq), "genomes": []})
            for id, seq in utilities.parse_fasta(args["fna"])
        ]
    )
    for genome, alns in zip(genomes, alignments):
        hits = id_blast_hits(alns, args["contig_aln"], args["contig_pid"])
        for contig in hits:
            contigs[contig]["hits"] += 1
            contigs[contig]["genomes"].append(genome[0])
    for id in contigs:
        hit_rate = contigs[id]["hits"] / float(len(alignments))
        contigs[id]["hit_rate"] = hit_rate
    return contigs


def flag_contigs(args, contigs):
    flagged = []
    for id in contigs:
        if contigs[id]["hit_rate"] > args["hit_rate"]:
            continue
        else:
            flagged.append(id)
    return flagged


def main(args):
    utilities.add_tmp_dir(args)
    utilities.check_input(args)
    utilities.check_dependencies(["mash"])
    if not os.path.exists(args["mash_sketch"]):
        sys.exit(f"\nError: mash sketch '{args['mash_sketch']}' not found\n")
    print("\u001b[1m" + "• Finding conspecific genomes in database" + "\u001b[0m")
    run_mash(args["mash_sketch"], args["fna"], args["tmp_dir"], args["threads"])
    genomes = find_conspecific(args["tmp_dir"], args["mash_dist"], args["exclude"])
    print(f"  {len(genomes)} genomes within {args['mash_dist']} mash-dist")
    out = f"{args['tmp_dir']}/conspecific.list"
    with open(out, "w") as f:
        f.write("genome_id\tmash_dist\n")
        for genome_id, mash_dist in genomes:
            f.write(genome_id + "\t" + str(mash_dist) + "\n")
    print(f"  list of genomes: {out}")
    print(f"  mash output: {args['tmp_dir']}/mash.dist")
    if len(genomes) < args["min_genomes"]:
        sys.exit("\nError: insufficient number of conspecific genomes\n")
    if len(genomes) > args["max_genomes"]:
        print(
            "\u001b[1m"
            + f"\n• Selecting top {args['max_genomes']} most-similar genomes"
            + "\u001b[0m"
        )
        genomes = genomes[0 : args["max_genomes"]]
        out = f"{args['tmp_dir']}/conspecific_subset.list"
        with open(out, "w") as f:
            f.write("genome_id\tmash_dist\n")
            for genome_id, mash_dist in genomes:
                f.write(genome_id + "\t" + str(mash_dist) + "\n")
        print(f"  list of genomes: {out}")
    print(
        "\u001b[1m"
        + "\n• Performing pairwise alignment of contigs in bin to database genomes"
        + "\u001b[0m"
    )
    alignments = align_contigs(args, genomes)
    num_alns = sum(len(_.split("\n")) for _ in alignments)
    print(f"  total alignments: {num_alns}")
    print("\u001b[1m" + "\n• Summarizing alignments" + "\u001b[0m")
    contigs = find_contig_targets(args, genomes, alignments)
    out = f"{args['tmp_dir']}/contig_hits.tsv"
    with open(out, "w") as f:
        f.write("contig_id\tlength\talignment_rate\n")
        for contig, values in contigs.items():
            row = [contig, str(values["len"]), f"{values['hits']}/{len(genomes)}"]
            f.write("\t".join(row) + "\n")
    print(f"  contig features: {out}")
    print(
        "\u001b[1m"
        + "\n• Identifying contigs with no conspecific alignments"
        + "\u001b[0m"
    )
    flagged = flag_contigs(args, contigs)
    out = f"{args['tmp_dir']}/flagged_contigs"
    with open(out, "w") as f:
        for contig in flagged:
            f.write(contig + "\n")
    print(f"  {len(flagged)} flagged contigs: {out}")
