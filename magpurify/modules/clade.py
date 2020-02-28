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
import collections
import copy
import operator
import os
import sys
from magpurify import utilities

ranks = ["k", "p", "c", "o", "f", "g", "s"]
rank_names = {
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def fetch_args(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="clade-markers")
    parser.add_argument("fna", type=str, help="Path to input genome in FASTA format")
    parser.add_argument(
        "out", type=str, help="Output directory to store results and intermediate files",
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="Number of CPUs to use",
    )
    parser.add_argument(
        "--db",
        type=str,
        help="Path to reference database. By default, the MAGPURIFY environmental variable is used",
    )
    parser.add_argument(
        "--exclude_clades",
        nargs="+",
        type=str,
        help="List of clades to exclude (ex: s__Variovorax_sp_CF313)",
    )
    parser.add_argument(
        "--min_bin_fract",
        type=float,
        default=0.6,
        help="Min fraction of bin length supported by contigs that agree with consensus taxonomy",
    )
    parser.add_argument(
        "--min_contig_fract",
        type=float,
        default=0.75,
        help="Min fraction of classified contig length that agree with consensus taxonomy",
    )
    parser.add_argument(
        "--min_gene_fract",
        type=float,
        default=0.0,
        help="Min fraction of classified genes that agree with consensus taxonomy",
    )
    parser.add_argument(
        "--min_genes",
        type=float,
        default=None,
        help="Min number of genes that agree with consensus taxonomy (default=rank-specific-cutoffs)",
    )
    parser.add_argument(
        "--lowest_rank",
        choices=["s", "g", "f", "o", "c", "p", "k"],
        help="Lowest rank for bin classification",
    )


def read_ref_taxonomy(db_dir):
    ref_taxonomy = {}
    inpath = f"{db_dir}/clade-markers/taxonomy.tsv"
    for line in open(inpath):
        ref_id, taxonomy = line.rstrip().split()
        ref_taxonomy[ref_id] = taxonomy
    return ref_taxonomy


def flag_contigs(contigs, bin):
    for bin_taxon in bin.taxonomy:
        bin_rank = bin_taxon.split("__")[0]
        for contig in contigs.values():
            contig_taxa = [
                gene.taxa[bin_rank] for gene in contig.genes if gene.taxa[bin_rank]
            ]
            # contig has no gene hits at rank
            if len(contig_taxa) == 0:
                continue
            # flag contig as discordant
            elif bin_taxon not in contig_taxa:
                contig.flagged = True
            # contig already flagged as discordant at higher rank
            elif contig.flagged:
                continue
            # flag contig as concordant
            # this is debatable...contig could be flagged as concordant,
            # but at a higher rank than the bin in annotated to
            else:
                contig.flagged = False


class Gene:
    def __init__(self):
        self.id = None
        self.aln = None
        self.taxa = dict([(rank, None) for rank in ranks])


class Contig:
    def __init__(self):
        self.id = None
        self.length = None
        self.genes = []
        self.flagged = None
        self.cons_taxa = dict([(rank, None) for rank in ranks])

    def classify(self):
        for rank in ranks:
            taxa = [g.taxa[rank] for g in self.genes if g.taxa[rank]]
            if len(taxa) > 0:
                counts = collections.Counter(taxa).items()
                self.cons_taxa[rank] = sorted(
                    counts, key=operator.itemgetter(1), reverse=True
                )[0][0]


class Bin:
    def __init__(self):
        self.cons_taxon = None
        self.bin_fract = None
        self.contig_fract = None
        self.gene_fract = None
        self.taxonomy = [None]
        self.tagged_length = 0
        self.tagged_genes = 0

    def classify(
        self,
        contigs,
        min_bin_fract,
        min_contig_fract,
        min_gene_fract,
        min_genes,
        lowest_rank,
    ):
        for rank in ranks:
            count_genes = collections.defaultdict(int)
            for contig in contigs.values():
                for gene in contig.genes:
                    if gene.taxa[rank]:
                        count_genes[gene.taxa[rank]] += 1
            count_length = collections.defaultdict(int)
            for contig in contigs.values():
                contig_taxon = contig.cons_taxa[rank]
                if contig_taxon is not None:
                    count_length[contig_taxon] += contig.length
            if sum(count_length.values()) > 0:
                cons_taxon = sorted(
                    count_length.items(), key=operator.itemgetter(1), reverse=True
                )[0][0]
                gene_fract = count_genes[cons_taxon] / float(sum(count_genes.values()))
                contig_fract = count_length[cons_taxon] / float(
                    sum(count_length.values())
                )
                bin_fract = sum(count_length.values()) / sum(
                    c.length for c in contigs.values()
                )
            else:
                cons_taxon = "NA"
                gene_fract = 0.0
                contig_fract = 0.0
                bin_fract = 0.0
            if bin_fract < min_bin_fract:
                continue
            elif contig_fract < min_contig_fract:
                continue
            elif gene_fract < min_gene_fract:
                continue
            elif count_genes[cons_taxon] < max(min_genes[rank], 1):
                continue
            else:
                self.cons_taxon = cons_taxon
                self.gene_fract = gene_fract
                self.contig_fract = contig_fract
                self.bin_fract = bin_fract
                self.tagged_genes = sum(count_genes.values())
                self.tagged_length = sum(count_length.values())
            if lowest_rank and rank == lowest_rank:
                break


def main(args):
    utilities.add_tmp_dir(args)
    utilities.check_input(args)
    utilities.check_database(args)
    print("\u001b[1m" + "• Reading database info" + "\u001b[0m")
    ref_taxonomy = read_ref_taxonomy(args["db"])
    taxon_to_taxonomy = {}
    for taxonomy in set(ref_taxonomy.values()):
        for taxon in taxonomy.split("|"):
            taxon_to_taxonomy[taxon] = taxonomy
    min_pid = {"k": 57, "p": 77, "c": 82, "o": 86, "f": 87, "g": 91, "s": 96}
    if args["min_genes"] is not None:
        args["min_genes"] = dict([(r, args["min_genes"]) for r in ranks])
    else:
        args["min_genes"] = {
            "k": 237,
            "p": 44,
            "c": 30,
            "o": 24,
            "f": 22,
            "g": 20,
            "s": 19,
        }
    print("\u001b[1m" + "\n• Calling genes with Prodigal" + "\u001b[0m")
    utilities.run_prodigal(args["fna"], args["tmp_dir"])
    print(f"  all genes: {args['tmp_dir']}/genes.[ffn|faa]")
    print(
        "\u001b[1m"
        + "\n• Performing pairwise alignment of genes against MetaPhlan2 database of clade-specific genes"
        + "\u001b[0m"
    )
    utilities.run_lastal(args["db"], args["tmp_dir"], args["threads"])
    print(f"  alignments: {args['tmp_dir']}/genes.m8")

    print("\u001b[1m" + "\n• Finding top hits to database" + "\u001b[0m")
    genes = {}
    for aln in utilities.parse_last(args["tmp_dir"] + "/genes.m8"):
        # clade exclusion
        ref_taxa = ref_taxonomy[aln["tid"]].split("|")
        if args["exclude_clades"] and any(
            [taxon in ref_taxa for taxon in args["exclude_clades"]]
        ):
            continue
        # initialize gene
        if aln["qid"] not in genes:
            genes[aln["qid"]] = Gene()
            genes[aln["qid"]].id = aln["qid"]
            genes[aln["qid"]].contig_id = aln["qid"].rsplit("_", 1)[0]

        # get top alignments
        if genes[aln["qid"]].aln is None:
            genes[aln["qid"]].aln = aln
            genes[aln["qid"]].ref_taxa = ref_taxa
        elif float(aln["score"]) > float(genes[aln["qid"]].aln["score"]):
            genes[aln["qid"]].ref_taxa = ref_taxa
    print("  %s genes with a database hit" % len(genes))
    print("\u001b[1m" + "\n• Classifying genes at each taxonomic rank" + "\u001b[0m")
    counts = {}
    for gene in genes.values():
        for ref_taxon in gene.ref_taxa:
            rank = ref_taxon.split("__")[0]
            if rank not in counts:
                counts[rank] = 0
            if rank == "t":
                continue
            elif float(gene.aln["pid"]) < min_pid[rank]:
                continue
            elif gene.aln["qcov"] < 0.4:
                continue
            elif gene.aln["tcov"] < 0.4:
                continue
            gene.taxa[rank] = ref_taxon
            counts[rank] += 1
    for rank in ranks:
        print(f"  {rank_names[rank]}: {counts[rank]} classified genes")
    print("\u001b[1m" + "\n• Taxonomically classifying contigs" + "\u001b[0m")
    contigs = {}
    for id, seq in utilities.parse_fasta(args["fna"]):
        contigs[id] = Contig()
        contigs[id].id = id
        contigs[id].length = len(seq)
    # aggregate hits by contig
    for gene in genes.values():
        contigs[gene.contig_id].genes.append(gene)
    # classify contigs at each level
    for contig in contigs.values():
        contig.classify()
    # summarize
    counts = {}
    for contig in contigs.values():
        for rank, taxon in contig.cons_taxa.items():
            if rank not in counts:
                counts[rank] = 0
            if taxon is not None:
                counts[rank] += 1
    print("  total contigs: %s" % len(contigs))
    for rank in ranks:
        print(f"  {rank_names[rank]}: {counts[rank]} classified contigs")

    print("\u001b[1m" + "\n• Taxonomically classifying genome" + "\u001b[0m")
    bin = Bin()
    bin.classify(
        contigs,
        args["min_bin_fract"],
        args["min_contig_fract"],
        args["min_gene_fract"],
        args["min_genes"],
        args["lowest_rank"],
    )
    print(f"  consensus taxon: {bin.cons_taxon}")
    print("\u001b[1m" + "\n• Identifying taxonomically discordant contigs" + "\u001b[0m")
    if bin.cons_taxon is not None:
        bin.rank_index = (
            taxon_to_taxonomy[bin.cons_taxon].split("|").index(bin.cons_taxon)
        )
        bin.taxonomy = taxon_to_taxonomy[bin.cons_taxon].split("|")[
            0 : bin.rank_index + 1
        ]
        flag_contigs(contigs, bin)
    flagged = []
    for contig in contigs.values():
        if contig.flagged:
            flagged.append(contig.id)
    out = f"{args['tmp_dir']}/flagged_contigs"
    print(f"  {len(flagged)} flagged contigs: {out}")
    with open(out, "w") as f:
        for contig in flagged:
            f.write(contig + "\n")

