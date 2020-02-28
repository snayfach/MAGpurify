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
import csv
import os
import sys
from collections import Counter
from magpurify import utilities


def fetch_args(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="phylo-markers")
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
        help="Path to reference database. By default, the MAGPURIFYDB environmental variable is used",
    )
    parser.add_argument(
        "--continue",
        action="store_true",
        default=False,
        help="Go straight to quality estimation and skip all previous steps",
    )
    parser.add_argument(
        "--max_target_seqs",
        type=int,
        default=1,
        help="Maximum number of targets reported in BLAST table",
    )
    parser.add_argument(
        "--cutoff_type",
        choices=["strict", "sensitive", "none"],
        default="strict",
        help="Use strict or sensitive %%ID cutoff for taxonomically annotating genes",
    )
    parser.add_argument(
        "--seq_type",
        choices=["dna", "protein", "both", "either"],
        default="protein",
        help="Choose to search genes versus DNA or protein database",
    )
    parser.add_argument(
        "--hit_type",
        choices=["all_hits", "top_hit"],
        default="top_hit",
        help="Transfer taxonomy of all hits or top hit per gene",
    )
    parser.add_argument(
        "--exclude_clades",
        type=str,
        nargs="+",
        help="List of clades to exclude (ex: s__1300164.4)",
    )
    parser.add_argument(
        "--bin_fract",
        type=float,
        default=0.7,
        help="Min fraction of genes in bin that agree with consensus taxonomy for bin annotation",
    )
    parser.add_argument(
        "--contig_fract",
        type=float,
        default=1.0,
        help="Min fraction of genes in that disagree with bin taxonomy for filtering",
    )
    parser.add_argument(
        "--allow_noclass",
        action="store_true",
        default=False,
        help="Allow a bin to be unclassfied and flag any classified contigs",
    )


def extract_homologs(tmp_dir):
    # create outdir
    seqs_dir = f"{tmp_dir}/markers"
    if not os.path.isdir(seqs_dir):
        os.makedirs(seqs_dir)
    # fetch best hits from hmmsearch
    gene_to_aln = utilities.fetch_hmm_best_hits(f"{tmp_dir}/phyeco.hmmsearch")
    # open output files
    outfiles = {}
    marker_ids = set([aln["qacc"] for aln in list(gene_to_aln.values())])
    for marker_id in marker_ids:
        outfiles[marker_id] = {}
        outfiles[marker_id]["ffn"] = open(f"{seqs_dir}/{marker_id}.ffn", "w")
        outfiles[marker_id]["faa"] = open(f"{seqs_dir}/{marker_id}.faa", "w")
    # write seqs
    for ext in ["ffn", "faa"]:
        in_path = f"{tmp_dir}/genes.{ext}"
        for id, seq in utilities.parse_fasta(in_path):
            if id in gene_to_aln:
                marker_id = gene_to_aln[id]["qacc"]
                seq = seq.rstrip("*")
                outfiles[marker_id][ext].write(">" + id + "\n" + seq + "\n")
    # close files
    for marker_id in outfiles:
        outfiles[marker_id]["ffn"].close()
        outfiles[marker_id]["faa"].close()


def align_homologs(db_dir, tmp_dir, seq_type, threads):
    aln_dir = f"{tmp_dir}/alns"
    if not os.path.exists(aln_dir):
        os.makedirs(aln_dir)
    seq_files = os.listdir(f"{tmp_dir}/markers")
    for file_index, seq_file in enumerate(seq_files):
        marker_id, ext = seq_file.split(".")
        if ext == "ffn" and seq_type == "protein":
            continue
        elif ext == "faa" and seq_type == "dna":
            continue
        if ext == "faa":
            program = "blastp" if ext == "faa" else "blastn"
            db_path = f"{db_dir}/phylo-markers/{program}/{marker_id}"
            query_path = f"{tmp_dir}/markers/{marker_id}.{ext}"
            out_path = f"{tmp_dir}/alns/{marker_id}.{ext}.m8"
            utilities.run_blastp(db_path, query_path, out_path, threads)
        else:
            program = "blastn"
            db_path = f"{db_dir}/phylo-markers/{program}/{marker_id}"
            query_path = f"{tmp_dir}/markers/{marker_id}.{ext}"
            out_path = f"{tmp_dir}/alns/{marker_id}.{ext}.m8"
            utilities.run_blastp(db_path, query_path, out_path, threads)


class Bin:
    def __init__(self):
        self.cons_taxon = None
        self.cons_fract = None
        self.cons_count = None
        self.genes = None

    def exclude_clades(self, clades):
        # loop over genes in bin
        for gene in list(self.genes.values()):
            # keep track of which annotations to exclude for each gene
            exclude_indexes = []
            # loop over annotations for each gene
            for index, annotation in enumerate(gene.annotations):
                is_match = any([c in annotation.taxon for c in clades])
                if is_match:
                    exclude_indexes.append(index)
            # delete annotations
            for index in exclude_indexes[::-1]:
                del gene.annotations[index]

    def only_keep_top_hits(self):
        for gene in list(self.genes.values()):
            if len(gene.annotations) > 1:
                max_score = max([a.score for a in gene.annotations])
                gene.annotations = [a for a in gene.annotations if a.score == max_score]

    def classify_taxonomy(self, allow_noclass=False, min_fraction=0.5):
        ranks = ["s", "g", "f", "o", "c", "p"]
        for rank_index, rank in enumerate(ranks):
            # get list of taxa across all genes
            # never count a taxon 2x per gene
            gene_taxa = []
            for gene in list(self.genes.values()):
                # do not count unclassified
                if not allow_noclass:
                    taxa = list(
                        set(
                            [
                                a.taxon[rank_index]
                                for a in gene.annotations
                                if a.taxon[rank_index]
                            ]
                        )
                    )
                # count unclassified
                else:
                    taxa = list(set([a.taxon[rank_index] for a in gene.annotations]))
                gene_taxa += taxa
            # skip rank where there are no annotations
            if len(gene_taxa) == 0:
                continue
            # get most common annotation
            cons_taxon, cons_count = Counter(gene_taxa).most_common()[0]
            # determine if enough genes have consensus annotation
            cons_fract = cons_count / float(len(self.genes))
            if cons_fract >= min_fraction:
                self.cons_taxon = cons_taxon
                self.cons_fract = cons_fract
                self.cons_count = cons_count
                self.rank = rank
                self.rank_index = rank_index
                return


class Marker:
    def __init__(self):
        self.id = None


class Gene:
    def __init__(self):
        self.id = None
        self.marker = None
        self.contig = None
        self.annotations = []


class Annotation:
    def __init__(self):
        self.taxon = [None] * 6
        self.score = None

    def add_taxon(self, genome_taxon, rank_index):
        for i in range(rank_index, 6):
            self.taxon[i] = genome_taxon.split(";")[i]


class Contig:
    def __init__(self):
        self.id = None
        self.genes = []
        self.genes_agree = 0
        self.genes_conflict = 0
        self.flagged = None

    def compare_taxonomy(self, bin):
        for gene in self.genes:
            agrees = []
            for annotation in gene.annotations:
                agrees.append(bin.cons_taxon == annotation.taxon[bin.rank_index])
            if any(agrees):
                self.genes_agree += 1
            else:
                self.genes_conflict += 1

    def flag(self, min_fract):
        if len(self.genes) == 0:
            self.flagged = None
        elif self.genes_conflict / float(len(self.genes)) >= min_fract:
            self.flagged = True
        else:
            self.flagged = False


def flag_contigs(db_dir, tmp_dir, args):
    # step 0. read in reference data files
    # cutoffs
    cutoffs = {}
    cutoffs_path = f"{db_dir}/phylo-markers/max_fscores.tsv"
    reader = csv.DictReader(open(cutoffs_path), delimiter="\t")
    for r in reader:
        key = (r["marker_id"], r["seq_type"], r["score_type"], r["taxlevel"])
        value = {"sensitive": r["cutoff_lower"], "strict": r["cutoff_upper"], "none": 0.0}
        cutoffs[key] = value
    # taxonomy
    taxonomy = {}
    taxonomy_path = f"{db_dir}/phylo-markers/genome_taxonomy.tsv"
    reader = csv.DictReader(open(taxonomy_path), delimiter="\t")
    for r in reader:
        taxonomy[r["genome_id"]] = r["taxonomy"]
    # clustered seqs
    clusters = {}
    for type in ["ffn", "faa"]:
        clusters[type] = {}
        for file in os.listdir(f"{db_dir}/phylo-markers/{type}"):
            if file.split(".")[-1] == "uc":
                with open(f"{db_dir}/phylo-markers/{type}/{file}") as f:
                    for l in f:
                        v = l.rstrip().split()
                        rep_id = v[-1]
                        seq_id = v[-2]
                        if v[0] == "S":
                            clusters[type][seq_id] = [seq_id]
                        elif v[0] == "H":
                            clusters[type][rep_id].append(seq_id)
    # step 1. determine if bin is archaea or bacteria; initialize domain-level markers
    # to do: normalize counts by site of marker gene sets
    marker_ids = set([])
    counts = {"bacteria": 0, "archaea": 0}
    for aln_file in os.listdir(f"{tmp_dir}/alns"):
        marker_id, seq_type, ext = aln_file.split(".")
        marker_ids.add(marker_id)
    for marker_id in marker_ids:
        if "B" in marker_id:
            counts["bacteria"] += 1
        elif "A" in marker_id:
            counts["archaea"] += 1
    domain = "bacteria" if counts["bacteria"] >= counts["archaea"] else "archaea"
    markers = {}
    for marker_id in marker_ids:
        if domain == "bacteria" and "B" in marker_id:
            markers[marker_id] = Marker()
            markers[marker_id].id = marker_id
            markers[marker_id].genes = []
        elif domain == "archaea" and "A" in marker_id:
            markers[marker_id] = Marker()
            markers[marker_id].id = marker_id
            markers[marker_id].genes = []
    # step 2. initialize marker genes found in bin
    bin = Bin()
    bin.genes = {}
    hmm_path = f"{tmp_dir}/phyeco.hmmsearch"
    for gene_id, aln in list(utilities.fetch_hmm_best_hits(hmm_path).items()):
        if aln["qacc"] not in markers:
            continue
        gene = Gene()
        gene.id = gene_id
        gene.contig = gene_id.rsplit("_", 1)[0]
        gene.marker = aln["qacc"]
        bin.genes[gene_id] = gene
        markers[aln["qacc"]].genes.append(gene)
    # annotate genes
    #    fetch all non-redundant taxonomic annotations for each gene
    seq_types = None
    if args["seq_type"] in ["both", "either"]:
        seq_types = ["ffn", "faa"]
    elif args["seq_type"] == "protein":
        seq_types = ["faa"]
    else:
        seq_types = ["ffn"]
    for seq_type in seq_types:
        for marker_id in markers:
            aln_path = tmp_dir + "/alns/" + marker_id + "." + seq_type + ".m8"
            for aln in utilities.parse_blast(aln_path):
                # fetch all unique taxonomies for target sequence
                # a sequence can have multiple taxonomies if it was clustered with another sequence
                genome_taxa = []
                for target_id in clusters[seq_type][aln["tname"]]:
                    genome_id = target_id.split("_")[0]
                    if genome_id not in taxonomy:
                        continue
                    elif taxonomy[genome_id] in genome_taxa:
                        continue
                    else:
                        genome_taxa.append(taxonomy[genome_id])
                # loop over ranks; stop when gene has been annotated
                for rank_index, rank in enumerate(["s", "g", "f", "o", "c", "p"]):
                    # decide to use ffn or faa at species level
                    if args["seq_type"] == "either":
                        if seq_type == "ffn" and rank != "s":
                            continue
                        elif seq_type == "faa" and rank == "s":
                            continue
                    # get minimum % identity cutoff for transfering taxonomy
                    #   if cutoff_type is None, indicates that no cutoff should be used
                    min_pid = cutoffs[marker_id, seq_type, "pid", rank][
                        args["cutoff_type"]
                    ]
                    if float(aln["pid"]) < float(min_pid):
                        continue
                    # add taxonomy
                    for genome_taxon in genome_taxa:
                        annotation = Annotation()
                        annotation.add_taxon(genome_taxon, rank_index)
                        annotation.score = float(aln["bitscore"])
                        bin.genes[aln["qname"]].annotations.append(annotation)
                    # stop when gene has been annotated at lowest rank
                    break
    # optionally remove annotations matching <exclude_clades>
    if args["exclude_clades"] is not None:
        bin.exclude_clades(args["exclude_clades"])
    # optionally take top hit only
    if args["hit_type"] == "top_hit":
        bin.only_keep_top_hits()
    # create None annotations for unannotated genes
    for gene in bin.genes.values():
        if len(gene.annotations) == 0:
            gene.annotations.append(Annotation())
    # classify bin
    bin.classify_taxonomy(args["allow_noclass"], args["bin_fract"])
    # flag contigs with discrepant taxonomy
    bin.contigs = {}
    for id, seq in utilities.parse_fasta(args["fna"]):
        bin.contigs[id] = Contig()
        bin.contigs[id].id = id
        bin.contigs[id].length = len(seq)
    for gene in bin.genes.values():
        bin.contigs[gene.contig].genes.append(gene)
    if bin.cons_taxon is not None:
        for contig in bin.contigs.values():
            contig.compare_taxonomy(bin)
            contig.flag(args["contig_fract"])
    # write results
    flagged_contigs = []
    for contig in bin.contigs.values():
        if contig.flagged:
            flagged_contigs.append(contig.id)
    return flagged_contigs


def main(args):
    utilities.add_tmp_dir(args)
    utilities.check_input(args)
    utilities.check_dependencies(["prodigal", "hmmsearch", "blastp", "blastn"])
    utilities.check_database(args)
    print("\u001b[1m" + "• Calling genes with Prodigal" + "\u001b[0m")
    utilities.run_prodigal(args["fna"], args["tmp_dir"])
    print(f"  all genes: {args['tmp_dir']}/genes.[ffn|faa]")
    print("\u001b[1m" + "\n• Identifying PhyEco phylogenetic marker genes with HMMER" + "\u001b[0m")
    utilities.run_hmmsearch(args["db"], args["tmp_dir"], args["tmp_dir"], args["threads"])
    extract_homologs(args["tmp_dir"])
    print(f"  hmm results: {args['tmp_dir']}/phyeco.hmmsearch")
    print(f"  marker genes: {args['tmp_dir']}/markers")
    print("\u001b[1m" + "\n• Performing pairwise BLAST alignment of marker genes against database" + "\u001b[0m")
    align_homologs(args["db"], args["tmp_dir"], args["seq_type"], args["threads"])
    print(f"  blast results: {args['tmp_dir']}/alns")
    print("\u001b[1m" + "\n• Finding taxonomic outliers" + "\u001b[0m")
    flagged = flag_contigs(args["db"], args["tmp_dir"], args)
    out = f"{args['tmp_dir']}/flagged_contigs"
    print(f"  {len(flagged)} flagged contigs: {out}")
    with open(out, "w") as f:
        for contig in flagged:
            f.write(contig + "\n")
