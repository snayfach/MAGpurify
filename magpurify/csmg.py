#!/usr/bin/env python

import sys, os, copy, collections, operator
from . import utility
import argparse

ranks = ['k', 'p', 'c', 'o', 'f', 'g', 's']
rank_names = {
    'k': 'kingdom',
    'p': 'phylum',
    'c': 'class',
    'o': 'order',
    'f': 'family',
    'g': 'genus',
    's': 'species',
}


def fetch_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        description="MAGpurify: clade-markers module: find taxonomic discordant contigs using db of clade-specific marker genes",
    )
    parser.add_argument('program', help=argparse.SUPPRESS)
    parser.add_argument('fna', type=str, help="""Path to input genome in FASTA format""")
    parser.add_argument(
        'out',
        type=str,
        help="""Output directory to store results and intermediate files""",
    )
    parser.add_argument(
        '-t',
        dest='threads',
        type=int,
        default=1,
        help="""Number of CPUs to use (default=1)""",
    )
    parser.add_argument(
        '-d',
        dest='db',
        type=str,
        help="""Path to reference database
By default, the MAGPURIFY environmental variable is used""",
    )
    parser.add_argument(
        '-e',
        '--exclude_clades',
        type=str,
        help="""Comma separated list of clades to exclude (ex: s__Variovorax_sp_CF313)""",
    )
    parser.add_argument(
        '-b',
        '--min_bin_fract',
        type=float,
        default=0.6,
        help="""Min fraction of bin length supported by contigs that agree with consensus taxonomy (default=0.6)""",
    )
    parser.add_argument(
        '-c',
        '--min_contig_fract',
        type=float,
        default=0.75,
        help="""Min fraction of classified contig length that agree with consensus taxonomy (default=0.75)""",
    )
    parser.add_argument(
        '-g',
        '--min_gene_fract',
        type=float,
        default=0.0,
        help="""Min fraction of classified genes that agree with consensus taxonomy (default=0.0)""",
    )
    parser.add_argument(
        '-m',
        '--min_genes',
        type=float,
        default=None,
        help="""Min number of genes that agree with consensus taxonomy (default=rank-specific-cutoffs)""",
    )
    parser.add_argument(
        '-l',
        '--lowest_rank',
        choices=['s', 'g', 'f', 'o', 'c', 'p', 'k'],
        help="""Lowest rank for bin classification""",
    )

    args = vars(parser.parse_args())
    return args


def add_defaults(args):
    args['lowest_rank'] = None
    args['min_genes'] = None
    args['min_gene_fract'] = 0.0
    args['min_contig_fract'] = 0.75
    args['min_bin_fract'] = 0.6
    args['exclude_clades'] = None


def read_ref_taxonomy(db_dir):
    ref_taxonomy = {}
    inpath = '%s/clade-markers/taxonomy.tsv' % db_dir
    for line in open(inpath):
        ref_id, taxonomy = line.rstrip().split()
        ref_taxonomy[ref_id] = taxonomy
    return ref_taxonomy


def flag_contigs(contigs, bin):
    for bin_taxon in bin.taxonomy:
        bin_rank = bin_taxon.split('__')[0]
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
                cons_taxon = 'NA'
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


def main():

    args = fetch_args()
    utility.add_tmp_dir(args)
    utility.check_input(args)
    utility.check_database(args)

    print("\n## Reading database info")
    ref_taxonomy = read_ref_taxonomy(args['db'])
    taxon_to_taxonomy = {}
    for taxonomy in set(ref_taxonomy.values()):
        for taxon in taxonomy.split('|'):
            taxon_to_taxonomy[taxon] = taxonomy
    min_pid = {'k': 57, 'p': 77, 'c': 82, 'o': 86, 'f': 87, 'g': 91, 's': 96}
    if args['min_genes'] is not None:
        args['min_genes'] = dict([(r, args['min_genes']) for r in ranks])
    else:
        args['min_genes'] = {
            'k': 237,
            'p': 44,
            'c': 30,
            'o': 24,
            'f': 22,
            'g': 20,
            's': 19,
        }

    print("\n## Calling genes with Prodigal")
    utility.run_prodigal(args['fna'], args['tmp_dir'])
    print("   all genes: %s/genes.[ffn|faa]" % args['tmp_dir'])

    print(
        "\n## Performing pairwise alignment of genes against MetaPhlan2 db of clade-specific genes"
    )
    utility.run_lastal(args['db'], args['tmp_dir'], args['threads'])
    print("   alignments: %s/genes.m8" % args['tmp_dir'])

    print("\n## Finding top hits to db")
    genes = {}
    for aln in utility.parse_last(args['tmp_dir'] + '/genes.m8'):

        # clade exclusion
        ref_taxa = ref_taxonomy[aln['tid']].split('|')
        if args['exclude_clades'] and any(
            [taxon in ref_taxa for taxon in args['exclude_clades'].split(',')]
        ):
            continue

        # initialize gene
        if aln['qid'] not in genes:
            genes[aln['qid']] = Gene()
            genes[aln['qid']].id = aln['qid']
            genes[aln['qid']].contig_id = aln['qid'].rsplit('_', 1)[0]

        # get top alignments
        if genes[aln['qid']].aln is None:
            genes[aln['qid']].aln = aln
            genes[aln['qid']].ref_taxa = ref_taxa
        elif float(aln['score']) > float(genes[aln['qid']].aln['score']):
            genes[aln['qid']].ref_taxa = ref_taxa
    print("   %s genes with a database hit" % len(genes))

    print("\n## Classifying genes at each taxonomic rank")
    counts = {}
    for gene in genes.values():
        for ref_taxon in gene.ref_taxa:
            rank = ref_taxon.split('__')[0]
            if rank not in counts:
                counts[rank] = 0
            if rank == 't':
                continue
            elif float(gene.aln['pid']) < min_pid[rank]:
                continue
            elif gene.aln['qcov'] < 0.4:
                continue
            elif gene.aln['tcov'] < 0.4:
                continue
            gene.taxa[rank] = ref_taxon
            counts[rank] += 1
    for rank in ranks:
        print("   %s: %s classified genes" % (rank_names[rank], counts[rank]))

    print("\n## Taxonomically classifying contigs")
    contigs = {}
    for id, seq in utility.parse_fasta(args['fna']):
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
    print("   total contigs: %s" % len(contigs))
    for rank in ranks:
        print("   %s: %s classified contigs" % (rank_names[rank], counts[rank]))

    print("\n## Taxonomically classifying genome")
    bin = Bin()
    bin.classify(
        contigs,
        args['min_bin_fract'],
        args['min_contig_fract'],
        args['min_gene_fract'],
        args['min_genes'],
        args['lowest_rank'],
    )
    print("   consensus taxon: %s" % bin.cons_taxon)

    print("\n## Identifying taxonomically discordant contigs")
    if bin.cons_taxon is not None:
        bin.rank_index = (
            taxon_to_taxonomy[bin.cons_taxon].split('|').index(bin.cons_taxon)
        )
        bin.taxonomy = taxon_to_taxonomy[bin.cons_taxon].split('|')[
            0 : bin.rank_index + 1
        ]
        flag_contigs(contigs, bin)
    flagged = []
    for contig in contigs.values():
        if contig.flagged:
            flagged.append(contig.id)
    out = '%s/flagged_contigs' % args['tmp_dir']
    print("   flagged contigs: %s" % out)
    with open(out, 'w') as f:
        for contig in flagged:
            f.write(contig + '\n')

