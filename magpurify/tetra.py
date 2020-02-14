#!/usr/bin/env python

import argparse
import itertools
import os
import sys
import Bio.Seq
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from . import utility


def fetch_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        description="MAGpurify: tetra-freq module: find contigs with outlier tetranucleotide frequency",
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
        '--cutoff', type=float, default=0.06, help="""Cutoff (default=0.06)"""
    )
    args = vars(parser.parse_args())
    return args


def add_defaults(args):
    args['cutoff'] = 0.06


def init_kmers():
    tetra = {}
    for i in itertools.product("ACGT", repeat=4):
        kmer_fwd = ''.join(i)
        kmer_rev = utility.reverse_complement(kmer_fwd)
        if kmer_fwd in tetra:
            continue
        elif kmer_rev in tetra:
            continue
        else:
            tetra[kmer_fwd] = 0
    return tetra


class Contig:
    def __init__(self):
        pass


def main():

    args = fetch_args()
    utility.add_tmp_dir(args)
    utility.check_input(args)
    utility.check_dependencies(['blastn'])

    print("\n## Counting tetranucleotides")
    # init data
    kmer_counts = init_kmers()
    contigs = {}
    for id, seq in utility.parse_fasta(args['fna']):
        contig = Contig()
        contig.id = id
        contig.seq = str(seq)
        contig.kmers = kmer_counts.copy()
        contigs[id] = contig

    # count kmers
    for contig in contigs.values():
        start, stop, step = 0, 4, 1
        while stop <= len(contig.seq):
            kmer_fwd = contig.seq[start:stop]
            if kmer_fwd in kmer_counts:
                contig.kmers[kmer_fwd] += 1
            else:
                kmer_rev = utility.reverse_complement(kmer_fwd)
                contig.kmers[kmer_rev] += 1
            start += step
            stop += step

    print("\n## Normalizing counts")
    for contig in contigs.values():
        total = float(sum(contig.kmers.values()))
        for kmer, count in contig.kmers.items():
            if total > 0:
                contig.kmers[kmer] = 100 * count / total
            else:
                contig.kmers[kmer] = 0.00

    print("\n## Performing PCA")
    df = pd.DataFrame(dict([(c.id, c.kmers) for c in contigs.values()]))
    pca = PCA(n_components=1)
    pca.fit(df)
    pc1 = pca.components_[0]

    print(
        "\n## Computing per-contig deviation from the mean along the first principal component"
    )
    mean_pc = np.mean(pc1)
    for contig_id, contig_pc in zip(list(df.columns), pc1):
        contigs[contig_id].pc = contig_pc
        contigs[contig_id].values = {}
        contigs[contig_id].values['delta'] = abs(contig_pc - mean_pc)

    print("\n## Identifying outlier contigs")
    flagged = []
    for contig in contigs.values():
        if contig.values['delta'] > args['cutoff']:
            flagged.append(contig.id)
    out = '%s/flagged_contigs' % args['tmp_dir']
    print("   flagged contigs: %s" % out)
    with open(out, 'w') as f:
        for contig in flagged:
            f.write(contig + '\n')
