#!/usr/bin/env python

import argparse
import os
import sys
from operator import itemgetter
from . import utility


def fetch_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=argparse.SUPPRESS,
        description="MAGpurify: clean module: remove flagged contigs from input genome",
    )
    parser.add_argument('program', help=argparse.SUPPRESS)
    parser.add_argument('fna', type=str, help="""Path to input genome in FASTA format""")
    parser.add_argument(
        'out',
        type=str,
        help="""Output directory to store results and intermediate files""",
    )
    args = vars(parser.parse_args())
    return args


def main():
    args = fetch_args()
    utility.check_input(args)
    print("\n## Reading genome bin")
    bin = {}
    for id, seq in utility.parse_fasta(args['fna']):
        bin[id] = seq
    bin_length = round(sum(len(_) for _ in bin.values()) / 1000, 2)
    print("   genome length: %s contigs, %s Kbp" % (len(bin), bin_length))

    print("\n## Reading flagged contigs")
    flagged_contigs = []
    programs = [
        'phylo-markers',
        'clade-markers',
        'conspecific',
        'tetra-freq',
        'gc-content',
        'coverage',
        'known-contam',
    ]
    for program in programs:
        path = '%s/%s/flagged_contigs' % (args['out'], program)
        if not os.path.exists(path):
            print("   %s: no output file found" % program)
        else:
            contigs = [_.rstrip() for _ in open(path)]
            bases = round(sum(len(bin[id]) for id in contigs) / 1000, 2)
            flagged_contigs += contigs
            print("   %s: %s contigs, %s Kbp" % (program, len(contigs), bases))
    flagged_contigs = list(set(flagged_contigs))
    flagged_length = round(sum(len(bin[id]) for id in flagged_contigs) / 1000, 2)
    print("\n## Removing flagged contigs")
    clean = bin.copy()
    for id in flagged_contigs:
        del clean[id]
    clean_length = round(sum(len(_) for _ in clean.values()) / 1000, 2)
    print("   removed: %s contigs, %s Kbp" % (len(flagged_contigs), flagged_length))
    print("   remains: %s contigs, %s Kbp" % (len(clean), clean_length))
    out = '%s/cleaned_bin.fna' % args['out']
    with open(out, 'w') as f:
        for id, seq in clean.items():
            f.write('>' + id + '\n' + seq + '\n')
    print("   cleaned bin: %s" % out)