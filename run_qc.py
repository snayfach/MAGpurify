#!/usr/bin/env python

import sys
import argparse


def get_program():
    if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
        print('MAGpurify: Identify and remove incorrectly binned contigs from metagenome-assembled genomes')
        print('')
        print('Usage: run_qc.py <command> [options]')
        print('')
        print('Commands:')
        print('    phylo-markers: find taxonomic discordant contigs using db of phylogenetic marker genes')
        print('    clade-markers: find taxonomic discordant contigs using db of clade-specific marker genes')
        print('      conspecific: find contigs that fail to align to closely related genomes')
        print('       tetra-freq: find contigs with outlier tetranucleotide frequency')
        print('       gc-content: find contigs with outlier gc content')
        print('     known-contam: find contigs that match a database of known contaminants')
        print('        clean-bin: remove identified contigs from bin')
        print('')
        print('Note: use run_qc.py <command> -h to view usage for a specific command')
        quit()
    elif sys.argv[1] not in [
        'phylo-markers',
        'clade-markers',
        'conspecific',
        'tetra-freq',
        'gc-content',
        'known-contam',
        'clean-bin',
    ]:
        sys.exit("\nError: Unrecognized command: '%s'\n" % sys.argv[1])
        quit()
    else:
        return sys.argv[1]


if __name__ == "__main__":

    program = get_program()

    if program == 'phylo-markers':
        from magpurify import uscmg

        uscmg.main()

    if program == 'clade-markers':
        from magpurify import csmg

        csmg.main()

    if program == 'conspecific':
        from magpurify import conspecific

        conspecific.main()

    if program == 'tetra-freq':
        from magpurify import tetra

        tetra.main()

    if program == 'gc-content':
        from magpurify import gc

        gc.main()

    if program == 'known-contam':
        from magpurify import contam

        contam.main()

    if program == 'clean-bin':
        from magpurify import clean

        clean.main()

