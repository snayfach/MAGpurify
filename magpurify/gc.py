#!/usr/bin/env python

import sys
import utility
import numpy as np
import argparse
import os

def fetch_args():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="MAGpurify: gc-content module: find contigs with outlier gc content"
	)
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('fna', type=str,
		help="""Path to input genome in FASTA format""")
	parser.add_argument('out', type=str,
		help="""Output directory to store results and intermediate files""")
	parser.add_argument('-t', dest='threads', type=int, default=1,
		help="""Number of CPUs to use (default=1)""")
	parser.add_argument('-d', dest='db', type=str,
		help="""Path to reference database
By default, the MAGPURIFY environmental variable is used""")
	parser.add_argument('--cutoff', type=float, default=15.75,
		help="""Cutoff (default=15.75)""")
	args = vars(parser.parse_args())
	return args

def add_defaults(args):
	args['cutoff'] = 15.75

def compute_gc(dna):
	counts = [dna.upper().count(base) for base in list('ACGT')]
	if sum(counts) > 0:
		return round(100*(counts[1]+counts[2])/float(sum(counts)),2)
	else:
		return 0.0

class Contig:
	def __init__(self):
		pass

def main():

	args = fetch_args()
	utility.add_tmp_dir(args)
	utility.check_input(args)
	utility.check_database(args)
	
	print "\n## Computing mean genome-wide GC content"
	contigs = {}
	for id, seq in utility.parse_fasta(args['fna']):
		contig = Contig()
		contig.id = id
		contig.seq = str(seq)
		contig.gc = compute_gc(seq)
		contigs[id] = contig
	mean = np.mean([c.gc for c in contigs.values()])
	std = np.std([c.gc for c in contigs.values()])
	
	print "\n## Computing per-contig deviation from mean"
	for contig in contigs.values():
		contig.values = {}
		contig.values['delta'] = abs(contig.gc - mean)
		contig.values['percent'] = 100 * abs(contig.gc - mean)/mean
		contig.values['z-score'] = abs(contig.gc - mean)/std
		
	print "\n## Identifying outlier contigs"
	flagged = []
	for contig in contigs.values():
		if contig.values['delta'] > args['cutoff']:
			flagged.append(contig.id)
	out = '%s/flagged_contigs' % args['tmp_dir']
	print ("   flagged contigs: %s" % out)
	with open(out, 'w') as f:
		for contig in flagged:
			f.write(contig+'\n')



