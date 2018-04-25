#!/usr/bin/env python

import sys
from utility import *
import numpy as np

def fetch_args():
	import argparse
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument('fna_path', type=str,
		help="""Path to input genome in FASTA format""")
	parser.add_argument('out_dir', type=str,
		help="""Output directory to store results and intermediate files""")
	parser.add_argument('-t', dest='threads', type=int, default=1,
		help="""Number of CPUs to use (default=1)""")
	parser.add_argument('-d', dest='db', type=str,
		help="""Path to reference database
By default, the IMAG_DB environmental variable is used""")
	parser.add_argument('--cutoff', type=float, default=15.75,
		help="""Cutoff (default=15.75)""")
	args = vars(parser.parse_args())
	if (not args['db']
			and 'IMAG_DB' in os.environ):
		args['db'] = os.environ['IMAG_DB']
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

def main(args):

	add_defaults(args)
	
	# store contigs
	contigs = {}
	for rec in parse_fasta(args['fna_path']):
		contig = Contig()
		contig.id = rec.id
		contig.seq = str(rec.seq)
		contig.gc = compute_gc(contig.seq)
		contigs[rec.id] = contig

	# compute mean and std
	mean = np.mean([c.gc for c in contigs.values()])
	std = np.std([c.gc for c in contigs.values()])
	
	# compute deviation
	for contig in contigs.values():
		contig.values = {}
		contig.values['delta'] = abs(contig.gc - mean)
		contig.values['percent'] = 100 * abs(contig.gc - mean)/mean
		contig.values['z-score'] = abs(contig.gc - mean)/std
		
	# flag contigs
	flagged = []
	for contig in contigs.values():
		if contig.values['delta'] > args['cutoff']:
			flagged.append(contig.id)

	return flagged

if __name__ == "__main__":
	
	args = fetch_args()

	check_input(args)
	check_output(args)
	check_dependencies([])
	check_database(args)

	main(args)

