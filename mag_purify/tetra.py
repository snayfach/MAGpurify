#!/usr/bin/env python

import sys, Bio.Seq
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from utility import *

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
By default, the IMAGEN_DB environmental variable is used""")
	parser.add_argument('--cutoff', type=float, default=0.06,
		help="""Cutoff (default=0.06)""")

	args = vars(parser.parse_args())
	if (not args['db']
			and 'IMAG_DB' in os.environ):
		args['db'] = os.environ['IMAG_DB']
	return args

def add_defaults(args):
	args['cutoff'] = 0.06

def init_kmers():
	tetra = {}
	for b1 in list('ACGT'):
		for b2 in list('ACGT'):
			for b3 in list('ACGT'):
				for b4 in list('ACGT'):
					kmer_fwd = ''.join([b1, b2, b3, b4])
					kmer_rev = str(Bio.Seq.Seq(kmer_fwd).reverse_complement())
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

def main(args):

	add_defaults(args)
	
	tmp_dir = '%s/tmp/tnf' % args['out_dir']
	if not os.path.exists(tmp_dir):
		os.makedirs(tmp_dir)

	kmer_counts = init_kmers()
	
	# init data
	contigs = {}
	for rec in parse_fasta(args['fna_path']):
		contig = Contig()
		contig.id = rec.id
		contig.seq = str(rec.seq)
		contig.kmers = kmer_counts.copy()
		contigs[rec.id] = contig

	# count kmers
	for contig in contigs.values():
		start, stop, step = 0, 4, 1
		while stop <= len(contig.seq):
			kmer_fwd = contig.seq[start:stop]
			kmer_rev = str(Bio.Seq.Seq(kmer_fwd).reverse_complement())
			if kmer_fwd in kmer_counts:
				contigs[rec.id].kmers[kmer_fwd] += 1
			elif kmer_rev in kmer_counts:
				contigs[rec.id].kmers[kmer_rev] += 1
			start += step
			stop += step

	# normalize
	for contig in contigs.values():
		total = float(sum(contig.kmers.values()))
		for kmer, count in contig.kmers.items():
			if total > 0:
				contig.kmers[kmer] = 100*count/total
			else:
				contig.kmers[kmer] = 0.00

	# pca
	df = pd.DataFrame(dict([(c.id, c.kmers) for c in contigs.values()]))
	pca = PCA(n_components=1)
	pca.fit(df)
	pc1 = pca.components_[0]

	# compute mean and std
	mean_pc = np.mean(pc1)
	std_pc = np.std(pc1)

	# store data
	for contig_id, contig_pc in zip(list(df.columns), pc1):
		contigs[contig_id].pc = contig_pc
		contigs[contig_id].values = {}
		contigs[contig_id].values['zscore'] = abs(contig_pc - mean_pc)/std_pc if std_pc > 0 else 0.0
		contigs[contig_id].values['delta'] = abs(contig_pc - mean_pc)
		contigs[contig_id].values['percent'] = 100*abs(contig_pc - mean_pc)/mean_pc

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
	check_dependencies(['blastn'])
	check_database(args)
	
	main(args)




