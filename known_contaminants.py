#!/usr/bin/env python

import sys
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
	parser.add_argument('--pid', type=float, default=98,
		help="""Minimum % identity to reference (default=98)""")
	parser.add_argument('--evalue', type=float, default=1e-5,
		help="""Maximum evalue (default=1e-5)""")
	parser.add_argument('--qcov', type=float, default=25,
		help="""Minimum percent query coverage (default=25)""")
	args = vars(parser.parse_args())
	if (not args['db']
			and 'IMAGEN_DB' in os.environ):
		args['db'] = os.environ['IMAGEN_DB']
	return args

def run_blastn(query, db, out, threads, qcov=25, pid=98, evalue=1e-5):
	cmd = "blastn "
	cmd += "-query %s " % query
	cmd += "-db %s " % db
	cmd += "-out %s " % out
	cmd += "-outfmt '6 std qlen slen' "
	cmd += "-max_target_seqs 1 "
	cmd += "-max_hsps 1 "
	cmd += "-qcov_hsp_perc %s " % qcov
	cmd += "-perc_identity %s " % pid
	cmd += "-evalue %s " % evalue
	cmd += "-num_threads %s " % threads
	print cmd
	run_process(cmd)

def flag_contigs(args):
	flagged = set([])
	for target in ['hg38', 'phix']:
		db = '%s/contaminants/%s/%s' % (args['db'], target, target)
		out = '%s/%s.m8' % (args['out_dir'], target)
		run_blastn(args['fna_path'], db, out, args['threads'], args['qcov'], args['pid'], args['evalue'])
		for r in parse_blast(out):
			flagged.add(r['qname'])
	return flagged


if __name__ == "__main__":
	
	args = fetch_args()

	check_input(args)
	check_output(args)
	check_dependencies(['blastn'])
	check_database(args)

	flagged = flag_contigs(args)

	for contig_id in flagged:
		sys.stdout.write(contig_id+'\n')


