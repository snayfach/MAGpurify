#!/usr/bin/env python

import os
import sys
import utility
import argparse

def fetch_args():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="MAGpurify: known-contam module: find contigs that match a database of known contaminants"
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
By default, the IMAGEN_DB environmental variable is used""")
	parser.add_argument('--pid', type=float, default=98,
		help="""Minimum % identity to reference (default=98)""")
	parser.add_argument('--evalue', type=float, default=1e-5,
		help="""Maximum evalue (default=1e-5)""")
	parser.add_argument('--qcov', type=float, default=25,
		help="""Minimum percent query coverage (default=25)""")
	args = vars(parser.parse_args())
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
	utility.run_process(cmd)

def main():
	
	args = fetch_args()
	
	utility.add_tmp_dir(args)
	utility.check_input(args)
	utility.check_dependencies(['blastn'])
	utility.check_database(args)
	
	tmp_dir = '%s/%s' % (args['out'], args['program'])
	if not os.path.exists(args['tmp_dir']):
		os.makedirs(args['tmp_dir'])

	print ("\n## Searching database with BLASTN")
	for target in ['hg38', 'phix']:
		db = '%s/known-contam/%s/%s' % (args['db'], target, target)
		out = '%s/%s.m8' % (args['tmp_dir'], target)
		run_blastn(args['fna'], db, out, args['threads'], args['qcov'], args['pid'], args['evalue'])

	print ("\n## Identifying contigs with hits to db")
	flagged = set([])
	for target in ['hg38', 'phix']:
		out = '%s/%s.m8' % (args['tmp_dir'], target)
		for r in utility.parse_blast(out):
			flagged.add(r['qname'])
	flagged = list(flagged)
	out = '%s/flagged_contigs' % args['tmp_dir']
	print ("   flagged contigs: %s" % out)
	with open(out, 'w') as f:
		for contig in flagged:
			f.write(contig+'\n')




