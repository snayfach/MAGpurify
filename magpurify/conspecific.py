#!/usr/bin/env python

import os
import utility
import sys
import argparse
from operator import itemgetter

def fetch_args():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="MAGpurify: conspecific module: identify contigs that fail to align to closely related genomes"
	)
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('fna', type=str,
		help="""Path to input genome in FASTA format""")
	parser.add_argument('out', type=str,
		help="""Output directory to store results and intermediate files""")
	parser.add_argument('--threads', type=int, default=1, metavar='INT',
		help="""Number of CPUs to use (default=1)""")
	parser.add_argument('--mash-sketch', type=str, metavar='PATH', required=True,
		help="""Path to Mash sketch of reference genomes""")
	parser.add_argument('--mash-dist', type=float, default=0.05, metavar='FLOAT',
		help="Mash distance to reference genomes (default=0.05)")
	parser.add_argument('--max-genomes', type=int, default=25, metavar='INT',
		help="Max number of genomes to use (default=25)")
	parser.add_argument('--min-genomes', type=int, default=1, metavar='INT',
		help="Min number of genomes to use (default=1)")
	parser.add_argument('--contig-aln', type=float, default=0.50, metavar='FLOAT',
		help="Minimum fraction of contig aligned to reference (default=0.50)")
	parser.add_argument('--contig-pid', type=float, default=95.0, metavar='FLOAT',
		help="Minimum percent identity of contig aligned to reference (default=95.0)")
	parser.add_argument('--hit-rate', type=float, default=0.00, metavar='FLOAT',
		help="Hit rate for flagging contigs (default=0.00)")
	parser.add_argument('--exclude', default='', metavar='STR',
		help="Comma-separated list of references to exclude")
	args = vars(parser.parse_args())
	args['exclude'] = args['exclude'].split(',')
	return args

def run_mash(mash_sketch, fna_path, tmp_dir, threads=1):
	out_path = '%s/mash.dist' % tmp_dir
	command = "mash dist -p %s -d 0.25 %s %s > %s" % (threads, fna_path, mash_sketch, out_path)
	out, err = utility.run_process(command)
	with open(tmp_dir+'/id_map.tsv', 'w') as f:
		for id, rec in enumerate(utility.parse_mash(out_path)):
			f.write(str(id)+'\t'+rec['target']+'\n')

def find_conspecific(tmp_dir, max_dist, exclude, max_targets=10):

	targets = []
	for rec in utility.parse_mash('%s/mash.dist' % tmp_dir):
		if rec['query'] == rec['target']:
			continue
		elif rec['dist'] > max_dist:
			continue
		elif rec['pvalue'] > 1e-3:
			continue
		elif rec['target'] in exclude:
			continue
		else:
			targets.append([rec['target'], rec['dist']])

	targets = sorted(targets, key=itemgetter(1))

	return targets[0:max_targets]

def blastn(query, target, outdir, id):
	out_path = outdir+'/'+id+'.m8'
	if not os.path.exists(out_path):
		cmd = "blastn -outfmt '6 std qlen slen' "
		cmd += "-max_target_seqs 1 -max_hsps 1 "
		cmd += "-query %s -subject %s " % (query, target)
		out, err = utility.run_process(cmd)
		open(out_path, 'w').write(out)
	else:
		out = open(out_path).read()
	return out

def id_blast_hits(blast_out, min_aln, min_pid):
	blast_hits = set([])
	for rec in utility.parse_blast(blast_out, type='string'):
		if rec['qcov'] < min_aln:
			continue
		elif rec['pid'] < min_pid:
			continue
		else:
			blast_hits.add(rec['qname'])
	return blast_hits

def align_contigs(args, genomes):
	target_to_id = {}
	for line in open(args['tmp_dir']+'/id_map.tsv'):
		id, target = line.split()
		target_to_id[target] = id
	alignments = []
	for genome_path, mash_dist in genomes:
		blast_out = blastn(args['fna'], genome_path, args['tmp_dir'], target_to_id[target])
		alignments.append(blast_out)
	return alignments
		

def find_contig_targets(args, alignments):

	contigs = dict([(id,{'hits':0,'len':len(seq)}) for id, seq in utility.parse_fasta(args['fna'])])

	for alns in alignments:
		hits = id_blast_hits(alns, args['contig_aln'], args['contig_pid'])
		for contig in hits:
			contigs[contig]['hits'] += 1
	
	for id in contigs:
		hit_rate = contigs[id]['hits']/float(len(alignments))
		contigs[id]['hit_rate'] = hit_rate

	return contigs

def flag_contigs(args, contigs):
	flagged = []
	for id in contigs:
		if contigs[id]['hit_rate'] > args['hit_rate']:
			continue
		else:
			flagged.append(id)
	return flagged

def main():

	args = fetch_args()
	
	utility.add_tmp_dir(args)
	utility.check_input(args)
	utility.check_dependencies(['mash'])
	if not os.path.exists(args['mash_sketch']):
		sys.exit("\nError: mash sketch '%s' not found\n" % args['mash_sketch'])

	print ("\n## Finding conspecific genomes in database")
	run_mash(args['mash_sketch'], args['fna'], args['tmp_dir'], args['threads'])
	genomes = find_conspecific(args['tmp_dir'], args['mash_dist'], args['exclude'], args['max_genomes'])

	if len(genomes) < args['min_genomes']:
		sys.stderr.write("Genome targets: %s (FAIL)\n" % len(genomes))
		sys.stderr.write("Flagged contigs: 0\n")
		sys.stderr.write("Flagged length: 0\n")
		return

	print ("\n## Performing pairwise alignment of contigs")
	alignments = align_contigs(args, genomes)
	
	print ("\n## Summarizing pairwise alignments")
	contigs = find_contig_targets(args, alignments)

	print ("\n## Identifying contigs with no conspecific alignments")
	flagged = flag_contigs(args, contigs)
	out = '%s/flagged_contigs' % args['tmp_dir']
	print ("   flagged contigs: %s" % out)
	with open(out, 'w') as f:
		for contig in flagged:
			f.write(contig+'\n')

	print ("\n## Reporting results")
	total_seqs = len(contigs)
	total_len = sum([c['len'] for c in contigs.values()])
	flagged_seqs = 0
	flagged_len = 0
	for id in flagged:
		#sys.stdout.write(id+"\n")
		flagged_seqs += 1
		flagged_len += contigs[id]['len']
	sys.stderr.write("Genome targets: %s (PASS)\n" % len(genomes))
	sys.stderr.write("Flagged contigs: %s/%s\n" % (flagged_seqs, total_seqs))
	sys.stderr.write("Flagged length: %s/%s\n" % (flagged_len, total_len))




