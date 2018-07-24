#!/usr/bin/env python

import os
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

#	parser.add_argument('query')
#	parser.add_argument('outdir')
	parser.add_argument('--mash_dist', type=float, default=0.05,
		help="Mash distance to reference genomes (default=0.05)")
	parser.add_argument('--max_genomes', type=int, default=1000,
		help="Max number of genomes to use (default=1000)")
	parser.add_argument('--min_genomes', type=int, default=1,
		help="Min number of genomes to use (default=1)")
	parser.add_argument('--contig_aln', type=float, default=0.50,
		help="Minimum fraction of contig aligned to reference (default=0.50)")
	parser.add_argument('--contig_pid', type=float, default=95.0,
		help="Minimum percent identity of contig aligned to reference (default=95.0)")
	parser.add_argument('--hit_rate', type=float, default=0.00,
		help="Hit rate for flagging contigs (default=0.00)")
	parser.add_argument('--exclude', default='',
		help="Comma-separated list of references to exclude")
	args = vars(parser.parse_args())
	args['exclude'] = args['exclude'].split(',')
	return args

def run_mash(db_dir, fna_path, tmp_dir, threads=1):
	db_path = '%s/sketch.msh' % db_dir
	out_path = '%s/mash.dist' % tmp_dir
	command = "mash dist -p %s -d 0.25 %s %s > %s" % (threads, fna_path, db_path, out_path)
	out, err = run_process(command)

	with open(tmp_dir+'/id_map.tsv', 'w') as f:
		for id, rec in enumerate(parse_mash(out_path)):
			f.write(str(id)+'\t'+rec['target']+'\n')

def find_conspecific(tmp_dir, max_dist, exclude, max_targets=10):

	targets = []
	for rec in parse_mash('%s/mash.dist' % tmp_dir):
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

	from operator import itemgetter
	targets = sorted(targets, key=itemgetter(1))

	return targets[0:max_targets]

def blastn(query, target, outdir, id):
	out_path = outdir+'/'+id+'.m8'
	if not os.path.exists(out_path):
		cmd = "blastn -outfmt '6 std qlen slen' "
		cmd += "-max_target_seqs 1 -max_hsps 1 "
		cmd += "-query %s -subject %s " % (query, target)
		out, err = run_process(cmd)
		open(out_path, 'w').write(out)
	else:
		out = open(out_path).read()
	return out

def id_blast_hits(blast_out, min_aln, min_pid):
	blast_hits = set([])
	for rec in BlastIO.parse_out(blast_out):
		if rec['qcov'] < min_aln:
			continue
		elif rec['pid'] < min_pid:
			continue
		else:
			blast_hits.add(rec['qname'])
	return blast_hits

def find_contig_targets():
	
	target_to_id = {}
	for line in open(args['tmp_dir']+'/id_map.tsv'):
		id, target = line.split()
		target_to_id[target] = id
	
	contigs = dict([(rec.id,{'hits':0,'len':len(rec.seq)}) for rec in parse_fasta(args['fna_path'])])

	for genome_path, mash_dist in genomes:
	
		blast_out = blastn(args['fna_path'], genome_path, args['tmp_dir'], target_to_id[target])
		
		hits = id_blast_hits(blast_out, args['contig_aln'], args['contig_pid'])
		print blast_out
		quit()
		
		for contig in hits:
			contigs[contig]['hits'] += 1
	
	for id in contigs:
		hit_rate = contigs[id]['hits']/float(len(genomes))
		contigs[id]['hit_rate'] = hit_rate

	return contigs

def flag_contigs():
	flagged = []
	for id in contigs:
		if contigs[id]['hit_rate'] > args['hit_rate']:
			continue
		else:
			flagged.append(id)
	return flagged

if __name__ == "__main__":

	args = fetch_args()
	
	tmp_dir = '%s/tmp/conspecific' % args['out_dir']
	if not os.path.exists(tmp_dir):
		os.makedirs(tmp_dir)

	check_input(args)
	check_output(args)
	check_dependencies(['mash'])
	check_database(args)

	run_mash(args['db'], args['fna_path'], args['tmp_dir'], args['threads'])
	
	genomes = find_conspecific(args['tmp_dir'], args['mash_dist'], args['exclude'], args['max_genomes'])

	if len(genomes) < args['min_genomes']:
		sys.stderr.write("Genome targets: %s (FAIL)\n" % len(genomes))
		sys.stderr.write("Flagged contigs: 0\n")
		sys.stderr.write("Flagged length: 0\n")
	else:
		contigs = find_contig_targets()
		total_seqs = len(contigs)
		total_len = sum([c['len'] for c in contigs.values()])
		flagged = flag_contigs()
		flagged_seqs = 0
		flagged_len = 0
		for id in flagged:
			sys.stdout.write(id+"\n")
			flagged_seqs += 1
			flagged_len += contigs[id]['len']
		sys.stderr.write("Genome targets: %s (PASS)\n" % len(genome_targets))
		sys.stderr.write("Flagged contigs: %s/%s\n" % (flagged_seqs, total_seqs))
		sys.stderr.write("Flagged length: %s/%s\n" % (flagged_len, total_len))





