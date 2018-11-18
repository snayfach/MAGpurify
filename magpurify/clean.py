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
		description="MAGpurify: clean module: remove flagged contigs from input genome"
	)
	parser.add_argument('program', help=argparse.SUPPRESS)
	parser.add_argument('fna', type=str,
		help="""Path to input genome in FASTA format""")
	parser.add_argument('out', type=str,
		help="""Output directory to store results and intermediate files""")
	args = vars(parser.parse_args())
	return args

def main():

	args = fetch_args()
	
	utility.check_input(args)
	
	print("\n## Reading genome bin")
	bin = {}
	for id, seq in utility.parse_fasta(args['fna']):
		bin[id] = seq
	bin_length = round(sum([len(_) for _ in bin.values()])/1000.0,2)
	print("   genome length: %s contigs, %s Kbp" % (len(bin), bin_length))

	print("\n## Reading flagged contigs")
	flagged_contigs = []
	programs = ['phylo-markers', 'clade-markers', 'conspecific', 'tetra-freq', 'gc-content', 'known-contam']
	for program in programs:
		path = '%s/%s/flagged_contigs' % (args['out'], program)
		if not os.path.exists(path):
			print("   %s: no output file found" % program)
		else:
			contigs = [_.rstrip() for _ in open(path)]
			bases = round(sum([len(bin[id]) for id in contigs])/1000.0,2)
			flagged_contigs += contigs
			print("   %s: %s contigs, %s Kbp" % (program, len(contigs), bases))
	flagged_contigs = list(set(flagged_contigs))
	flagged_length = round(sum([len(bin[id]) for id in flagged_contigs])/1000.0,2)

	print("\n## Removing flagged contigs")
	clean = bin.copy()
	for id in flagged_contigs:
		del clean[id]
	clean_length = round(sum([len(_) for _ in clean.values()])/1000.0,2)
	print("   removed: %s contigs, %s Kbp" % (len(flagged_contigs), flagged_length))
	print("   remains: %s contigs, %s Kbp" % (len(clean), clean_length))

	out = '%s/cleaned_bin.fna' % args['out']
	with open(out, 'w') as f:
		for id, seq in clean.items():
			f.write('>'+id+'\n'+seq+'\n')
	print("   cleaned bin: %s" % out)


#	check_input(args)
#	check_database(args)
#
#	contigs = {}
#	for rec in parse_fasta(args['fna_path']):
#		contig = Contig()
#		contig.id = rec.id
#		contig.seq = str(rec.seq)
#		contig.length = len(rec.seq)
#		contig.flags = []
#		contigs[rec.id] = contig
#
#	if args['modules'] and 'uscmg' in args['modules']:
#		from magpurify import uscmg
#		check_dependencies(['prodigal', 'blastp', 'hmmsearch'])
#		sys.stdout.write("discordant taxonomy based on universal-single-copy marker-genes\n")
#		flagged = magpurify.uscmg.main(args)
#		for id in flagged:
#			contigs[id].flags.append('uscmg')
#		sys.stdout.write("   %s contigs flagged\n" % len(flagged))
#
#	if args['modules'] and 'csmg' in args['modules']:
#		from magpurify import csmg
#		check_dependencies(['prodigal', 'lastal'])
#		sys.stdout.write("discordant taxonomy based on clade-specific marker-genes\n")
#		flagged = csmg.main(args)
#		for id in flagged:
#			contigs[id].flags.append('csmg')
#		sys.stdout.write("   %s contigs flagged\n" % len(flagged))
#
#	if args['modules'] and 'contdb' in args['modules']:
#		from magpurify import contam
#		check_dependencies(['blastn'])
#		sys.stdout.write("hits to database of contaminants\n")
#		flagged = contam.main(args)
#		for id in flagged:
#			contigs[id].flags.append('contdb')
#		sys.stdout.write("   %s contigs flagged\n" % len(flagged))
#
#	if args['modules'] and 'tnf' in args['modules']:
#		check_dependencies([])
#		sys.stdout.write("outlier tetranucleotide frequency\n")
#		from magpurify import tetra
#		flagged = tetra.main(args)
#		for id in flagged:
#			contigs[id].flags.append('tnf')
#		sys.stdout.write("   %s contigs flagged\n" % len(flagged))
#
#	if args['modules'] and 'gc' in args['modules']:
#		check_dependencies([])
#		sys.stdout.write("outlier gc content\n")
#		from magpurify import gc
#		flagged = gc.main(args)
#		for id in flagged:
#			contigs[id].flags.append('gcc')
#		sys.stdout.write("   %s contigs flagged\n" % len(flagged))
#
#	with open('%s/cleaned.fa' % args['out_dir'], 'w') as f:
#		for contig in contigs.values():
#			if len(contig.flags) > 0:
#				continue
#			f.write('>'+contig.id+'\n'+contig.seq+'\n')
#
#	with open('%s/summary.tsv' % args['out_dir'], 'w') as f:
#		fields = ['contig_id', 'length', 'flags']
#		f.write('\t'.join(fields)+'\n')
#		for contig in contigs.values():
#			values = [contig.id, str(contig.length), ','.join(contig.flags)]
#			f.write('\t'.join(values)+'\n')
