#!/usr/bin/env python


import os, sys, csv
import utility
from collections import Counter
import argparse

def parse_args():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS,
		description="MAGpurify: phylo-markers module: find taxonomic discordant contigs using db of phylogenetic marker genes"
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
By default, the MAGPURIFYDB environmental variable is used""")
	parser.add_argument('-c', dest='continue', action='store_true', default=False,
		help="""Go straight to quality estimation and skip all previous steps""")
	parser.add_argument('--max_target_seqs', type=int, default=1,
		help="""Maximum number of targets reported in BLAST table (default=1)""")
	parser.add_argument('--cutoff_type', choices=['strict', 'sensitive', 'none'], default='strict',
		help="""Use strict or sensitive %%ID cutoff for taxonomically annotating genes (default=strict)""")
	parser.add_argument('--seq_type', choices=['dna', 'protein', 'both', 'either'], default='protein',
		help="""Choose to search genes versus DNA or protein database (default=protein)""")
	parser.add_argument('--hit_type', choices=['all_hits', 'top_hit'], default='top_hit',
		help="""Transfer taxonomy of all hits or top hit per gene (default=top_hit)""")
	parser.add_argument('--exclude_clades', type=str,
		help="""Comma separated list of clades to exclude (ex: s__1300164.4)""")
	parser.add_argument('--bin_fract', type=float, default=0.7,
		help="""Min fraction of genes in bin that agree with consensus taxonomy for bin annotation (default=0.7)""")
	parser.add_argument('--contig_fract', type=float, default=1.0,
		help="""Min fraction of genes in that disagree with bin taxonomy for filtering (default=1.0)""")
	parser.add_argument('--allow_noclass', action='store_true', default=False,
		help="""Allow a bin to be unclassfied and flag any classified contigs (default=False)""")
	args = vars(parser.parse_args())
	return args

def add_defaults(args):
	args['continue'] = False
	args['max_target_seqs'] = 1
	args['cutoff_type'] = 'strict'
	args['seq_type'] = 'protein'
	args['hit_type'] = 'top_hit'
	args['exclude_clades'] = None
	args['bin_fract'] = 0.7
	args['contig_fract'] = 1.0
	args['allow_noclass'] = False

def extract_homologs(tmp_dir):
	# create outdir
	seqs_dir = "%s/seqs" % tmp_dir
	if not os.path.isdir(seqs_dir):
		os.makedirs(seqs_dir)
	# fetch best hits from hmmsearch
	gene_to_aln = utility.fetch_hmm_best_hits("%s/phyeco.hmmsearch" % tmp_dir)
	# open output files
	outfiles = {}
	marker_ids = set([aln['qacc'] for aln in gene_to_aln.values()])
	for marker_id in marker_ids:
		outfiles[marker_id] = {}
		outfiles[marker_id]['ffn'] = open("%s/%s.ffn" % (seqs_dir, marker_id), "w")
		outfiles[marker_id]['faa'] = open("%s/%s.faa" % (seqs_dir, marker_id), "w")
	# write seqs
	for ext in ['ffn', 'faa']:
		in_path = "%s/genes.%s" % (tmp_dir, ext)
		for id, seq in utility.parse_fasta(in_path):
			if id in gene_to_aln:
				marker_id = gene_to_aln[id]['qacc']
				seq = seq.upper().rstrip('*')
				outfiles[marker_id][ext].write('>'+id+'\n'+seq+'\n')
	# close files
	for marker_id in outfiles:
		outfiles[marker_id]['ffn'].close()
		outfiles[marker_id]['faa'].close()

def align_homologs(db_dir, tmp_dir, seq_type, threads):
	aln_dir = "%s/alns" % tmp_dir
	if not os.path.exists(aln_dir):
		os.makedirs(aln_dir)
	seq_files = os.listdir("%s/seqs" % tmp_dir)
	for file_index, seq_file in enumerate(seq_files):
		marker_id, ext = seq_file.split('.')
		if ext == 'ffn' and seq_type == 'protein':
			continue
		elif ext == 'faa' and seq_type == 'dna':
			continue
		if ext == 'faa':
			program = 'blastp' if ext == 'faa' else 'blastn'
			db_path = '%s/uscmg/%s/%s' % (db_dir, program, marker_id)
			query_path = '%s/seqs/%s.%s' % (tmp_dir, marker_id, ext)
			out_path = '%s/alns/%s.%s.m8' % (tmp_dir, marker_id, ext)
			utility.run_blastp(db_path, query_path, out_path, threads)
		else:
			program = 'blastn'
			db_path = '%s/uscmg/%s/%s' % (db_dir, program, marker_id)
			query_path = '%s/seqs/%s.%s' % (tmp_dir, marker_id, ext)
			out_path = '%s/alns/%s.%s.m8' % (tmp_dir, marker_id, ext)
			utility.run_blastp(db_path, query_path, out_path, threads)

class Bin:
	def __init__(self):
		self.cons_taxon = None
		self.cons_fract = None
		self.cons_count = None
		self.genes = None
	
	def exclude_clades(self, clades):
	
		# loop over genes in bin
		for gene in self.genes.values():

			# keep track of which annotations to exclude for each gene
			exclude_indexes = []
			
			# loop over annotations for each gene
			for index, annotation in enumerate(gene.annotations):
				is_match = any([c in annotation.taxon for c in clades.split(',')])
				if is_match:
					exclude_indexes.append(index)
			
			# delete annotations
			for index in exclude_indexes[::-1]:
				del gene.annotations[index]
	
	def only_keep_top_hits(self):
		for gene in self.genes.values():
			if len(gene.annotations) > 1:
				max_score = max([a.score for a in gene.annotations])
				gene.annotations = [a for a in gene.annotations if a.score == max_score]
	
	def classify_taxonomy(self, allow_noclass=False, min_fraction=0.5):
		"""
		determine taxon and rank of bin based on lowest rank
		where at least <min_fraction> of genes match the consensus annotation
		"""
		ranks = ['s','g','f','o','c','p']
		for rank_index, rank in enumerate(ranks):
			
			# get list of taxa across all genes
			# never count a taxon 2x per gene
			gene_taxa = []
			for gene in self.genes.values():
				# do not count unclassified
				if not allow_noclass:
					taxa = list(set([a.taxon[rank_index] for a in gene.annotations if a.taxon[rank_index]]))
				# count unclassified
				else:
					taxa = list(set([a.taxon[rank_index] for a in gene.annotations]))
				gene_taxa += taxa

			# skip rank where there are no annotations
			if len(gene_taxa) == 0:
				continue

			# get most common annotation
			cons_taxon, cons_count = Counter(gene_taxa).most_common()[0]

			# determine if enough genes have consensus annotation
			cons_fract = cons_count/float(len(self.genes))
			if cons_fract >= min_fraction:
				self.cons_taxon = cons_taxon
				self.cons_fract = cons_fract
				self.cons_count = cons_count
				self.rank = rank
				self.rank_index = rank_index
				return

class Marker:
	def __init__(self):
		self.id = None

class Gene:
	def __init__(self):
		self.id = None
		self.marker = None
		self.contig = None
		self.annotations = []

class Annotation:
	def __init__(self):
		self.taxon = [None]*6
		self.score = None
	
	def add_taxon(self, genome_taxon, rank_index):
		for i in range(rank_index, 6):
			self.taxon[i] = genome_taxon.split(';')[i]

class Contig:
	def __init__(self):
		self.id = None
		self.genes = []
		self.genes_agree = 0
		self.genes_conflict = 0
		self.flagged = None
	
	def compare_taxonomy(self, bin):
		for gene in self.genes:
			agrees = []
			for annotation in gene.annotations:
				agrees.append(bin.cons_taxon == annotation.taxon[bin.rank_index])
			if any(agrees):
				self.genes_agree += 1
			else:
				self.genes_conflict += 1

	def flag(self, min_fract):
		if len(self.genes) == 0:
			self.flagged = None
		elif self.genes_conflict/float(len(self.genes)) >= min_fract:
			self.flagged = True
		else:
			self.flagged = False

def flag_contigs(db_dir, tmp_dir, args):
	
	# step 0. read in reference data files
	# cutoffs
	cutoffs = {}
	cutoffs_path = '%s/uscmg/max_fscores.tsv' % db_dir
	reader = csv.DictReader(open(cutoffs_path), delimiter="\t")
	for r in reader:
		key = (r['marker_id'], r['seq_type'], r['score_type'], r['taxlevel'])
		value = {'sensitive':r['cutoff_lower'], 'strict':r['cutoff_upper'], 'none':0.0}
		cutoffs[key] = value
	# taxonomy
	taxonomy = {}
	taxonomy_path = '%s/uscmg/genome_taxonomy.tsv' % db_dir
	reader = csv.DictReader(open(taxonomy_path), delimiter="\t")
	for r in reader:
		taxonomy[r['genome_id']] = r['taxonomy']
	# clustered seqs
	clusters = {}
	for type in ['ffn', 'faa']:
		clusters[type] = {}
		for file in os.listdir('%s/uscmg/%s' % (db_dir, type)):
			if file.split('.')[-1] == 'uc':
				with open('%s/uscmg/%s/%s' % (db_dir, type, file)) as f:
					for l in f:
						v = l.rstrip().split()
						rep_id = v[-1]
						seq_id = v[-2]
						if v[0] == 'S':
							clusters[type][seq_id] = [seq_id]
						elif v[0] == 'H':
							clusters[type][rep_id].append(seq_id)

	# step 1. determine if bin is archaea or bacteria; initialize domain-level markers
	# to do: normalize counts by site of marker gene sets
	marker_ids = set([])
	counts = {'bacteria':0, 'archaea':0}
	for aln_file in os.listdir('%s/alns' % tmp_dir):
		marker_id, seq_type, ext = aln_file.split(".")
		marker_ids.add(marker_id)
	for marker_id in marker_ids:
		if 'B' in marker_id:
			counts['bacteria'] += 1
		elif 'A' in marker_id:
			counts['archaea'] += 1
	domain = 'bacteria' if counts['bacteria'] >= counts['archaea'] else 'archaea'
	markers = {}
	for marker_id in marker_ids:
		if domain == 'bacteria' and 'B' in marker_id:
			markers[marker_id] = Marker()
			markers[marker_id].id = marker_id
			markers[marker_id].genes = []
		elif domain == 'archaea' and 'A' in marker_id:
			markers[marker_id] = Marker()
			markers[marker_id].id = marker_id
			markers[marker_id].genes = []

	# step 2. initialize marker genes found in bin
	bin = Bin()
	bin.genes = {}
	hmm_path = "%s/phyeco.hmmsearch" % tmp_dir
	for gene_id, aln in utility.fetch_hmm_best_hits(hmm_path).items():
		if aln['qacc'] not in markers:
			continue
		gene = Gene()
		gene.id = gene_id
		gene.contig = gene_id.rsplit('_', 1)[0]
		gene.marker = aln['qacc']
		bin.genes[gene_id] = gene
		markers[aln['qacc']].genes.append(gene)

	# annotate genes
	#    fetch all non-redundant taxonomic annotations for each gene
	seq_types = None
	if args['seq_type'] in ['both', 'either']:
		seq_types = ['ffn', 'faa']
	elif args['seq_type'] == 'protein':
		seq_types = ['faa']
	else:
		seq_types = ['ffn']

	for seq_type in seq_types:
		for marker_id in markers:
			aln_path = tmp_dir+"/alns/"+marker_id+"."+seq_type+".m8"
			for aln in utility.parse_blast(aln_path):
		
				# fetch all unique taxonomies for target sequence
				# a sequence can have multiple taxonomies if it was clustered with another sequence
				genome_taxa = []
				for target_id in clusters[seq_type][aln['tname']]:
					genome_id = target_id.split('_')[0]
					if genome_id not in taxonomy:
						continue
					elif taxonomy[genome_id] in genome_taxa:
						continue
					else:
						genome_taxa.append(taxonomy[genome_id])
				
				# loop over ranks; stop when gene has been annotated
				for rank_index, rank in enumerate(['s','g','f','o','c','p']):
					
					# decide to use ffn or faa at species level
					if args['seq_type'] == 'either':
						if seq_type == 'ffn' and rank != 's':
							continue
						elif seq_type == 'faa' and rank == 's':
							continue
					
					# get minimum % identity cutoff for transfering taxonomy
					#   if cutoff_type is None, indicates that no cutoff should be used
					min_pid = cutoffs[marker_id, seq_type, 'pid', rank][args['cutoff_type']]

					if float(aln['pid']) < float(min_pid):
						continue
					
					# add taxonomy
					for genome_taxon in genome_taxa:
						annotation = Annotation()
						annotation.add_taxon(genome_taxon, rank_index)
						annotation.score = float(aln['bitscore'])
						bin.genes[aln['qname']].annotations.append(annotation)
					
					# stop when gene has been annotated at lowest rank
					break

	# optionally remove annotations matching <exclude_clades>
	if args['exclude_clades'] is not None:
		bin.exclude_clades(args['exclude_clades'])

	# optionally take top hit only
	if args['hit_type'] == 'top_hit':
		bin.only_keep_top_hits()

	# create None annotations for unannotated genes
	for gene in bin.genes.values():
		if len(gene.annotations) == 0:
			gene.annotations.append(Annotation())

	# classify bin
	bin.classify_taxonomy(args['allow_noclass'], args['bin_fract'])
	
	# flag contigs with discrepant taxonomy
	bin.contigs = {}
	for id, seq in utility.parse_fasta(args['fna']):
		bin.contigs[id] = Contig()
		bin.contigs[id].id = id
		bin.contigs[id].length = len(seq)

	for gene in bin.genes.values():
		bin.contigs[gene.contig].genes.append(gene)

	if bin.cons_taxon is not None:
		for contig in bin.contigs.values():
			contig.compare_taxonomy(bin)
			contig.flag(args['contig_fract'])

	# write results
	flagged_contigs = []
	for contig in bin.contigs.values():
		if contig.flagged:
			flagged_contigs.append(contig.id)

	return flagged_contigs

def main():

	args = parse_args()
	utility.add_tmp_dir(args)
	utility.check_input(args)
	utility.check_dependencies(['prodigal', 'hmmsearch', 'blastp', 'blastn'])
	utility.check_database(args)
	
	print ("\n## Calling genes with Prodigal")
	utility.run_prodigal(args['fna'], args['tmp_dir'])

	print ("\n## Identifying PhyEco phylogenetic marker genes with HMMER")
	utility.run_hmmsearch(args['db'], args['tmp_dir'], args['tmp_dir'], args['threads'])
	extract_homologs(args['tmp_dir'])

	print ("\n## Performing pairwise BLAST alignment of marker genes against database")
	align_homologs(args['db'], args['tmp_dir'], args['seq_type'], args['threads'])

	print ("\n## Finding taxonomic outliers")
	flagged = flag_contigs(args['db'], args['tmp_dir'], args)
	out = '%s/flagged_contigs' % args['tmp_dir']
	print ("   flagged contigs: %s" % out)
	with open(out, 'w') as f:
		for contig in flagged:
			f.write(contig+'\n')







