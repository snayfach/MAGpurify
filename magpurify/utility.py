
import sys, os, Bio.SeqIO

def add_tmp_dir(args):
	tmp_dir = '%s/%s' % (args['out'], args['program'])
	if not os.path.exists(tmp_dir):
		os.makedirs(tmp_dir)
	args['tmp_dir'] = tmp_dir

def check_input(args):
	if not os.path.exists(args['fna']):
		error = "\nInput file not found: %s\n" % args['fna']
		sys.exit(error)
	
def check_dependencies(programs):
	for program in programs:
		if not exists_on_env_path(program):
			error = "\nRequired program '%s' not found\n" % program
			error += "Make sure this program has been installed and added to your PATH\n"
			sys.exit(error)

def exists_on_env_path(program):
	""" Check whether program exists in PATH and is executable """
	for dir in os.environ["PATH"].split(os.pathsep):
		fpath = dir+'/'+program
		if (os.path.exists(fpath) and
				os.access(fpath, os.X_OK)):
			return True
	return False

def check_database(args):
	if args['db'] is None and 'MAGPURIFYDB' in os.environ:
		args['db'] = os.environ['MAGPURIFYDB']
	if args['db'] is None:
		error = "\nError: No reference database specified\n"
		error += "Use the flag --db to specify a database,\n"
		error += "Or set the MAGPURIFYDB environmental variable: export MAGPURIFYDB=/path/to/MAGpurify_db_v1.0.0\n"
		sys.exit(error)
	if not os.path.isdir(args['db']):
		error = "\nError: Specified reference database does not exist: %s\n" % args['db']
		error += "\nCheck that you've entered the path correctly and the database exists\n"
		sys.exit(error)

def run_process(command):
	""" Capture stdout, stderr. Check unix exit code and exit if non-zero """
	import subprocess as sp
	process = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
	out, err = process.communicate()
	if process.returncode != 0:
		err_message = "\nError encountered executing:\n%s\n\nError message:\n%s\n" % (command, err)
		sys.exit(err_message)
	return out, err

def parse_fasta(fpath):
	return Bio.SeqIO.parse(fpath, "fasta")

def parse_last(inpath):
	fields = ['qid', 'tid', 'pid', 'aln', 'mis', 'gap', 'qstart', 'qend',
			  'tstart', 'tend', 'eval', 'score', 'qlen', 'tlen']
	for line in open(inpath):
		if line[0] == '#':
			continue
		else:
			values = line.rstrip().split()
			d = dict([(f,v) for f,v in zip(fields, values)])
			d['qcov'] = float(d['aln'])/float(d['qlen'])
			d['tcov'] = float(d['aln'])/float(d['tlen'])
			yield d
			
def parse_blast(input, type='file'):
	formats = [('qname', str),
			   ('tname', str),
			   ('pid', float),
			   ('aln', int),
			   ('mis', int),
			   ('gap', int),
			   ('qstart', int),
			   ('qend', int),
			   ('tstart', int),
			   ('tend', int),
			   ('evalue', float),
			   ('bitscore', float),
		  	   ('qlen', int),
			   ('tlen', int)]
	if type=='file':
		lines = open(input).read().rstrip('\n').split('\n')
	else:
		lines = input.rstrip('\n').split('\n')
	if lines == ['']:
		return
	for line in lines:
		values = line.split('\t')
		record = dict([(f[0], f[1](v)) for f,v in zip(formats, values)])
		record['qcov'] = 100*record['aln']/float(record['qlen'])
		record['tcov'] = 100*record['aln']/float(record['tlen'])
		yield record

def parse_mash(fpath):
	fields = ['query', 'target', 'dist', 'pvalue', 'fraction']
	formats = [str, str, float, float, str]
	out = open(fpath).read()
	if len(out) > 0:
		lines = out.rstrip('\n').split('\n')
		for line in lines:
			values = line.split()
			rec = dict([(f,m(v)) for f,m,v in zip(fields, formats, values)])
			yield rec

def parse_hmmsearch(fpath):
	with open(fpath) as infile:
		fields = ['tname','tacc','tlen','qname','qacc','qlen','evalue','score',
				  'bias','ndom','tdom','c-evalue','i-evalue','domscore','dombias',
				  'hmmfrom','hmmto','alifrom','alito','envfrom','envto','prob','tdesc']
		formts = [str,str,int,str,str,int,float,float,float,int,int,float,
				  float,float,float,int,int,int,int,int,int,float,str]
		for line in infile:
			if line[0] == '#':
				continue
			values = line.rstrip('\n').split(None, 22)
			r = dict([(field, format(value)) for field, format, value in zip(fields, formts, values)])
			r['tcov'] = float(r['alito'] - r['alifrom'] + 1)/r['tlen'] # target is gene
			r['qcov'] = float(r['hmmto'] - r['hmmfrom'] + 1)/r['qlen'] # query is hmm
			yield r

def fetch_hmm_best_hits(fpath):
	gene_to_aln = {}
	for aln in parse_hmmsearch(fpath):
		if aln['tname'] not in gene_to_aln:
			gene_to_aln[aln['tname']] = aln
		elif aln['score'] > gene_to_aln[aln['tname']]['score']:
			gene_to_aln[aln['tname']] = aln
	return gene_to_aln

def run_prodigal(fna_path, out_dir):
	cmd = "prodigal "
	cmd += "-i %s " % fna_path # input fna
	cmd += "-a %s/genes.faa " % out_dir # protein seqs
	cmd += "-d %s/genes.ffn " % out_dir # nucleotide seqs
	cmd += "-o %s/genes.out " % out_dir # prodigal output
	cmd += "> /dev/null" # prodigal output
	out, err = run_process(cmd)

def run_lastal(db_dir, out_dir, threads=1, seed_freq=10):
	cmd = "lastal -p BLOSUM62 -P 1 -f blasttab+ "
	cmd += "-m %s " % seed_freq
	cmd += "%s/clade-markers/markers.faa " % db_dir
	cmd += "%s/genes.faa " % out_dir
	cmd += "-P %s " % threads
	cmd += "> %s/genes.m8 " % out_dir
	out, err = run_process(cmd)

def run_hmmsearch(db_dir, in_path, out_dir, threads=1):
	cmd = "hmmsearch "
	cmd += "--noali "
	cmd += "--domtblout %s/phyeco.hmmsearch " % out_dir
	cmd += "--cpu %s " % threads
	cmd += "--cut_ga "
	cmd += "%s/phylo-markers/PhyEco.hmm " % db_dir
	cmd += "%s/genes.faa " % in_path
	out, err = run_process(cmd)

def run_blastp(db_path, query_path, out_path, threads=1, max_targets=1, qcov=40):
	cmd = "blastp "
	cmd += "-db %s " % db_path
	cmd += "-query %s " % query_path
	cmd += "-out %s " % out_path
	cmd += "-outfmt '6 std qlen slen' "
	cmd += "-num_threads %s " % threads
	cmd += "-max_target_seqs %s " % max_targets
	cmd += "-qcov_hsp_perc %s " % qcov
	cmd += "-max_hsps 1 "
	out, err = run_process(cmd)

def run_blastn(db_path, query_path, out_path, threads=1, max_targets=1, qcov=40):
	cmd = "blastn -task blastn "
	cmd += "-db %s " % db_path
	cmd += "-query %s " % query_path
	cmd += "-out %s " % out_path
	cmd += "-outfmt '6 std qlen slen' "
	cmd += "-num_threads %s " % threads
	cmd += "-max_target_seqs %s " % max_targets
	cmd += "-qcov_hsp_perc %s " % qcov
	cmd += "-max_hsps 1 "
	out, err = run_process(cmd)

def parse_fasta(path):
	with open(path) as file:
		try: id = next(file).split()[0].lstrip('>')
		except: return
		seq = ''
		for line in file:
			if line[0]=='>':
				yield id, seq
				try: id = line.split()[0].lstrip('>')
				except: return
				seq = ''
			else:
				seq += line.rstrip()
		yield id, seq


