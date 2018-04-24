
import sys, os

def check_input(args):
	if not os.path.exists(args['fna_path']):
		error = "Input file not found: %s" % args['fna_path']
		sys.exit(error)
	
def check_output(args):
	if not os.path.isdir(args['out_dir']):
		os.makedirs(args['out_dir'])

def check_dependencies(programs):
	for program in programs:
		if not exists_on_env_path(program):
			error = "\nRequired program '%s' not found\n" % program
			error += "Make sure this program has been installed and added to your PATH"
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
	if args['db'] is None:
		error = "\nError: No reference database specified\n"
		error += "Use the flag -d to specify a database,\n"
		error += "Or set the IMAGEN_DB environmental variable: export IMAGEN_DB=/path/to/midas/db\n"
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

def parse_blast(fpath):
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
	lines = open(fpath).read().rstrip('\n').split('\n')
	for line in lines:
		values = line.split('\t')
		record = dict([(f[0], f[1](v)) for f,v in zip(formats, values)])
		record['qcov'] = 100*record['aln']/float(record['qlen'])
		record['tcov'] = 100*record['aln']/float(record['tlen'])
		yield record


