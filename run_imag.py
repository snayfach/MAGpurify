#!/usr/bin/env python

import sys
import imag
from imag.utility import *

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
	args = vars(parser.parse_args())
	if (not args['db']
			and 'IMAGEN_DB' in os.environ):
		args['db'] = os.environ['IMAGEN_DB']
	return args

if __name__ == "__main__":
	
	args = fetch_args()

	check_input(args)
	check_output(args)
	check_dependencies(['blastn'])
	check_database(args)

	contigs = {}
	
	sys.stdout.write("discordant taxonomy based on universal-single-copy marker-genes\n")
	from imag import uscmg
	contigs['uscmg'] = uscmg.main(args)
	sys.stdout.write("   %s contigs flagged\n" % len(contigs['uscmg']))
	
	sys.stdout.write("discordant taxonomy based on clade-specific marker-genes\n")
	from imag import csmg
	contigs['csmg'] = csmg.main(args)
	sys.stdout.write("   %s contigs flagged\n" % len(contigs['csmg']))

	sys.stdout.write("hits to database of contaminants\n")
	from imag import contam
	contigs['contam'] = contam.main(args)
	sys.stdout.write("   %s contigs flagged\n" % len(contigs['contam']))
	
	sys.stdout.write("outlier tetranucleotide frequency\n")
	from imag import tetra
	contigs['tetra'] = tetra.main(args)
	sys.stdout.write("   %s contigs flagged\n" % len(contigs['tetra']))

	sys.stdout.write("outlier gc content\n")
	from imag import gc
	contigs['gc_content'] = gc.main(args)
	sys.stdout.write("   %s contigs flagged\n" % len(contigs['gc_content']))








