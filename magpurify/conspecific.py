#!/usr/bin/env python

import argparse
import os
import sys
from operator import itemgetter
from . import utility


def fetch_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=argparse.SUPPRESS,
        description="MAGpurify: conspecific module: identify contigs that fail to align to closely related genomes",
    )
    parser.add_argument('program', help=argparse.SUPPRESS)
    parser.add_argument('fna', type=str, help="""Path to input genome in FASTA format""")
    parser.add_argument(
        'out',
        type=str,
        help="""Output directory to store results and intermediate files""",
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help="""Number of CPUs to use (default=1)""",
    )
    parser.add_argument(
        '--mash-sketch',
        type=str,
        required=True,
        help="""Path to Mash sketch of reference genomes""",
    )
    parser.add_argument(
        '--mash-dist',
        type=float,
        default=0.05,
        help="Mash distance to reference genomes (default=0.05)",
    )
    parser.add_argument(
        '--max-genomes',
        type=int,
        default=25,
        help="Max number of genomes to use (default=25)",
    )
    parser.add_argument(
        '--min-genomes',
        type=int,
        default=1,
        help="Min number of genomes to use (default=1)",
    )
    parser.add_argument(
        '--contig-aln',
        type=float,
        default=0.50,
        help="Minimum fraction of contig aligned to reference (default=0.50)",
    )
    parser.add_argument(
        '--contig-pid',
        type=float,
        default=95.0,
        help="Minimum percent identity of contig aligned to reference (default=95.0)",
    )
    parser.add_argument(
        '--hit-rate',
        type=float,
        default=0.00,
        help="Hit rate for flagging contigs (default=0.00)",
    )
    parser.add_argument(
        '--exclude',
        default='',
        help="Comma-separated list of references to exclude",
    )
    args = vars(parser.parse_args())
    args['exclude'] = args['exclude'].split(',')
    return args


def run_mash(mash_sketch, fna_path, tmp_dir, threads=1):
    out_path = '%s/mash.dist' % tmp_dir
    command = "mash dist -p %s -d 0.25 %s %s > %s" % (
        threads,
        fna_path,
        mash_sketch,
        out_path,
    )
    out, err = utility.run_process(command)
    with open(tmp_dir + '/id_map.tsv', 'w') as f:
        for id, rec in enumerate(utility.parse_mash(out_path)):
            f.write(str(id) + '\t' + rec['target'] + '\n')


def find_conspecific(tmp_dir, max_dist, exclude):
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
    return targets


def blastn(query, target, outdir, id):
    out_path = outdir + '/' + id + '.m8'
    if not os.path.exists(out_path):
        cmd = "blastn -outfmt '6 std qlen slen' "
        cmd += "-max_target_seqs 1 -max_hsps 1 "
        cmd += "-query %s -subject %s " % (query, target)
        out, err = utility.run_process(cmd)
        out = out.decode("utf-8")
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
    for line in open(args['tmp_dir'] + '/id_map.tsv'):
        id, target = line.split()
        target_to_id[target] = id
    alignments = []
    for genome_path, mash_dist in genomes:
        blast_out = blastn(
            args['fna'], genome_path, args['tmp_dir'], target_to_id[target]
        )
        alignments.append(blast_out)
    return alignments


def find_contig_targets(args, genomes, alignments):
    contigs = dict(
        [
            (id, {'hits': 0, 'len': len(seq), 'genomes': []})
            for id, seq in utility.parse_fasta(args['fna'])
        ]
    )
    for genome, alns in zip(genomes, alignments):
        hits = id_blast_hits(alns, args['contig_aln'], args['contig_pid'])
        for contig in hits:
            contigs[contig]['hits'] += 1
            contigs[contig]['genomes'].append(genome[0])
    for id in contigs:
        hit_rate = contigs[id]['hits'] / float(len(alignments))
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
    print("\n## Finding conspecific genomes in database")
    run_mash(args['mash_sketch'], args['fna'], args['tmp_dir'], args['threads'])
    genomes = find_conspecific(args['tmp_dir'], args['mash_dist'], args['exclude'])
    print("   %s genomes within %s mash-dist" % (len(genomes), args['mash_dist']))
    out = '%s/conspecific.list' % args['tmp_dir']
    with open(out, 'w') as f:
        f.write('genome_id\tmash_dist\n')
        for genome_id, mash_dist in genomes:
            f.write(genome_id + '\t' + str(mash_dist) + '\n')
    print("   list of genomes: %s" % (out))
    print("   mash output: %s/mash.dist" % args['tmp_dir'])
    if len(genomes) < args['min_genomes']:
        sys.exit("\nError: insufficient number of conspecific genomes\n")
    if len(genomes) > args['max_genomes']:
        print("\n## Selecting top %s most-similar genomes" % args['max_genomes'])
        genomes = genomes[0 : args['max_genomes']]
        out = '%s/conspecific_subset.list' % args['tmp_dir']
        with open(out, 'w') as f:
            f.write('genome_id\tmash_dist\n')
            for genome_id, mash_dist in genomes:
                f.write(genome_id + '\t' + str(mash_dist) + '\n')
        print("   list of genomes: %s" % (out))
    print("\n## Performing pairwise alignment of contigs in bin to database genomes")
    alignments = align_contigs(args, genomes)
    num_alns = sum(len(_.split('\n')) for _ in alignments)
    print("   total alignments: %s" % num_alns)
    print("\n## Summarizing alignments")
    contigs = find_contig_targets(args, genomes, alignments)
    out = '%s/contig_hits.tsv' % args['tmp_dir']
    with open(out, 'w') as f:
        f.write('contig_id\tlength\talignment_rate\n')
        for contig, values in contigs.items():
            row = [contig, str(values['len']), '%s/%s' % (values['hits'], len(genomes))]
            f.write('\t'.join(row) + '\n')
    print("   contig features: %s" % out)
    print("\n## Identifying contigs with no conspecific alignments")
    flagged = flag_contigs(args, contigs)
    out = f"{args['tmp_dir']}/flagged_contigs"
    with open(out, 'w') as f:
        for contig in flagged_contigs:
            f.write(contig + '\n')
    print(f"   {len(flagged)} flagged contigs: {out}")
