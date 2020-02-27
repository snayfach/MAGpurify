# -*- coding: utf-8 -*-
#
#   This file is part of the magpurify package, available at:
#   https://github.com/snayfach/MAGpurify
#
#   Magpurify is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <https://www.gnu.org/licenses/>.

import os
import subprocess as sp
import sys
from Bio import Seq, SeqIO


def add_tmp_dir(args):
    tmp_dir = f"{args['out']}/{args['program']}"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    args["tmp_dir"] = tmp_dir


def check_input(args):
    if not os.path.exists(args["fna"]):
        error = f"\nInput file not found: {args['fna']}\n"
        sys.exit(error)


def check_dependencies(programs):
    for program in programs:
        if not exists_on_env_path(program):
            error = f"\nRequired program '{program}' not found\n"
            error += "Make sure this program has been installed and added to your PATH\n"
            sys.exit(error)


def exists_on_env_path(program):
    " Check whether program exists in PATH and is executable"
    for dir in os.environ["PATH"].split(os.pathsep):
        fpath = dir + "/" + program
        if os.path.exists(fpath) and os.access(fpath, os.X_OK):
            return True
    return False


def check_database(args):
    if args["db"] is None and "MAGPURIFYDB" in os.environ:
        args["db"] = os.environ["MAGPURIFYDB"]
    if args["db"] is None:
        error = "\nError: No reference database specified\n"
        error += "Use the --db argument to specify a database,\n"
        error += "Or set the MAGPURIFYDB environmental variable: export MAGPURIFYDB=/path/to/MAGpurify_db_v1.0.0\n"
        sys.exit(error)
    if not os.path.isdir(args["db"]):
        error = f"\nError: Specified reference database does not exist: {args['db']}\n"
        error += (
            "\nCheck that you've entered the path correctly and the database exists\n"
        )
        sys.exit(error)


def reverse_complement(sequence):
    rc = Seq.Seq(sequence).reverse_complement()
    return str(rc)


def run_process(command):
    process = sp.Popen(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    out, err = process.communicate()
    if process.returncode != 0:
        err_message = (
            f"\nError encountered executing:\n{command}\n\nError message:\n{err}\n"
        )
        sys.exit(err_message)
    return out, err


def parse_last(inpath):
    fields = [
        "qid",
        "tid",
        "pid",
        "aln",
        "mis",
        "gap",
        "qstart",
        "qend",
        "tstart",
        "tend",
        "eval",
        "score",
        "qlen",
        "tlen",
    ]
    for line in open(inpath):
        if line[0] == "#":
            continue
        else:
            values = line.rstrip().split()
            d = dict([(f, v) for f, v in zip(fields, values)])
            d["qcov"] = float(d["aln"]) / float(d["qlen"])
            d["tcov"] = float(d["aln"]) / float(d["tlen"])
            yield d


def parse_blast(input, type="file"):
    formats = [
        ("qname", str),
        ("tname", str),
        ("pid", float),
        ("aln", int),
        ("mis", int),
        ("gap", int),
        ("qstart", int),
        ("qend", int),
        ("tstart", int),
        ("tend", int),
        ("evalue", float),
        ("bitscore", float),
        ("qlen", int),
        ("tlen", int),
    ]
    if type == "file":
        lines = open(input).read().rstrip("\n").split("\n")
    else:
        lines = input.rstrip("\n").split("\n")
    if lines == [""]:
        return
    for line in lines:
        values = line.split("\t")
        record = dict([(f[0], f[1](v)) for f, v in zip(formats, values)])
        record["qcov"] = 100 * record["aln"] / float(record["qlen"])
        record["tcov"] = 100 * record["aln"] / float(record["tlen"])
        yield record


def parse_mash(fpath):
    fields = ["query", "target", "dist", "pvalue", "fraction"]
    formats = [str, str, float, float, str]
    out = open(fpath).read()
    if len(out) > 0:
        lines = out.rstrip("\n").split("\n")
        for line in lines:
            values = line.split()
            rec = dict([(f, m(v)) for f, m, v in zip(fields, formats, values)])
            yield rec


def parse_hmmsearch(fpath):
    with open(fpath) as infile:
        fields = [
            "tname",
            "tacc",
            "tlen",
            "qname",
            "qacc",
            "qlen",
            "evalue",
            "score",
            "bias",
            "ndom",
            "tdom",
            "c-evalue",
            "i-evalue",
            "domscore",
            "dombias",
            "hmmfrom",
            "hmmto",
            "alifrom",
            "alito",
            "envfrom",
            "envto",
            "prob",
            "tdesc",
        ]
        formts = [
            str,
            str,
            int,
            str,
            str,
            int,
            float,
            float,
            float,
            int,
            int,
            float,
            float,
            float,
            float,
            int,
            int,
            int,
            int,
            int,
            int,
            float,
            str,
        ]
        for line in infile:
            if line[0] == "#":
                continue
            values = line.rstrip("\n").split(None, 22)
            r = dict(
                [
                    (field, format(value))
                    for field, format, value in zip(fields, formts, values)
                ]
            )
            r["tcov"] = float(r["alito"] - r["alifrom"] + 1) / r["tlen"]  # target is gene
            r["qcov"] = float(r["hmmto"] - r["hmmfrom"] + 1) / r["qlen"]  # query is hmm
            yield r


def fetch_hmm_best_hits(fpath):
    gene_to_aln = {}
    for aln in parse_hmmsearch(fpath):
        if aln["tname"] not in gene_to_aln:
            gene_to_aln[aln["tname"]] = aln
        elif aln["score"] > gene_to_aln[aln["tname"]]["score"]:
            gene_to_aln[aln["tname"]] = aln
    return gene_to_aln


def run_coverm(bam_list, out_dir):
    cmd = "coverm contig "
    cmd += "--contig-end-exclusion 75 "
    cmd += "--min-read-percent-identity 0.97 "
    cmd += f"--bam-files {' '.join(bam_list)} "
    cmd += f"> {out_dir}/coverage.tsv"
    out, err = run_process(cmd)


def run_prodigal(fna_path, out_dir):
    cmd = "prodigal "
    cmd += f"-i {fna_path} "  # input fna
    cmd += f"-a {out_dir}/genes.faa "  # protein seqs
    cmd += f"-d {out_dir}/genes.ffn "  # nucleotide seqs
    cmd += f"-o {out_dir}/genes.out "  # prodigal output
    cmd += "> /dev/null"  # prodigal output
    out, err = run_process(cmd)


def run_lastal(db_dir, out_dir, threads=1, seed_freq=10):
    cmd = "lastal -p BLOSUM62 -P 1 -f blasttab+ "
    cmd += f"-m {seed_freq} "
    cmd += f"{db_dir}/clade-markers/markers.faa "
    cmd += f"{out_dir}/genes.faa "
    cmd += f"-P {threads} "
    cmd += f"> {out_dir}/genes.m8 "
    out, err = run_process(cmd)


def run_hmmsearch(db_dir, in_path, out_dir, threads=1):
    cmd = "hmmsearch "
    cmd += "--noali "
    cmd += f"--domtblout {out_dir}/phyeco.hmmsearch "
    cmd += f"--cpu {threads} "
    cmd += "--cut_ga "
    cmd += f"{db_dir}/phylo-markers/PhyEco.hmm "
    cmd += f"{in_path}/genes.faa "
    out, err = run_process(cmd)


def run_blastp(db_path, query_path, out_path, threads=1, max_targets=1, qcov=40):
    cmd = "blastp "
    cmd += f"-db {db_path} "
    cmd += f"-query {query_path} "
    cmd += f"-out {out_path} "
    cmd += "-outfmt '6 std qlen slen' "
    cmd += f"-num_threads {threads} "
    cmd += f"-max_target_seqs {max_targets} "
    cmd += f"-qcov_hsp_perc {qcov} "
    cmd += "-max_hsps 1 "
    out, err = run_process(cmd)


def run_blastn(db_path, query_path, out_path, threads=1, max_targets=1, qcov=40):
    cmd = "blastn -task blastn "
    cmd += f"-db {db_path} "
    cmd += f"-query {query_path} "
    cmd += f"-out {out_path} "
    cmd += "-outfmt '6 std qlen slen' "
    cmd += f"-num_threads {threads} "
    cmd += f"-max_target_seqs {max_targets} "
    cmd += f"-qcov_hsp_perc {qcov} "
    cmd += "-max_hsps 1 "
    out, err = run_process(cmd)


def parse_fasta(path):
    with open(path) as file:
        for record in SeqIO.parse(file, "fasta"):
            id = record.id
            seq = str(record.seq).upper()
            yield id, seq
