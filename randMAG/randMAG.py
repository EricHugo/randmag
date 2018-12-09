#!/usr/bin/env python

"""
randMAG is a script designed to split whole genomes into contigs
based on a distribution of sizes. Randomly reduce the completeness
based on given percentage, as well as simulate contamination by
introducing contigs from other simulated MAGs.

Produces table of created MAGs with % completeness, contamination,
and source genome.
"""

import argparse
import sys
import re
import os
import random
import string
import numpy as np
from Bio import SeqIO
from scipy import stats
from fitdist import FitDist

def _worker(seq, dist_param, min_length, comp=1.0, contam=0.0):
    name = os.path.basename(seq)
    print(name)
    # find full length of genome
    genome, genome_length = get_seq_length(seq)
    # split into contig given sizes in loop
    contigs = {len(contig): contig for contig in split_contigs(genome, dist_param)}
    #print(contigs)
    # remove contigs based on desired completeness
    # return dict of genome with list of contig sizes produced
    ## to be used for contamination simulation
    #print(genome_length)
    return contigs, name.split('.')[0]

def get_seq_length(seq):
    seq_fasta = [seqi.seq for seqi in SeqIO.parse(seq, "fasta")]
    if len(seq_fasta) > 1:
        print("Warning, more than one contig in genome", file=sys.stderr)
        print(len(seq_fasta))
    # only takes the length of first contig, should be chromosome
    seq_length = len(seq_fasta[0])
    #print(type(seq_fasta[0]))
    #print(seq_fasta[0][0:20])
    return seq_fasta, seq_length

def split_contigs(seq, params, min_length=300):
    # get the parameters here
    shape, _, scale = params
    while seq[0]:
        size = round(np.random.gamma(shape, scale))
        #print(size)
        contig = seq[0][0:size-1]
        seq[0] = seq[0][size-1:]
        if len(seq[0]) < min_length:
            contig = contig + seq[0]
            seq[0] = ''
        yield contig

def output_randcontigs(name, genome):
    """Writes shuffled contigs to file"""
    randname = ''.join([random.choice(string.ascii_lowercase)
                        for _ in range(8)])
    with open(name + '.' + randname + ".fna", "w") as handle:
        # necessary to put keys into list before shuffle
        shuffled_headers = list(genome.keys())
        random.shuffle(shuffled_headers)
        for header in shuffled_headers:
            handle.write(">" + str(header) + "\n")
            handle.write(str(genome[header]) + "\n")

def main():
    parser = argparse.ArgumentParser(description="""
                                    Simulate MAGs by splitting whole
                                    genomes."""
                                    )
    parser.add_argument("genome_tab", help="Genomes in .fna in list format")
    parser.add_argument("distribution", help="Set of lengths to base the "\
                                             "distribution on.")
    parser.add_argument("-c", "--completeness", default="1.0-1.0",
                        help="Range of completeness levels to be produced. "\
                             "Default=1.0-1.0")
    parser.add_argument("-r", "--contamination", default="0.0-0.0",
                        help="Range of contamination to be included in "\
                             "produced MAGs. Default=0.0-0.0")
    parser.add_argument("-n", "--num", default=False, type=int,
                        help="Number of randomised MAGs to produce. Produces "\
                             "one randomised per provided genome by default.")
    args = parser.parse_args()

    with open(args.genome_tab) as seq_file:
        input_seqs = [seq.strip() for seq in seq_file
                      if not re.match('#|\n', seq)]

    with open(args.distribution) as dist_file:
        distances = [int(dist.strip()) for dist in dist_file
                    if not re.match('#|\n', dist)]
    # here get distribution and parameters
    f = FitDist(distances)
    distribution = f.find_best_dist()
    f.histogram(distribution)
    params = stats.gamma.fit(distances)
    #print(params)
    # feed into worker with random sampling

    if not args.num:
        args.num = len(input_seqs)
    while args.num > 0:
        for seq in input_seqs:
            contigs, name = _worker(seq, params, min(distances))
            output_randcontigs(name, contigs)
            args.num -= 1
            #print(args.num)
            if args.num <= 0:
                break

if __name__ == "__main__":
    main()
