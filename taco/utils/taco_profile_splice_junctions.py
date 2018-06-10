#!/usr/bin/env python
'''
TACO: Multi-sample transcriptome assembly from RNA-Seq

Utility script that profiles the splice junctions in a BED file

Requirements: pysam library
'''
import os
import sys
import argparse
import logging
import operator
from collections import Counter

from taco.lib.base import Strand
from taco.lib.transfrag import Transfrag

from taco.lib.pysam.cfaidx import FastaFile


__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2018"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "MIT"
__version__ = "0.7.3"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


rev_comp_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N': 'N'}


def dna_reverse_complement(seq):
    return ''.join(rev_comp_dict[x] for x in reversed(seq))


def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument('genome_fasta_file')
    parser.add_argument('bed_file')
    args = parser.parse_args()

    # check args
    if not os.path.exists(args.genome_fasta_file):
        parser.error('genome fasta file %s not found' % args.genome_fasta_file)
    if not os.path.exists(args.bed_file):
        parser.error('bed file %s not found' % args.bed_file)
    logging.info('genome fasta file: %s' % args.genome_fasta_file)
    logging.info('bed file: %s' % args.bed_file)

    # process bed file to get junctions
    logging.info('Reading Junctions')
    splice_juncs = set()
    fasta_fh = FastaFile(args.genome_fasta_file)
    with open(args.bed_file) as bed_fh:
        for line in bed_fh:
            t = Transfrag.from_bed(line)
            if t.chrom not in fasta_fh:
                continue
            for start, end in t.iterintrons():
                splice_juncs.add((t.chrom, start, end, t.strand))
    logging.info('Read %d Junctions' % (len(splice_juncs)))

    logging.info('Profiling Splice Motifs')
    motif_counter = Counter()
    for chrom, start, end, strand in splice_juncs:
        s = fasta_fh.fetch(chrom, start, start + 2)
        s += fasta_fh.fetch(chrom, end - 2, end)
        if strand == Strand.NEG:
            s = dna_reverse_complement(s)
        motif_counter[s] += 1
    fasta_fh.close()

    # report statistics
    total = sum(motif_counter.values())
    print '\t'.join(['motif', 'count', 'frac'])
    for motif, count in motif_counter.most_common():
        print '\t'.join([motif, str(count), str(float(count) / total)])
    logging.info('Done')


if __name__ == '__main__':
    sys.exit(main())
