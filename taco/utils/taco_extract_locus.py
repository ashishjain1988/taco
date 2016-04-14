#!/usr/bin/env python
'''
TACO: Multi-sample transcriptome assembly from RNA-Seq

Utility script that extracts the GTF features for an individual locus.
'''
import os
import sys
import argparse
import logging
import subprocess

from taco.lib.run import Results
from taco.lib.assemble import parse_locus_index


__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.3"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument('taco_dir',
                        help='path to taco run output directory')
    parser.add_argument('locus_id',
                        help='locus id to fetch')
    args = parser.parse_args()

    output_dir = args.taco_dir
    locus_id = args.locus_id
    if not os.path.exists(output_dir):
        parser.error('TACO output directory %s not found' % output_dir)
    if not os.path.isdir(output_dir):
        parser.error('%s is not a directory' % output_dir)
    r = Results(output_dir)

    logging.info('output directory: %s' % output_dir)
    logging.info('locus id: %s' % locus_id)

    # search for locus
    for locus in parse_locus_index(r.locus_index_file):
        if locus.name == locus_id:
            break

    logging.info('locus: %s features: %d coords: %s:%d-%d' %
                 (locus.name, locus.num_lines, locus.chrom,
                  locus.start, locus.end))

    # extract locus
    with open(r.transfrags_bed_file) as fh:
        fh.seek(locus.filepos)
        for i in xrange(locus.num_lines):
            print fh.next(),

    logging.info('done.')

if __name__ == '__main__':
    sys.exit(main())
