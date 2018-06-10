#!/usr/bin/env python
'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
'''
import sys
import logging

from taco.lib.run import Run

__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2018"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "MIT"
__version__ = "0.7.3"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"

PROFILE = False
PROFILE_FILENAME = 'taco.prof'


def main():
    # instantiate from command line
    R = Run.create()
    R.log_args()
    logging.info('Samples: %d' % (len(R.samples)))
    # merge gtf files
    msg = 'Aggregating GTF files'
    if R.status.aggregate:
        logging.info('[SKIPPING] %s' % msg)
    else:
        logging.info(msg)
        R.aggregate()
    # index loci
    msg = 'Indexing Loci'
    if R.status.index_loci:
        logging.info('[SKIPPING] %s' % msg)
    else:
        logging.info(msg)
        R.index_loci()
    # assemble
    msg = 'Assembling GTF files'
    if R.status.assemble:
        logging.info('[SKIPPING] %s' % msg)
    else:
        logging.info(msg)
        R.assemble()
    return 0


if __name__ == '__main__':
    if PROFILE:
        import cProfile
        sys.exit(cProfile.run('main()', filename=PROFILE_FILENAME))
    else:
        return_value = main()
        if (return_value != 0):
            sys.exit(return_value)
