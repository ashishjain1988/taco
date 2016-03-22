'''
TACO: Transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer
'''
import os
import logging
import collections
import multiprocessing
import re

from base import Sample, TacoError
from gtf import GTF, GTFError
from stats import scoreatpercentile
from caggregate import caggregate
from csort import csort

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "mkiyer@umich.edu"
__status__ = "Development"

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [atoi(c) for c in re.split('(\d+)', text)]

def aggregate(samples, ref_gtf_file, gtf_expr_attr, tmp_dir,
              output_gtf_file, stats_file):
    '''
    Aggregate/merge individual sample GTF files
    '''

    num_of_files = 0

    # aggregate ref gtf
    '''
    if ref_gtf_file is not None:
        tmp_ref_gtf_file = "transcripts." + str(num_of_files) + ".gtf"
        num_of_files += 1
        tmp_ref_file = os.path.join(tmp_dir, tmp_ref_gtf_file)
        logging.debug('Reference: %s' % ref_gtf_file)
        caggregate(ref_gtf_file, str(Sample.REF_ID), gtf_expr_attr,
                   tmp_ref_file, str(True))
    '''

    # aggregate sample gtfs
    input_filenames = [sample.gtf_file for sample in samples]
    total_cpu_count = multiprocessing.cpu_count()

    caggregate(input_filenames, gtf_expr_attr, tmp_dir, str(False), total_cpu_count)
    sys.exit(1)

    for sample in samples:
        # Create new GTF transcript file names to prevent duplicates
        tmp_gtf_file = "transcripts." + str(num_of_files) + ".gtf"
        num_of_files += 1
        tmp_gtf_fullpath = os.path.join(tmp_dir, tmp_gtf_file)

        logging.debug('Sample: %s %s' % (sample._id, sample.gtf_file))
        caggregate(sample.gtf_file, str(sample._id), gtf_expr_attr,
                   tmp_gtf_fullpath, stats_file, str(False))
    sys.exit(1)
    print("NEED TO TAKE CARE OF REF_GTF_FILE")
    print("Need TO MAKE GOOD CPU SPLITS FOR AGGREGATE VS SORT")
    sys.exit(1)

    # merge sorted gtfs
    logging.info("Sorting GTF")
    aggregated_gtf_files = sorted([os.path.join(tmp_dir, s) for s in os.listdir(tmp_dir) if (".DS" not in s)], key=natural_keys)
    csort(aggregated_gtf_files, tmp_dir, multiprocessing.cpu_count())
