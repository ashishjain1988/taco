'''
TACO: Transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2016
'''
import os
import collections
import subprocess
import logging

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2015"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.5"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def check_parallel_sort(num_processes):
    # check if parallel sort exists (GNU Coreutils 8.6+)
    parallel_sort = False
    with open(os.devnull, "w") as fnull:
        cmdline = 'echo "2 1" | %s --parallel=2'
        if subprocess.call(cmdline % 'gsort', stdout=fnull, stderr=fnull, shell=True) == 0:
            parallel_sort = 'gsort'
        if subprocess.call(cmdline % 'sort', stdout=fnull, stderr=fnull, shell=True) == 0:
            parallel_sort = 'sort'

    if not parallel_sort:
        logging.warning('Command line "sort" command does not support '
                        '--parallel flag. For improved performance, consider '
                        'upgrading/installing the latest GNU coreutils to '
                        'enable parallel sort.')
        args = ['sort']
    else:
        logging.debug('Command line "sort" supports --parallel flag')
        args = [parallel_sort, '--parallel=%d' % num_processes]
    return args


def sort_bed(input_file, output_file, num_processes=1, tmp_dir=None):
    if num_processes > 1:
        args = check_parallel_sort(num_processes)
    else:
        args = ['sort']
    if tmp_dir is not None:
        args.extend(["-T", tmp_dir])
    args.extend(["-k1,1", "-k2,2n", input_file])
    myenv = os.environ.copy()
    myenv["LC_ALL"] = "C"
    return subprocess.call(args, stdout=open(output_file, "w"), env=myenv)


def merge_bed(input_files, output_file, num_processes=1, tmp_dir=None):
    if num_processes > 1:
        args = check_parallel_sort()
    else:
        args = ['sort']
    if tmp_dir is not None:
        args.extend(['-T', tmp_dir])
    args.extend(['-k1,1', '-k2,2n', '-m'])
    args.extend(input_files)
    myenv = os.environ.copy()
    myenv["LC_ALL"] = "C"
    return subprocess.call(args, stdout=open(output_file, "w"), env=myenv)
