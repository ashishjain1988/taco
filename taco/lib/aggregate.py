'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer
'''
import os
import logging
import collections
import shutil
from multiprocessing import Process, JoinableQueue

from taco.lib.pysam.cfaidx import FastaFile
from bed import sort_bed, merge_bed
from batch_sort import merge_files
from base import Exon, Strand, Sample, Results, TacoError
from transfrag import Transfrag
from gtf import GTF, GTFError
from caggregate import caggregate


__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.5"
__maintainer__ = "Yashar Niknafs"
__email__ = "mkiyer@umich.edu"
__status__ = "Development"


DNA_COMPLEMENT_DICT = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N': 'N'}
SPLICE_MOTIFS_ALLOWED = {'GTAG', 'GCAG'}

def dna_reverse_complement(seq):
    return ''.join(DNA_COMPLEMENT_DICT[x] for x in reversed(seq.upper()))

def is_splice_motif_allowed(genome_fasta_fh, chrom_str, start, end, strand):
    s = genome_fasta_fh.fetch(chrom_str, start, start + 2)
    s += genome_fasta_fh.fetch(chrom_str, end - 2, end)

    if (strand == "-"):
        s = dna_reverse_complement(s)

    if s not in SPLICE_MOTIFS_ALLOWED:
        return False
    else:
        return True

def aggregate_sample(sample, gtf_expr_attr, is_ref, min_length, min_expr,
                     filter_splice_juncs, genome_fasta_fh,
                     bed_file_name, filtered_bed_file_name):
    logging.debug('Aggregate sample %s: %s' % (sample._id, sample.gtf_file))
    caggregate(sample.gtf_file, str(sample._id), gtf_expr_attr, str(is_ref), bed_file_name, filtered_bed_file_name, float(min_length), float(min_expr), str(filter_splice_juncs), genome_fasta_fh, is_splice_motif_allowed)

def aggregate_worker(input_queue, args, output_dir):
    results = Results(output_dir)
    # create temp directories
    tmp_dir = os.path.join(results.output_dir, 'tmp')
    if not os.path.exists(tmp_dir):
        logging.debug('\tcreating tmp dir %s' % (tmp_dir))
        os.makedirs(tmp_dir)
    # create set of unsorted results
    tmp_results = Results(tmp_dir)
    # setup genome fasta file
    genome_fasta_fh = None
    if args.filter_splice_juncs and args.ref_genome_fasta_file:
        genome_fasta_fh = FastaFile(args.ref_genome_fasta_file)
    # process samples via input queue
    while True:
        sample = input_queue.get()
        if sample is None:
            break
        aggregate_sample(sample,
                         gtf_expr_attr=args.gtf_expr_attr,
                         is_ref=(sample._id == Sample.REF_ID),
                         min_length=args.filter_min_length,
                         min_expr=args.filter_min_expr,
                         filter_splice_juncs=args.filter_splice_juncs,
                         genome_fasta_fh=genome_fasta_fh,
                         bed_file_name=tmp_results.transfrags_bed_file,
                         filtered_bed_file_name=tmp_results.transfrags_filtered_bed_file
                         )
        input_queue.task_done()
    input_queue.task_done()
    # cleanup and close files
    if genome_fasta_fh:
        genome_fasta_fh.close()

    # sort output files
    logging.debug('Sorting aggregated files: "%s"' % (output_dir))
    # sort bed file
    logging.debug('\ttransfrags bed file: %s' % (results.output_dir))
    retcode = sort_bed(tmp_results.transfrags_bed_file,
                       results.transfrags_bed_file,
                       num_processes=1,
                       tmp_dir=tmp_dir)
    if retcode != 0:
        raise TacoError('Error running linux sort')
    os.remove(tmp_results.transfrags_bed_file)
    # sort filtered bed file
    logging.debug('\tfiltered transfrags bed file: %s' % (results.output_dir))
    retcode = sort_bed(tmp_results.transfrags_filtered_bed_file,
                       results.transfrags_filtered_bed_file,
                       num_processes=1,
                       tmp_dir=tmp_dir)
    os.remove(tmp_results.transfrags_filtered_bed_file)
    if retcode != 0:
        raise TacoError('Error running linux sort')
    # remove temporary directories
    logging.debug('\t%s cleaning up' % (results.output_dir))
    os.rmdir(tmp_dir)


def aggregate_parallel(samples, args, results):
    '''
    Process and aggregate GTF input files

    samples: list of Sample objects
    args: from Argparse module. command-line arguments to configure the
          assembly process
    results: Results object containing input and output filenames
    '''
    logging.info('Aggregating in parallel using %d processes' %
                 (args.num_processes))

    if args.filter_splice_juncs and args.ref_genome_fasta_file:
        # test opening FastaFile
        logging.info('Indexing reference genome fasta file (if necessary)')
        fasta_fh = FastaFile(args.ref_genome_fasta_file)
        fasta_fh.close()

    # create queue
    input_queue = JoinableQueue(maxsize=args.num_processes * 2)
    # start worker processes
    procs = []
    worker_results = []
    for i in xrange(args.num_processes):
        worker_id = 'aggregate_worker%03d' % i
        worker_dir = os.path.join(results.tmp_dir, worker_id)
        if not os.path.exists(worker_dir):
            os.makedirs(worker_dir)
        worker_results.append(Results(worker_dir))
        p = Process(target=aggregate_worker,
                    args=(input_queue, args, worker_dir))
        p.start()
        procs.append(p)

    # reference gtf
    if args.ref_gtf_file is not None:
        logging.debug('Reference: %s' % args.ref_gtf_file)
        input_queue.put(Sample(args.ref_gtf_file, Sample.REF_ID))
    # parse samples
    for sample in samples:
        input_queue.put(sample)
    for p in procs:
        input_queue.put(None)
    # close input queue
    input_queue.join()
    input_queue.close()
    # join worker processes
    for p in procs:
        p.join()

    # merge output files
    logging.info('Merging aggregated files')
    logging.debug('\tmerging bed files')
    retcode = merge_bed(input_files=[r.transfrags_bed_file for r in worker_results],
                        output_file=results.transfrags_bed_file,
                        num_processes=args.num_processes,
                        tmp_dir=results.tmp_dir)
    if retcode != 0:
        raise TacoError('Error running linux merge')

    logging.debug('\tmerging filtered bed files')
    retcode = merge_bed(input_files=[r.transfrags_filtered_bed_file for r in worker_results],
                        output_file=results.transfrags_filtered_bed_file,
                        num_processes=args.num_processes,
                        tmp_dir=results.tmp_dir)
    if retcode != 0:
        raise TacoError('Error running linux merge')

    # cleanup worker data
    logging.info('Removing temporary files')
    def shutil_error_callback(func, path, excinfo):
        logging.error('Error removing tmp files path=%s message=%s' %
                      (path, excinfo))
    for r in worker_results:
        shutil.rmtree(r.output_dir, onerror=shutil_error_callback)
    logging.info('Aggregate done')
    return 0
