'''
TACO: Transcriptome meta-assembly from RNA-Seq
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


def parse_gtf(gtf_iter, sample_id, gtf_expr_attr, is_ref):
    '''
    returns list of Transfrag objects
    '''
    t_dict = collections.OrderedDict()
    total_expr = 0.0
    cur_t_id = 1
    for gtf_line in gtf_iter:
        f = GTF.Feature.from_str(gtf_line)
        if f.feature == 'transcript':
            t_id = f.attrs[GTF.Attr.TRANSCRIPT_ID]
            if t_id in t_dict:
                raise GTFError("Transcript '%s' duplicate detected" % t_id)
            # rename transcript id
            new_t_id = "%s.%d" % (sample_id, cur_t_id)
            cur_t_id += 1
            # parse expression
            if is_ref:
                expr = 0.0
            else:
                if gtf_expr_attr not in f.attrs:
                    raise GTFError("GTF expression attribute '%s' not found" %
                                   (gtf_expr_attr))
                expr = float(f.attrs[gtf_expr_attr])
                total_expr += expr
            # create transfrag
            t = Transfrag(chrom=f.seqid,
                          strand=Strand.from_gtf(f.strand),
                          _id=new_t_id,
                          expr=float(expr),
                          is_ref=is_ref,
                          exons=None)
            t_dict[t_id] = t
        elif f.feature == 'exon':
            t_id = f.attrs[GTF.Attr.TRANSCRIPT_ID]
            if t_id not in t_dict:
                logging.error('Feature: "%s"' % str(f))
                raise GTFError("Transcript '%s' exon feature appeared in "
                               "gtf file prior to transcript feature" %
                               t_id)
            t = t_dict[t_id]
            t.exons.append(Exon(f.start, f.end))
    return t_dict.values(), total_expr


def aggregate_sample(sample, gtf_expr_attr, is_ref, min_length, min_expr,
                     filter_splice_juncs, genome_fasta_fh,
                     bed_fh, filtered_bed_fh, stats_fh):
    logging.debug('Aggregate sample %s: %s' % (sample._id, sample.gtf_file))
    # read all transcripts
    with open(sample.gtf_file) as fh:
        transcripts, total_expr = parse_gtf(fh, sample._id, gtf_expr_attr, is_ref)

    # track filtering stats
    nlength = 0
    nexpr = 0
    nsplice = 0
    for t in transcripts:
        # normalize expression
        if total_expr > 0:
            t.expr = 1.0e6 * t.expr / total_expr

        # check filter conditions
        keep = True
        if t.length < min_length:
            keep = False
            nlength += 1
        if t.expr < min_expr:
            keep = False
            nexpr += 1
        if filter_splice_juncs:
            # remove transfrags with non-canonical splice motifs
            for start, end in t.iterintrons():
                s = genome_fasta_fh.fetch(t.chrom, start, start + 2)
                s += genome_fasta_fh.fetch(t.chrom, end - 2, end)
                if t.strand == Strand.NEG:
                    s = dna_reverse_complement(s)
                if s not in SPLICE_MOTIFS_ALLOWED:
                    keep = False
                    nsplice += 1

        # write transcript to bed
        line = '\t'.join(t.to_bed())
        if keep:
            print >>bed_fh, line
        else:
            print >>filtered_bed_fh, line
    fields = [sample._id, len(transcripts), nlength, nexpr, nsplice]
    print >>stats_fh, '\t'.join(map(str, fields))


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
    # setup output files
    bed_fh = open(tmp_results.transfrags_bed_file, 'w')
    filtered_bed_fh = open(tmp_results.transfrags_filtered_bed_file, 'w')
    stats_fh = open(results.sample_stats_file, 'w')
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
                         bed_fh=bed_fh,
                         filtered_bed_fh=filtered_bed_fh,
                         stats_fh=stats_fh)
        input_queue.task_done()
    input_queue.task_done()
    # cleanup and close files
    bed_fh.close()
    filtered_bed_fh.close()
    stats_fh.close()
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

    logging.debug('\tmerging sample stats')
    def sort_key_field0(line):
        fields = line.split('\t', 1)
        return fields[0]
    stats_header = ['sample_id', 'num_transfrags', 'filtered_length',
                    'filtered_expr', 'filtered_splice']
    stats_header = '\t'.join(stats_header)
    merge_files(input_files=[r.sample_stats_file for r in worker_results],
                output_file=results.sample_stats_file,
                key=sort_key_field0,
                header=stats_header)
    # cleanup worker data
    logging.info('Removing temporary files')
    def shutil_error_callback(func, path, excinfo):
        logging.error('Error removing tmp files path=%s message=%s' %
                      (path, excinfo))
    for r in worker_results:
        shutil.rmtree(r.output_dir, onerror=shutil_error_callback)
    logging.info('Aggregate done')
    return 0
