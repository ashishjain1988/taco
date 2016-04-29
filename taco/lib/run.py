'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
'''
import os
import argparse
import logging
import pickle
import json
import shutil

from base import Sample, Results
from locus import Locus
from aggregate import aggregate_parallel
from assemble import assemble_parallel
from clocusindex import bed_index_loci

__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "GPL"
__version__ = "0.5.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class Args:
    VERBOSE = False
    NUM_PROCESSES = 1
    GTF_EXPR_ATTR = 'FPKM'
    FILTER_SPLICE_JUNCS = False
    FILTER_MIN_LENGTH = 200
    FILTER_MIN_EXPR = 0.5
    ISOFORM_FRAC = 0.05
    MAX_ISOFORMS = 0
    ASSEMBLE_UNSTRANDED = False
    CHANGE_POINT = True
    CHANGE_POINT_PVALUE = 0.01
    CHANGE_POINT_FOLD_CHANGE = 0.85
    CHANGE_POINT_TRIM = True
    PATH_GRAPH_KMAX = 0
    PATH_FRAC = 0.0
    MAX_PATHS = 0
    GUIDED_ASSEMBLY = False
    GUIDED_STRAND = False
    GUIDED_ENDS = False
    RESUME = False
    ASSEMBLE = None
    OUTPUT_DIR = 'output'
    PROG = 'taco'
    DESCRIPTION = 'TACO: Multi-sample transcriptome assembly from RNA-Seq'

    @staticmethod
    def create():
        parser = argparse.ArgumentParser(description=Args.DESCRIPTION)
        advgrp = parser.add_argument_group('Advanced Options',
                                           '(recommend leaving at their '
                                           'default settings for most '
                                           'purposes)')
        # guidedgrp = parser.add_argument_group('Guided-Assembly Options',
        #                                       '(DO NOT USE: under '
        #                                       'development)')
        parser.add_argument("-o", "--output-dir", dest="output_dir",
                            metavar='DIR',
                            default=Args.OUTPUT_DIR,
                            help='directory where output files will be '
                            'stored (if already exists then --resume must '
                            'be specified) [default=%(default)s]')
        parser.add_argument('-p', '--num-processes', type=int,
                            metavar='N',
                            dest='num_processes',
                            default=Args.NUM_PROCESSES,
                            help='Run TACO in parallel with N '
                            'processes [default=%(default)s]')
        parser.add_argument('-v', '--verbose', dest='verbose',
                            action="store_true",
                            default=Args.VERBOSE,
                            help='Enabled detailed logging '
                            '(for debugging)')
        parser.add_argument('--resume', dest='resume',
                            action='store_true',
                            default=Args.RESUME,
                            help='Resumes an existing run that may have '
                            'ended prematurely. Specify the location of the '
                            'run using the -o/--output-dir option.')
        parser.add_argument('--assemble',
                            dest='assemble',
                            metavar='BED',
                            default=Args.ASSEMBLE,
                            help='Assemble transfrags produced by a previous '
                            'TACO run, bypassing the GTF aggregation '
                            'step. Accepts a taco-formatted BED file.')
        parser.add_argument('sample_file', nargs='?')
        parser.add_argument('--gtf-expr-attr',
                            dest='gtf_expr_attr',
                            default=Args.GTF_EXPR_ATTR,
                            metavar='ATTR',
                            help='GTF attribute field containing '
                            'expression [default=%(default)s]')
        parser.add_argument('--filter-min-length',
                            dest='filter_min_length',
                            type=int, metavar='N',
                            default=Args.FILTER_MIN_LENGTH,
                            help='Filter input transfrags with length < N '
                            'prior to assembly [default=%(default)s]')
        parser.add_argument('--filter-min-expr',
                            dest='filter_min_expr',
                            type=float, metavar='X',
                            default=Args.FILTER_MIN_EXPR,
                            help='Filter input transfrags with transcripts '
                            'per milliion (TPM) < X prior to assembly '
                            '[default=%(default)s]')
        parser.add_argument('--filter-splice-juncs',
                            dest='filter_splice_juncs',
                            action='store_true',
                            default=Args.FILTER_SPLICE_JUNCS,
                            help='Filter input transfrags that possess '
                            'non-canonical splice motifs prior to assembly. '
                            'Splice motifs are GTAG and GCAG are allowed '
                            '[default=%(default)s]. Requires genome sequence '
                            'to be specified using --ref-genome-fasta.')
        parser.add_argument('--ref-genome-fasta',
                            dest='ref_genome_fasta_file',
                            default=None,
                            help='Reference genome sequence in FASTA format '
                            'needed to assess splice junction motif sequences. '
                            'Use in conjunction with --filter-splice-juncs.')
        parser.add_argument('--isoform-frac',
                            dest='isoform_frac', type=float, metavar='FRAC',
                            default=Args.ISOFORM_FRAC,
                            help='Report transcript isoforms with '
                            'expression fraction >= FRAC (0.0-1.0) relative '
                            'to the major isoform within each gene '
                            '[default=%(default)s]')
        parser.add_argument('--max-isoforms', type=int, metavar='N',
                            dest='max_isoforms',
                            default=Args.MAX_ISOFORMS,
                            help='Maximum isoforms to report for each '
                            'gene [default=%(default)s]')
        parser.add_argument('--assemble-unstranded',
                            dest='assemble_unstranded',
                            action='store_true',
                            default=Args.ASSEMBLE_UNSTRANDED,
                            help='Enable assembly of unstranded transfrags '
                            '[default=%(default)s]')
        parser.add_argument('--no-assemble-unstranded',
                            dest='assemble_unstranded',
                            action='store_false',
                            help='Disable assembly of unstranded transfrags')
        parser.add_argument('--change-point', dest='change_point',
                            action='store_true',
                            default=Args.CHANGE_POINT,
                            help='Enable change point detection [default='
                            '%(default)s]')
        parser.add_argument('--no-change-point', dest='change_point',
                            action='store_false',
                            help='Disable change point detection')
        advgrp.add_argument('--change-point-pvalue', type=float,
                            dest='change_point_pvalue',
                            default=Args.CHANGE_POINT_PVALUE,
                            metavar='<float 0.0-1.0>',
                            help='Mann-Whitney-U p-value threshold for '
                            'calling change points [default=%(default)s]')
        advgrp.add_argument('--change-point-fold-change', type=float,
                            dest='change_point_fold_change',
                            default=Args.CHANGE_POINT_FOLD_CHANGE,
                            metavar='<float 0.0-1.0>',
                            help='Fold change threshold between the means of '
                            'two putative segments on either side of a change '
                            'point. A value of 0.0 is the most strict '
                            'setting, effectively calling no change points. '
                            'Conversely, setting the value to 1.0 calls all'
                            'change points [default=%(default)s]')
        advgrp.add_argument('--change-point-trim', dest='change_point_trim',
                            action='store_true',
                            default=Args.CHANGE_POINT_TRIM,
                            help='Trim transfrags around change points '
                            '[default=%(default)s]')
        advgrp.add_argument('--no-change-point-trim', dest='change_point_trim',
                            action='store_false',
                            help='Disable trimming around change points')
        advgrp.add_argument('--path-kmax', dest='path_graph_kmax', type=int,
                            metavar='kmax',
                            default=Args.PATH_GRAPH_KMAX,
                            help='Limit optimization for choosing parameter k '
                            'for path graph (DeBruijn graph) to k <= kmax '
                            '[default=%(default)s]')
        advgrp.add_argument('--path-frac', type=float, metavar='X',
                            dest='path_frac',
                            default=Args.PATH_FRAC,
                            help='dynamic programming algorithm will stop '
                            'finding suboptimal paths when path expression '
                            'drops below a fraction X (0.0-1.0) of the total '
                            'locus expression [default=%(default)s]')
        advgrp.add_argument('--max-paths', type=int, metavar='N',
                            dest='max_paths',
                            default=Args.MAX_PATHS,
                            help='dynamic programming algorithm will stop '
                            'after finding N paths [default=%(default)s]')
        parser.add_argument('--ref-gtf', dest='ref_gtf_file',
                            metavar='<gtf_file>',
                            default=None,
                            help=argparse.SUPPRESS)
                            # help='Option reference GTF file of "true" '
                            # 'validated transcripts to facilitate guided '
                            # 'assembly and/or noise filtering')
        parser.add_argument('--guided-strand', dest='guided_strand',
                            action='store_true',
                            default=Args.GUIDED_STRAND,
                            help=argparse.SUPPRESS)
                            # help='Enable use of reference strand information '
                            # 'to help resolve unstranded transfrags (requires '
                            # 'reference GTF to be specified using --ref-gtf)')
        parser.add_argument('--guided-ends', dest='guided_ends',
                            action='store_true',
                            default=Args.GUIDED_ENDS,
                            help=argparse.SUPPRESS)
                            # help='Enable use of reference transcript start '
                            # 'and end sites during assembly (requires '
                            # 'reference GTF to be specified using --ref-gtf)')
        parser.add_argument('--guided-assembly', dest='guided_assembly',
                            action='store_true',
                            default=Args.GUIDED_ASSEMBLY,
                            help=argparse.SUPPRESS)
                            # help='Enable guided assembly (requires a '
                            # 'reference GTF to be specified using '
                            # '--ref-gtf)')
        return parser

    @staticmethod
    def load(filename):
        return pickle.load(open(filename))

    @staticmethod
    def dump(args, filename):
        pickle.dump(args, open(filename, 'wb'))

    @staticmethod
    def log(args, func=logging.info):
        logging.info('%s version %s' % (Args.PROG, __version__))
        spacer = '-' * 78
        fmt = '{:<35}{:<35}'
        func(spacer)
        func(fmt.format('verbose logging:', str(args.verbose)))
        func(fmt.format('num processes:', str(args.num_processes)))
        func(fmt.format('output directory:', str(args.output_dir)))
        func(fmt.format('filter min length:', args.filter_min_length))
        func(fmt.format('filter min expression:', args.filter_min_expr))
        func(fmt.format('filter splice juncs:', args.filter_splice_juncs))
        func(fmt.format('reference genome FASTA file:',
                        args.ref_genome_fasta_file))
        func(fmt.format('reference GTF file:', str(args.ref_gtf_file)))
        func(fmt.format('guided assembly mode:', str(args.guided_assembly)))
        func(fmt.format('guided strand mode:', str(args.guided_strand)))
        func(fmt.format('guided ends mode:', str(args.guided_ends)))
        func(fmt.format('GTF expression attribute:', args.gtf_expr_attr))
        func(fmt.format('isoform fraction:', args.isoform_frac))
        func(fmt.format('max_isoforms:', args.max_isoforms))
        func(fmt.format('assemble_unstranded:', args.assemble_unstranded))
        func(fmt.format('change point:', str(args.change_point)))
        func(fmt.format('change point pvalue:', str(args.change_point_pvalue)))
        func(fmt.format('change point fold change:',
                        str(args.change_point_fold_change)))
        func(fmt.format('change point trim:', str(args.change_point_trim)))
        func(fmt.format('path frac:', args.path_frac))
        func(fmt.format('max paths:', args.max_paths))

    @staticmethod
    def parse():
        parser = Args.create()
        args = parser.parse_args()
        if args.resume:
            # check if we are trying to resume a previous run
            if not os.path.exists(args.output_dir):
                parser.error("Output directory '%s' does not exist" %
                             args.output_dir)
            # check that basic files exist in directory
            results = Results(args.output_dir)
            can_resume = os.path.exists(results.args_file)
            can_resume = can_resume and os.path.exists(results.status_file)
            can_resume = can_resume and os.path.exists(results.sample_file)
            if not can_resume:
                parser.error("Output directory '%s' not valid" %
                             args.output_dir)
        else:
            if os.path.exists(args.output_dir):
                parser.error("Output directory '%s' already exists" %
                             args.output_dir)
            # check if we are trying bypass aggregate and assemble from a
            # preexisting BED file
            if args.assemble is not None:
                if not os.path.exists(args.assemble):
                    parser.error('(--assemble) BED file %s not found')
                args.assemble = os.path.abspath(args.assemble)

            if args.sample_file is None:
                parser.error('Sample file not specified. Run with --help for more info.')
            if not os.path.exists(args.sample_file):
                parser.error('Sample file %s not found' % (args.sample_file))
            args.sample_file = os.path.abspath(args.sample_file)

            if args.filter_min_length < 0:
                parser.error("filter_min_length < 0")
            if args.filter_min_expr < 0:
                parser.error("filter_min_expr < 0")
            if (args.isoform_frac < 0) or (args.isoform_frac > 1):
                parser.error("isoform_frac out of range (0.0-1.0)")
            if (args.max_isoforms < 0):
                parser.error("max_isoforms < 0")

            if args.path_graph_kmax < 0:
                parser.error("kmax must be >= 0")
            if args.max_paths < 0:
                parser.error("max_paths must be >= 0")
            if not (0 <= args.path_frac <= 1):
                parser.error("path_frac not in range (0.0-1.0)")

            if args.change_point:
                if not (0.0 <= args.change_point_pvalue <= 1.0):
                    parser.error('change point pvalue invalid')
                if not (0.0 <= args.change_point_fold_change <= 1.0):
                    parser.error('change point fold change invalid')

            if args.filter_splice_juncs:
                if not args.ref_genome_fasta_file:
                    parser.error('Filtering of splice junctions enabled '
                                 '(--filter-splice-juncs) but '
                                 'reference genome FASTA file not specified')
                if not os.path.exists(args.ref_genome_fasta_file):
                    parser.error('Reference genome FASTA file not found')
                args.ref_genome_fasta_file = os.path.abspath(args.ref_genome_fasta_file)

            if args.ref_gtf_file is not None:
                if not os.path.exists(args.ref_gtf_file):
                    parser.error("reference GTF file %s not found" %
                                 (args.ref_gtf_file))
                args.ref_gtf_file = os.path.abspath(args.ref_gtf_file)
            elif (args.guided_assembly or args.guided_ends or
                  args.guided_strand):
                parser.error('Guided assembly modes require a '
                             'reference GTF (--ref-gtf)')
        return args


class Status(object):
    FIELDS = ('create', 'aggregate', 'index_loci', 'assemble')

    def __init__(self):
        for f in Status.FIELDS:
            setattr(self, f, False)

    def write(self, filename):
        d = dict((f, getattr(self, f)) for f in Status.FIELDS)
        json.dump(d, open(filename, 'w'))

    @staticmethod
    def load(filename):
        d = json.load(open(filename))
        self = Status()
        for f in Status.FIELDS:
            setattr(self, f, d[f])
        return self


class Run(object):

    @staticmethod
    def create():
        self = Run()
        # parse command line args
        args = Args.parse()
        # setup logging
        if args.verbose:
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(level=level,
                            format="%(asctime)s pid=%(process)d "
                                   "%(levelname)s - %(message)s")

        if args.resume:
            self.results = Results(args.output_dir)
            self.status = Status.load(self.results.status_file)
            self.samples = Sample.parse_tsv(self.results.sample_file,
                                            header=True)
            self.args = Args.load(self.results.args_file)
        else:
            self.args = args
            self.results = Results(args.output_dir)
            self.status = Status()
            self.samples = Sample.parse_tsv(self.args.sample_file)

        # create output directories
        results = self.results
        if not os.path.exists(results.output_dir):
            logging.debug("Creating output directory '%s'" %
                          (results.output_dir))
            os.makedirs(results.output_dir)
        if not os.path.exists(results.tmp_dir):
            logging.debug("Creating tmp directory '%s'" % (results.tmp_dir))
            os.makedirs(results.tmp_dir)

        if not args.resume:
            # write command line args
            Args.dump(self.args, self.results.args_file)
            # write samples
            Sample.write_tsv(self.samples, self.results.sample_file)
            # bypass aggregate step
            if args.assemble is not None:
                # create a symbolic link to BED file
                os.symlink(args.assemble, self.results.transfrags_bed_file)
                self.status.aggregate = True
            # update status and write to file
            self.status.create = True
            self.status.write(self.results.status_file)
        return self

    def log_args(self):
        Args.log(self.args)

    def aggregate(self):
        '''
        Aggregate/merge individual sample GTF files
        '''
        aggregate_parallel(self.samples, self.args, self.results)
        # update status and write to file
        self.status.aggregate = True
        self.status.write(self.results.status_file)

    def index_loci(self):
        r = self.results
        retcode = bed_index_loci(r.transfrags_bed_file, r.locus_index_file)
        if retcode != 0:
            raise TacoError('Error indexing bed file')
        # update status and write to file
        self.status.index_loci = True
        self.status.write(self.results.status_file)

    def assemble(self):
        assemble_parallel(self.args, self.results, len(self.samples))
        # update status and write to file
        self.status.assemble = True
        self.status.write(self.results.status_file)
