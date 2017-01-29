'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
'''
import collections
import logging
import os

__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "GPL"
__version__ = "0.5.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


# assemblyline GTF attributes
class GTFAttr:
    GENE_ID = 'gene_id'
    TRANSCRIPT_ID = 'transcript_id'
    SAMPLE_ID = 'sid'
    LIBRARY_ID = 'lid'
    REF = 'ref'
    TEST = 'tst'
    PCTRANK = 'pct'
    CATEGORY = 'cat'
    ANN_REF_ID = 'aid'
    ANN_COV_RATIO = 'acr'
    ANN_INTRON_RATIO = 'air'
    SCORE = 'score'
    MEAN_SCORE = 'avgscore'
    MEAN_PCTRANK = 'avgpct'
    MEAN_RECURRENCE = 'avgrecur'
    RESOLVED_STRAND = 'resolvedstrand'
    LOG10LR = 'log10lr'

class Category(object):
    # constant transcript category values
    SAME_STRAND = 0
    OPP_STRAND = 1
    INTRONIC_SAME_STRAND = 2
    INTRONIC_OPP_STRAND = 3
    INTRONIC_AMBIGUOUS = 4
    INTERLEAVING = 5
    INTERGENIC = 6
    NUM_CATEGORIES = 7
    CATEGORIES = range(0, NUM_CATEGORIES)

    # used when comparing assemblies
    INTERLEAVING_SAME_STRAND = 8
    INTERLEAVING_OPP_STRAND = 9
    ENCOMPASSING_SAME_STRAND = 10
    ENCOMPASSING_OPP_STRAND = 11
    # used when a gene overlaps multiple independent genes in the SAME_STRAND category
    READ_THROUGH = 12

    # groups of categories
    INTRONIC_LIKE = set([INTRONIC_SAME_STRAND, INTRONIC_AMBIGUOUS])
    INTERGENIC_LIKE = set([OPP_STRAND, INTRONIC_OPP_STRAND, INTERLEAVING,
                           INTERGENIC])

    INT_TO_STR_DICT = {SAME_STRAND: 'same_strand',
                       OPP_STRAND: 'opp_strand',
                       INTRONIC_SAME_STRAND: 'intronic_same_strand',
                       INTRONIC_OPP_STRAND: 'intronic_opp_strand',
                       INTRONIC_AMBIGUOUS: 'intronic_ambiguous',
                       INTERLEAVING: 'interleaving',
                       INTERGENIC: 'intergenic',
                       INTERLEAVING_SAME_STRAND: 'interleaving_same_strand',
                       INTERLEAVING_OPP_STRAND: 'interleaving_opp_strand',
                       ENCOMPASSING_SAME_STRAND: 'encompassing_same_strand',
                       ENCOMPASSING_OPP_STRAND: 'encompassing_opp_strand',
                       READ_THROUGH: 'read_through'}

    STR_TO_INT_DICT = dict((v,k) for k,v in INT_TO_STR_DICT.items())

    @staticmethod
    def to_str(catint):
        return Category.INT_TO_STR_DICT[catint]
    @staticmethod
    def to_int(catstr):
        return Category.STR_TO_INT_DICT[catstr]


class TacoError(Exception):
    pass


Exon = collections.namedtuple('Exon', ['start', 'end'])


class Strand:
    POS = 0
    NEG = 1
    NA = 2

    STRANDS = [POS, NEG, NA]
    NAMES = ['pos', 'neg', 'none']
    FROM_GTF = {'+': POS, '-': NEG, '.': NA}
    TO_GTF = {POS: '+', NEG: '-', NA: '.'}
    FROM_BED = {'+': POS, '-': NEG, '.': NA}
    TO_BED = {POS: '+', NEG: '-', NA: '.'}

    @staticmethod
    def to_str(s):
        return Strand.NAMES[s]

    @staticmethod
    def from_gtf(s):
        return Strand.FROM_GTF[s]

    @staticmethod
    def to_gtf(s):
        return Strand.TO_GTF[s]

    @staticmethod
    def from_bed(s):
        return Strand.FROM_BED[s]

    @staticmethod
    def to_bed(s):
        return Strand.TO_BED[s]



class Sample(object):
    REF_ID = 'R'

    def __init__(self, gtf_file, _id):
        self.gtf_file = gtf_file
        self._id = _id

    @staticmethod
    def parse_tsv(filename, header=False, sep='\t'):
        cur_sample_id = 1
        gtf_files = set()
        ids = set()
        samples = []
        with open(filename) as f:
            if header:
                f.next()
            # table rows
            for line in f:
                fields = line.strip().split(sep)
                gtf_file = fields[0]
                if gtf_file in gtf_files:
                    m = "GTF file '%s' is not unique" % gtf_file
                    raise TacoError(m)
                if not Sample.gtf_valid(gtf_file):
                    m = "GTF file '%s' is not valid" % gtf_file
                    raise TacoError(m)
                if len(fields) > 1:
                    _id = fields[1]
                    if _id in ids:
                        m = "sample_id '%s' is not unique" % _id
                        raise TacoError(m)
                else:
                    _id = cur_sample_id
                    cur_sample_id += 1
                samples.append(Sample(gtf_file, _id))
        return samples

    @staticmethod
    def write_tsv(samples, filename, header=True, sep='\t'):
        with open(filename, 'w') as f:
            if header:
                print >>f, sep.join(['gtf', 'sample_id'])
            for s in samples:
                print >>f, sep.join([s.gtf_file, str(s._id)])

    @staticmethod
    def gtf_valid(gtf_file):
        if gtf_file is None:
            logging.error('GTF file %s is None' % (gtf_file))
            return False
        if not os.path.exists(gtf_file):
            logging.error('GTF file %s not found' % (gtf_file))
            return False
        return True


class Results(object):
    TMP_DIR = 'tmp'
    STATUS_FILE = 'status.json'
    ARGS_FILE = 'args.pickle'
    SAMPLE_FILE = 'samples.txt'
    TRANSFRAGS_BED_FILE = 'transfrags.bed'
    TRANSFRAGS_FILTERED_BED_FILE = 'transfrags.filtered.bed'
    SAMPLE_STATS_FILE = 'sample_stats.txt'
    LOCUS_INDEX_FILE = 'loci.txt'
    SPLICE_GRAPH_GTF_FILE = 'splice_graph.gtf'
    CHANGE_POINT_GTF_FILE = 'change_points.gtf'
    BEDGRAPH_FILES = ['expr.pos.bedgraph', 'expr.neg.bedgraph',
                      'expr.none.bedgraph']
    SPLICE_BED_FILE = 'splice_junctions.bed'
    PATH_GRAPH_STATS_FILE = 'path_graph_stats.txt'
    ASSEMBLY_GTF_FILE = 'assembly.gtf'
    ASSEMBLY_BED_FILE = 'assembly.bed'

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.tmp_dir = os.path.join(output_dir, Results.TMP_DIR)
        self.args_file = os.path.join(output_dir, Results.ARGS_FILE)
        self.status_file = os.path.join(output_dir, Results.STATUS_FILE)
        self.sample_file = os.path.join(output_dir, Results.SAMPLE_FILE)
        self.sample_stats_file = \
            os.path.join(output_dir, Results.SAMPLE_STATS_FILE)
        self.transfrags_bed_file = \
            os.path.join(output_dir, Results.TRANSFRAGS_BED_FILE)
        self.transfrags_filtered_bed_file = \
            os.path.join(output_dir, Results.TRANSFRAGS_FILTERED_BED_FILE)
        self.locus_index_file = \
            os.path.join(output_dir, Results.LOCUS_INDEX_FILE)
        self.splice_graph_gtf_file = \
            os.path.join(output_dir, Results.SPLICE_GRAPH_GTF_FILE)
        self.change_point_gtf_file = \
            os.path.join(output_dir, Results.CHANGE_POINT_GTF_FILE)
        self.bedgraph_files = [os.path.join(output_dir, x)
                               for x in Results.BEDGRAPH_FILES]
        self.splice_bed_file = \
            os.path.join(output_dir, Results.SPLICE_BED_FILE)
        self.path_graph_stats_file = \
            os.path.join(output_dir, Results.PATH_GRAPH_STATS_FILE)
        self.assembly_gtf_file = \
            os.path.join(output_dir, Results.ASSEMBLY_GTF_FILE)
        self.assembly_bed_file = \
            os.path.join(output_dir, Results.ASSEMBLY_BED_FILE)
