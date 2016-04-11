'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
import collections
import logging

from base import Exon, Strand, Sample
from gtf import GTF, GTFError

__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


class Transfrag(object):
    __slots__ = ('chrom', 'strand', 'exons', '_id', 'sample_id', 'expr',
                 'is_ref')

    def __init__(self, chrom=None, strand=None, _id=None, sample_id=None,
                 expr=0.0, is_ref=False, exons=None):
        self.chrom = chrom
        self.strand = Strand.NA if strand is None else strand
        self._id = _id
        self.sample_id = sample_id
        self.expr = expr
        self.is_ref = is_ref
        self.exons = [] if exons is None else exons

    @property
    def length(self):
        return sum((e.end - e.start) for e in self.exons)

    @property
    def start(self):
        return self.exons[0].start

    @property
    def end(self):
        return self.exons[-1].end

    def iterintrons(self):
        e1 = self.exons[0]
        for j in xrange(1, len(self.exons)):
            e2 = self.exons[j]
            yield e1.end, e2.start
            e1 = e2

    def itersplices(self):
        e1 = self.exons[0]
        for j in xrange(1, len(self.exons)):
            e2 = self.exons[j]
            yield e1.end
            yield e2.start
            e1 = e2

    @property
    def txstart(self):
        '''start of transcription (adjust for strand)'''
        if self.strand == Strand.NEG:
            return self.exons[-1][1]
        else:
            return self.exons[0][0]

    @property
    def txstop(self):
        '''end of transcription (adjust for strand)'''
        if self.strand == Strand.NEG:
            return self.exons[0][0]
        else:
            return self.exons[-1][1]

    def to_bed(self):
        tx_start = self.exons[0].start
        tx_end = self.exons[-1].end
        block_sizes = []
        block_starts = []
        for e in self.exons:
            block_starts.append(e.start - tx_start)
            block_sizes.append(e.end - e.start)
        # make bed fields
        fields = [self.chrom,
                  str(tx_start),
                  str(tx_end),
                  self._id,
                  str(self.expr),
                  Strand.to_gtf(self.strand),
                  '0',
                  '0',
                  '0',
                  str(len(self.exons)),
                  ','.join(map(str, block_sizes)),
                  ','.join(map(str, block_starts))]
        return fields

    @staticmethod
    def from_bed(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        tx_start = int(fields[1])
        _id = fields[3]
        is_ref = (_id.split('.')[0] == Sample.REF_ID)
        expr = float(fields[4])
        strand = Strand.from_bed(fields[5])
        num_exons = int(fields[9])
        block_sizes = map(int, fields[10].split(','))
        block_starts = map(int, fields[11].split(','))
        exons = []
        for i in xrange(num_exons):
            start = tx_start + block_starts[i]
            end = start + block_sizes[i]
            exons.append(Exon(start, end))
        return Transfrag(chrom=chrom,
                         strand=strand,
                         _id=_id,
                         expr=expr,
                         is_ref=is_ref,
                         exons=exons)

    def to_gtf(self):
        strand_str = Strand.to_gtf(self.strand)
        f = GTF.Feature()
        f.seqid = self.chrom
        f.source = 'taco'
        f.feature = 'transcript'
        f.start = self.start
        f.end = self.end
        f.score = 0.0
        f.strand = strand_str
        f.phase = '.'
        f.attrs = {GTF.Attr.TRANSCRIPT_ID: self._id,
                   GTF.Attr.SAMPLE_ID: self.sample_id,
                   GTF.Attr.EXPR: str(self.expr),
                   GTF.Attr.REF: str(int(self.is_ref))}
        yield f
        for e in self.exons:
            f = GTF.Feature()
            f.seqid = self.chrom
            f.source = 'taco'
            f.feature = 'exon'
            f.start = e.start
            f.end = e.end
            f.score = 0.0
            f.strand = strand_str
            f.phase = '.'
            f.attrs = {GTF.Attr.TRANSCRIPT_ID: self._id}
            yield f

    @staticmethod
    def from_gtf(f):
        '''GTF.Feature object to Transfrag'''
        return Transfrag(chrom=f.seqid,
                         strand=Strand.from_gtf(f.strand),
                         _id=f.attrs[GTF.Attr.TRANSCRIPT_ID],
                         sample_id=f.attrs.get(GTF.Attr.SAMPLE_ID, None),
                         expr=float(f.attrs.get(GTF.Attr.EXPR, 0.0)),
                         is_ref=bool(int(f.attrs.get(GTF.Attr.REF, '0'))),
                         exons=None)

    @staticmethod
    def parse_gtf(gtf_lines):
        '''
        returns list of Transfrag objects
        '''
        t_dict = collections.OrderedDict()
        for gtf_line in gtf_lines:
            f = GTF.Feature.from_str(gtf_line)
            if f.feature == 'transcript':
                t_id = f.attrs[GTF.Attr.TRANSCRIPT_ID]
                if t_id in t_dict:
                    raise GTFError("Transcript '%s' duplicate detected" % t_id)
                t = Transfrag.from_gtf(f)
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
        return t_dict
