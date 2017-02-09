
import argparse
import logging
import os
import sys
import operator
import collections
import shutil
import subprocess
import multiprocessing
import taco
import share
import pkgutil
import string
import glob
from taco.lib.bx.intersection import Interval, IntervalTree #1
from taco.lib.bx.cluster import ClusterTree #lx1
from taco.lib.base import Category, GTFAttr #3
from taco.lib.gtf_comp import sort_gtf, GTFFeature #4
from taco.lib.transcript import cmp_strand, parse_gtf, \
    strand_int_to_str, NO_STRAND, POS_STRAND, NEG_STRAND, \
    find_exon_boundaries, split_exons

# for nearest transcripts calculation
MAX_LOCUS_DIST = 100000000

GENCODE_CATEGORY_MAP = {'IG_C_gene': 'protein_coding',
                     'IG_D_gene': 'protein_coding',
                     'IG_J_gene': 'protein_coding',
                     'IG_V_gene': 'protein_coding',
                     'IG_LV_gene': 'protein_coding',
                     'TR_C_gene': 'protein_coding',
                     'TR_J_gene': 'protein_coding',
                     'TR_V_gene': 'protein_coding',
                     'TR_D_gene': 'protein_coding',
                     'TEC': 'protein_coding',
                     'nonsense_mediated_decay': 'protein_coding',
                     'non_stop_decay': 'protein_coding',
                     'retained_intron': 'protein_coding',
                     'protein_coding': 'protein_coding',
                     'ambiguous_orf': 'protein_coding',
                     'processed_transcript': 'protein_coding',
                     'Mt_rRNA': 'ncRNA',
                     'Mt_tRNA': 'ncRNA',
                     'miRNA': 'ncRNA',
                     'misc_RNA': 'ncRNA',
                     'rRNA': 'ncRNA',
                     'snRNA': 'ncRNA',
                     'snoRNA': 'ncRNA',
                     '3prime_overlapping_ncrna': 'ncRNA',
                     'lincRNA': 'lncRNA',
                     'sense_intronic': 'lncRNA',
                     'sense_overlapping': 'lncRNA',
                     'antisense': 'lncRNA',
                     'bidirectional_promoter_lncrna' : 'lncRNA',
                     'ribozyme' : 'lncRNA',
                     'macro_lncRNA': 'lncRNA',
                     'non_coding': 'lncRNA',
                     'bidirectional_promoter_lncRNA': 'lncRNA',
                     'scaRNA': 'ncRNA',
                     'sRNA': 'ncRNA',
                     'IG_pseudogene': 'pseudogene',
                     'IG_C_pseudogene': 'pseudogene',
                     'IG_J_pseudogene': 'pseudogene',
                     'IG_V_pseudogene': 'pseudogene',
                     'TR_V_pseudogene': 'pseudogene',
                     'TR_J_pseudogene': 'pseudogene',
                     'Mt_tRNA_pseudogene': 'pseudogene',
                     'tRNA_pseudogene': 'pseudogene',
                     'snoRNA_pseudogene': 'pseudogene',
                     'snRNA_pseudogene': 'pseudogene',
                     'scRNA_pseudogene': 'pseudogene',
                     'rRNA_pseudogene': 'pseudogene',
                     'misc_RNA_pseudogene': 'pseudogene',
                     'miRNA_pseudogene': 'pseudogene',
                     'pseudogene': 'pseudogene',
                     'processed_pseudogene': 'pseudogene',
                     'polymorphic_pseudogene': 'pseudogene',
                     'retrotransposed': 'pseudogene',
                     'transcribed_processed_pseudogene': 'pseudogene',
                     'transcribed_unprocessed_pseudogene': 'pseudogene',
                     'unitary_pseudogene': 'pseudogene',
                     'transcribed_unitary_pseudogene': 'pseudogene',
                     'translated_unprocessed_pseudogene': 'pseudogene',
                     'unprocessed_pseudogene': 'pseudogene'}


def gencode_category_map(x):
    if 'ncrna' in x.lower():
        out = 'ncRNA'
    elif 'lncrna' in x.lower():
        out = 'lncRNA'
    else:
        out = GENCODE_CATEGORY_MAP.get(x, 'other')
    return out




# attributes to include in TSV that is generated at the end
FULL_GTF_ATTRS = ['gene_id',
                  'tss_id',
                  'annotation',
                  'category',
                  'category_relative',
                  'cateogry_relative_detail',
                  'cpat_coding_prob',
                  'orf_size',
                  'ref_transcript_id',
                  'ref_gene_id',
                  'ref_gene_name',
                  'ref_gene_type',
                  'ref_length',
                  'shared_same_strand_bp',
                  'shared_opp_strand_bp',
                  'shared_introns',
                  'shared_splicing']

def wc(infile):
    p = subprocess.Popen('wc -l %s' % infile, shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()
    return int(out.split()[0])

# class + function to get metadata TSV from gtf
class TranscriptMetadata(object):
    def __init__(self, gtf_attrs=None):
        self.chrom = None
        self.start = 0
        self.end = 0
        self.strand = '.'
        self.num_exons = 0
        self.length = 0
        if gtf_attrs is not None:
            for attr in gtf_attrs:
                setattr(self, attr, '')

def get_gtf_metadata(gtf_file, gtf_attrs):
    if gtf_attrs is None:
        gtf_attrs = []
    if 'transcript_id' in gtf_attrs:
        gtf_attrs.remove('transcript_id')
    # read gtf file
    metadata_dict = {}
    for feature in GTFFeature.parse(open(gtf_file)):
        if feature.feature_type != "exon":
            continue
        t_id = feature.attrs["transcript_id"]
        if t_id not in metadata_dict:
            # instantiate new metadata
            m = TranscriptMetadata()
            m.chrom = feature.seqid
            m.strand = feature.strand
            m.start = feature.start
            m.end = feature.end
            for gtf_attr in gtf_attrs:
                setattr(m, gtf_attr, feature.attrs.get(gtf_attr, 'NA'))
            metadata_dict[t_id] = m
        else:
            m = metadata_dict[t_id]
        # update metadata
        m.start = feature.start if feature.start < m.start else m.start
        m.end = feature.end if feature.end > m.end else m.end
        m.length += (feature.end - feature.start)
        m.num_exons += 1
    return metadata_dict

# function to get bed from GTF for CPAT
def write_bed(chrom, name, strand, score, exons):
    assert all(exons[0].start < x.start for x in exons[1:])
    assert all(exons[-1].end > x.end for x in exons[:-1])
    tx_start = exons[0].start
    tx_end = exons[-1].end
    block_sizes = []
    block_starts = []
    for e in exons:
        block_starts.append(e.start - tx_start)
        block_sizes.append(e.end - e.start)
    # make bed fields
    fields = [chrom,
              str(tx_start),
              str(tx_end),
              str(name),
              str(score),
              strand_int_to_str(strand),
              str(tx_start),
              str(tx_start),
              '0',
              str(len(exons)),
              ','.join(map(str,block_sizes)) + ',',
              ','.join(map(str,block_starts)) + ',']
    return fields

class CompareData(object):
    __slots__ = ('has_ref', 'has_test', 'category')
    def __init__(self):
        self.has_ref = False
        self.has_test = False
        self.category = None

class GlobalStats(object):
    FIELDS = ('introns_both', 'introns_ref_only', 'introns_test_only',
              'patterns_both', 'patterns_ref_only', 'patterns_test_only',
              'cov_both', 'cov_ref_only', 'cov_test_only')

    def __init__(self):
        for field in GlobalStats.FIELDS:
            setattr(self, field, 0)
        self.introns_by_category = collections.defaultdict(lambda: 0)
        self.patterns_by_category = collections.defaultdict(lambda: 0)
        self.cov_by_category = collections.defaultdict(lambda: 0)

    def report(self):
        # print stats report
        introns_total = self.introns_both + self.introns_ref_only + self.introns_test_only
        patterns_total = self.patterns_both + self.patterns_ref_only + self.patterns_test_only
        cov_total = self.cov_both + self.cov_ref_only + self.cov_test_only
        lines = ["introns_total=%d" % (introns_total),
                 "introns_both=%d" % (self.introns_both),
                 "introns_ref_only=%d" % (self.introns_ref_only),
                 "introns_test_only=%d" % (self.introns_test_only),
                 "introns_precision=%f" % (self.introns_both / float(max(1,self.introns_both + self.introns_test_only))),
                 "introns_recall=%f" % (self.introns_both / float(max(1,self.introns_both + self.introns_ref_only))),
                 "patterns_total=%d" % (patterns_total),
                 "patterns_both=%d" % (self.patterns_both),
                 "patterns_ref_only=%d" % (self.patterns_ref_only),
                 "patterns_test_only=%d" % (self.patterns_test_only),
                 "patterns_precision=%f" % (self.patterns_both / float(max(1,self.patterns_both + self.patterns_test_only))),
                 "patterns_recall=%f" % (self.patterns_both / float(max(1,self.patterns_both + self.patterns_ref_only))),
                 "cov_total=%d" % (cov_total),
                 "cov_both=%d" % (self.cov_both),
                 "cov_ref_only=%d" % (self.cov_ref_only),
                 "cov_test_only=%d" % (self.cov_test_only),
                 "cov_precision=%f" % (self.cov_both / float(max(1,self.cov_both + self.cov_test_only))),
                 "cov_recall=%f" % (self.cov_both / float(max(1,self.cov_both + self.cov_ref_only)))]
        for k in sorted(self.introns_by_category):
            lines.append("introns %s=%d" % (k,self.introns_by_category[k]))
        for k in sorted(self.patterns_by_category):
            lines.append("patterns %s=%d" % (k,self.patterns_by_category[k]))
        for k in sorted(self.cov_by_category):
            lines.append("cov %s=%d" % (k,self.cov_by_category[k]))
        return '\n'.join(lines)

    @staticmethod
    def from_file(filename):
        self = GlobalStats()
        with open(filename) as f:
            f.next() # introns total
            self.introns_both = int(f.next().split('=')[1])
            self.introns_ref_only = int(f.next().split('=')[1])
            self.introns_test_only = int(f.next().split('=')[1])
            f.next() # prec
            f.next() # recall
            f.next() # patterns total
            self.patterns_both = int(f.next().split('=')[1])
            self.patterns_ref_only = int(f.next().split('=')[1])
            self.patterns_test_only = int(f.next().split('=')[1])
            f.next() # prec
            f.next() # recall
            f.next() # cov total
            self.cov_both = int(f.next().split('=')[1])
            self.cov_ref_only = int(f.next().split('=')[1])
            self.cov_test_only = int(f.next().split('=')[1])
        return self

    def compute(self, transcripts):
        intron_dict = collections.defaultdict(lambda: CompareData())
        node_dict = collections.defaultdict(lambda: CompareData())
        splicing_pattern_dict = collections.defaultdict(lambda: CompareData())
        # find the intron domains of the transcripts
        boundaries = find_exon_boundaries(transcripts)
        unstranded_transcripts = []
        for t in transcripts:
            if t.strand == NO_STRAND:
                unstranded_transcripts.append(t)
                continue
            # separate ref and nonref transcripts
            is_ref = bool(int(t.attrs[GTFAttr.REF]))
            # split exons that cross boundaries and get the
            # nodes in the transcript path
            for n in split_exons(t, boundaries):
                n = (t.strand, n[0], n[1])
                if is_ref:
                    node_dict[n].has_ref = True
                else:
                    d = node_dict[n]
                    d.has_test = True
                    d.category = t.attrs['category']
            splicing_pattern = []
            for start,end in t.iterintrons():
                n = (t.strand, start, end)
                if is_ref:
                    intron_dict[n].has_ref = True
                else:
                    d = intron_dict[n]
                    d.has_test = True
                    d.category = t.attrs['category']
                splicing_pattern.append(n)
            splicing_pattern = tuple(splicing_pattern)
            if len(splicing_pattern) > 0:
                if is_ref:
                    splicing_pattern_dict[splicing_pattern].has_ref = True
                else:
                    d = splicing_pattern_dict[splicing_pattern]
                    d.has_test = True
                    d.category = t.attrs['category']
        # handle unstranded transcripts
        for t in unstranded_transcripts:
            # separate ref and nonref transcripts
            is_ref = bool(int(t.attrs[GTFAttr.REF]))
            for n in split_exons(t, boundaries):
                found_node = False
                for strand in (POS_STRAND, NEG_STRAND):
                    sn = (strand, n[0], n[1])
                    if sn in node_dict:
                        if is_ref:
                            node_dict[sn].has_ref = True
                        else:
                            d = node_dict[sn]
                            d.has_test = True
                            d.category = t.attrs['category']
                        found_node = True
                if not found_node:
                    sn = (NO_STRAND, n[0], n[1])
                    if is_ref:
                        node_dict[sn].has_ref = True
                    else:
                        d = node_dict[sn]
                        d.has_test = True
                        d.category = t.attrs['category']
            introns = list(t.iterintrons())
            assert len(introns) == 0
        # compile statistics
        for d in intron_dict.itervalues():
            if d.has_ref and d.has_test:
                self.introns_both += 1
                self.introns_by_category[d.category] += 1
            elif d.has_ref:
                self.introns_ref_only += 1
            elif d.has_test:
                self.introns_test_only += 1
                self.introns_by_category[d.category] += 1
        for d in splicing_pattern_dict.itervalues():
            if d.has_ref and d.has_test:
                self.patterns_both += 1
                self.patterns_by_category[d.category] += 1
            elif d.has_ref:
                self.patterns_ref_only += 1
            elif d.has_test:
                self.patterns_test_only += 1
                self.patterns_by_category[d.category] += 1
        for n,d in node_dict.iteritems():
            strand, start, end = n
            length = end - start
            if d.has_ref and d.has_test:
                self.cov_both += length
                self.cov_by_category[d.category] += length
            elif d.has_ref:
                self.cov_ref_only += length
            elif d.has_test:
                self.cov_test_only += length
                self.cov_by_category[d.category] += length

class Match(object):
    def __init__(self):
        self.nodes = collections.defaultdict(lambda: [])
        self.introns = []
        self.splicing = False

class MatchStats(object):
    @staticmethod
    def header_fields():
        return ['transcript_id', 'gene_id', 'locus', 'length', 'num_introns',
                'ref_transcript_id', 'ref_gene_id', 'ref_orig_gene_id', 'ref_gene_name',
                'ref_source', 'ref_gene_type', 'ref_locus',
                'ref_length', 'ref_num_introns',
                'shared_same_strand_bp', 'shared_opp_strand_bp',
                'shared_introns', 'shared_splicing',
                'distance', 'category']

    def __init__(self):
        for field in MatchStats.header_fields():
            setattr(self, field, None)

    def __str__(self):
        fields = []
        for field in MatchStats.header_fields():
            fields.append(getattr(self, field))
        return '\t'.join(map(str,fields))

    def copy(self):
        other = MatchStats()
        for field in MatchStats.header_fields():
            setattr(other, field, getattr(self,field))
        return other

    def add_gtf_attributes(self, feature):
        attrs = ['ref_transcript_id', 'ref_gene_id',
                #'ref_orig_gene_id',
                 'ref_gene_name', 'ref_source', 'ref_gene_type',
                 #'ref_locus',
                 'ref_length', 'ref_num_introns',
                 'shared_same_strand_bp', 'shared_opp_strand_bp',
                 'shared_introns', 'shared_splicing',
                 #'distance',
                 'category']
        for attr in attrs:
            v = getattr(self, attr)
            feature.attrs[attr] = v

    @staticmethod
    def from_transcript(t, ref=None):
        self = MatchStats()
        self.transcript_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
        if GTFAttr.GENE_ID not in t.attrs.keys():
            self.gene_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
        else:
            self.gene_id = t.attrs[GTFAttr.GENE_ID]
        self.locus = '%s:%d-%d[%s]' % (t.chrom, t.start, t.end, strand_int_to_str(t.strand))
        self.length = t.length
        self.num_introns = len(t.exons) - 1
        if ref is not None:
            self.ref_transcript_id = ref.attrs[GTFAttr.TRANSCRIPT_ID]
            self.ref_gene_id = ref.attrs[GTFAttr.GENE_ID]
            self.ref_locus = '%s:%d-%d[%s]' % (ref.chrom, ref.start, ref.end, strand_int_to_str(ref.strand))
            self.ref_length = ref.length
            self.ref_num_introns = len(ref.exons) - 1
            self.ref_orig_gene_id = ref.attrs.get('orig_gene_id', self.ref_gene_id)
            self.ref_source = ref.attrs.get('source', 'NA')
            if 'gene_name' in ref.attrs:
                self.ref_gene_name = ref.attrs['gene_name']
            elif 'transcript_name' in ref.attrs:
                self.ref_gene_name = ref.attrs['transcript_name']
            else:
                self.ref_gene_name = self.ref_gene_id
            if 'gene_type' in ref.attrs:
                self.ref_gene_type = ref.attrs['gene_type']
            elif 'gene_biotype' in ref.attrs:
                self.ref_gene_type = ref.attrs['gene_biotype']
            elif 'transcript_type' in ref.attrs:
                self.ref_gene_type = ref.attrs['transcript_type']
            else:
                self.ref_gene_type = 'None'
        return self

    @staticmethod
    def choose_best(lst, transcript_id_to_source_dict=None):
        hits = []
        for m in lst:
            total_introns = m.num_introns + m.ref_num_introns
            if total_introns == 0:
                intron_frac = 0.0
            else:
                intron_frac = float(m.shared_introns) / (total_introns - m.shared_introns)
            same_strand_frac = float(m.shared_same_strand_bp) / (m.length + m.ref_length - m.shared_same_strand_bp)
            opp_strand_frac = float(m.shared_opp_strand_bp) / (m.length + m.ref_length - m.shared_opp_strand_bp)
            category_int = Category.to_int(m.category)
            hits.append((int(m.shared_splicing), intron_frac,
                         same_strand_frac, opp_strand_frac,
                         int(category_int == Category.INTRONIC_SAME_STRAND),
                         int(category_int == Category.INTRONIC_OPP_STRAND),
                         int(category_int == Category.INTERLEAVING_SAME_STRAND),
                         int(category_int == Category.INTERLEAVING_OPP_STRAND),
                         int(category_int == Category.ENCOMPASSING_SAME_STRAND),
                         int(category_int == Category.ENCOMPASSING_OPP_STRAND),
                         int(category_int == Category.INTERGENIC),
                         -abs(m.distance), m))

        # sort matches
        max_value = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        max_match = None
        if (transcript_id_to_source_dict == None):
            for hit in hits:
                if hit[:12] > max_value:
                    max_value = hit[:12]
                    max_match = hit[-1]
        else:
            # Use ENSEMBL in a tie situation
            max_value = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            max_match = None
            max_match_source = None
            for hit in hits:
                if hit[:12] >= max_value:
                    if (hit[:12] > max_value):
                        max_value = hit[:12]
                        max_match = hit[-1]
                        max_match_source = transcript_id_to_source_dict[hit[-1].ref_transcript_id]
                    elif (max_match != None) and (max_match_source != "ENSEMBL") and (transcript_id_to_source_dict[hit[-1].ref_transcript_id] == "ENSEMBL"):
                        max_value = hit[:12]
                        max_match = hit[-1]
                        max_match_source = "ENSEMBL"


        return max_match

    @staticmethod
    def sort_genome(lst):
        poslst = []
        strands = []
        for m in lst:
            chrom, startend = m.ref_locus[:-3].split(':')
            strands.append(m.ref_locus[-2])
            start, end = map(int, startend.split('-'))
            poslst.append((chrom, start, end, m))
        reverse = any(x == '-' for x in strands)
        poslst.sort(key=operator.itemgetter(0,1,2), reverse=reverse)
        return [x[3] for x in poslst]

    @staticmethod
    def consensus(lst, transcript_id_to_source_dict=None):
        if len(lst) == 0:
            return None
        # first check for read through transcripts involving multiple
        # reference genes
        same_strand_hits = collections.defaultdict(lambda: [])
        for m in lst:
            category_int = Category.to_int(m.category)
            if category_int == Category.SAME_STRAND:
                same_strand_hits[m.ref_gene_id].append(m)
        # no same strand matches so don't need to worry about
        # read-throughs or multiple gene types
        if len(same_strand_hits) == 0:
            return MatchStats.choose_best(lst)
        # get consensus match from same strand overlapping genes
        total_introns = lst[0].num_introns
        total_length = lst[0].length
        shared_introns = 0
        shared_same_strand_bp = 0
        hits = []
        for genelst in same_strand_hits.itervalues():
            m = MatchStats.choose_best(genelst, transcript_id_to_source_dict).copy()
            m.ref_gene_type = ','.join(sorted(set(m.ref_gene_type for m in genelst)))
            total_introns += m.ref_num_introns
            total_length += m.ref_length
            shared_introns += m.shared_introns
            shared_same_strand_bp += m.shared_same_strand_bp
            hits.append(m)
        # sort reference genes by position
        hits = MatchStats.sort_genome(hits)
        # make a new MatchStats object
        hit = hits[0].copy()
        hit.ref_transcript_id = ','.join(x.ref_transcript_id for x in hits)
        hit.ref_gene_id = ','.join(x.ref_gene_id for x in hits)
        hit.ref_orig_gene_id = ','.join(x.ref_orig_gene_id for x in hits)
        hit.ref_gene_name = ','.join(x.ref_gene_name for x in hits)
        hit.ref_source = ','.join(x.ref_source for x in hits)
        hit.ref_gene_type = ','.join(x.ref_gene_type for x in hits)
        hit.ref_locus = ','.join(x.ref_locus for x in hits)
        hit.ref_length = ','.join(str(x.ref_length) for x in hits)
        hit.ref_num_introns = ','.join(str(x.ref_num_introns) for x in hits)
        hit.shared_same_strand_bp = shared_same_strand_bp
        hit.shared_opp_strand_bp = 0
        hit.shared_introns = shared_introns
        hit.shared_splicing = any(m.shared_splicing for m in hits)
        hit.distance = 0
        if len(same_strand_hits) > 1:
            hit.category = Category.to_str(Category.READ_THROUGH)
        return hit


def compare_locus(transcripts):
    # store reference introns
    # (strand,start,end) -> ids (set)
    ref_intron_dict = collections.defaultdict(lambda: [])
    ref_node_dict = collections.defaultdict(lambda: [])
    ref_splicing_patterns = collections.defaultdict(lambda: [])
    ref_dict = {}
    # find the intron domains of the transcripts
    boundaries = find_exon_boundaries(transcripts)
    test_transcripts = []
    for t in transcripts:
        # print 'is_ref', t.attrs[GTFAttr.REF]
        # separate ref and nonref transcripts
        is_ref = bool(int(t.attrs[GTFAttr.REF]))
        if is_ref:
            # add to dict
            ref_id = t.attrs[GTFAttr.TRANSCRIPT_ID]
            ref_dict[ref_id] = t
            # split exons that cross boundaries and get the
            # nodes in the transcript path
            for n in split_exons(t, boundaries):
                ref_node_dict[n].append(t)
            # add to introns
            splicing_pattern = []
            for start,end in t.iterintrons():
                intron = (t.strand, start, end)
                ref_intron_dict[intron].append(t)
                splicing_pattern.append(intron)
            # add to splicing patterns
            if len(splicing_pattern) > 0:
                ref_splicing_patterns[tuple(splicing_pattern)].append(t)
        else:
            test_transcripts.append(t)
    # print test_transcripts
    # index introns for fast intersection
    intron_tree = IntervalTree()
    for intron, refs in ref_intron_dict.iteritems():
        strand, start, end = intron
        intron_tree.insert_interval(Interval(start,end,strand=strand,value=refs))
    # categorize transcripts
    for t in test_transcripts:
        # get transcript nodes and introns
        nodes = list(split_exons(t, boundaries))
        introns = []
        for start,end in t.iterintrons():
            introns.append((t.strand,start,end))
        splicing_pattern = tuple(introns)
        # keep list of all matching ref transcripts
        matches = collections.defaultdict(lambda: Match())
        # dict of reference transcripts -> category -> list of nodes
        for n in nodes:
            if n in ref_node_dict:
                # look for reference transcripts that share this node
                for ref in ref_node_dict[n]:
                    if cmp_strand(t.strand, ref.strand):
                        c = Category.SAME_STRAND
                    else:
                        c = Category.OPP_STRAND
                    ref_id = ref.attrs[GTFAttr.TRANSCRIPT_ID]
                    m = matches[ref_id]
                    m.nodes[c].append(n)
            # look for reference introns that overlap this node
            for hit in intron_tree.find(*n):
                if cmp_strand(t.strand, hit.strand):
                    c = Category.INTRONIC_SAME_STRAND
                else:
                    c = Category.INTRONIC_OPP_STRAND
                for ref in hit.value:
                    ref_id = ref.attrs[GTFAttr.TRANSCRIPT_ID]
                    m = matches[ref_id]
                    m.nodes[c].append(n)
        # dict of introns -> list of reference transcripts
        for intron in introns:
            if intron in ref_intron_dict:
                for ref in ref_intron_dict[intron]:
                    ref_id = ref.attrs[GTFAttr.TRANSCRIPT_ID]
                    m = matches[ref_id]
                    m.introns.append(intron)
        # check splicing pattern matches
        if len(splicing_pattern) > 0:
            if splicing_pattern in ref_splicing_patterns:
                for ref in ref_splicing_patterns[splicing_pattern]:
                    ref_id = ref.attrs[GTFAttr.TRANSCRIPT_ID]
                    m = matches[ref_id]
                    m.splicing = True
        # go through the matches for this transcript and determine
        # the transcript category
        match_stats = []
        for ref_id, m in matches.iteritems():
            ref = ref_dict[ref_id]
            # calculate coverage
            same_strand_bp = sum((n[1] - n[0]) for n in m.nodes[Category.SAME_STRAND])
            opp_strand_bp = sum((n[1] - n[0]) for n in m.nodes[Category.OPP_STRAND])
            # count shared introns
            num_shared_introns = len(m.introns)
            # decide category for this test/ref transcript pair
            if m.splicing or (num_shared_introns > 0) or (same_strand_bp > 0):
                c = Category.SAME_STRAND
            elif (opp_strand_bp > 0):
                c = Category.OPP_STRAND
            else:
                # count nodes of different types
                num_same_strand = len(m.nodes[Category.SAME_STRAND])
                num_opp_strand = len(m.nodes[Category.OPP_STRAND])
                num_intronic_same_strand = len(m.nodes[Category.INTRONIC_SAME_STRAND])
                num_intronic_opp_strand = len(m.nodes[Category.INTRONIC_OPP_STRAND])
                assert num_same_strand == 0
                assert num_opp_strand == 0
                num_intronic = (num_intronic_same_strand +
                                num_intronic_opp_strand)
                assert num_intronic > 0
                if (num_intronic == len(nodes)):
                    # completely intronic
                    if num_intronic_same_strand > 0:
                        c = Category.INTRONIC_SAME_STRAND
                    else:
                        c = Category.INTRONIC_OPP_STRAND
                else:
                    # interleaving means some nodes intronic and other intergenic
                    if num_intronic_same_strand > 0:
                        c = Category.INTERLEAVING_SAME_STRAND
                    else:
                        c = Category.INTERLEAVING_OPP_STRAND
            # create a match object
            ms = MatchStats.from_transcript(t, ref)
            ms.shared_same_strand_bp = same_strand_bp
            ms.shared_opp_strand_bp = opp_strand_bp
            ms.shared_introns = num_shared_introns
            ms.shared_splicing = m.splicing
            ms.category = Category.to_str(c)
            ms.distance = 0
            match_stats.append(ms)
        yield (t, match_stats)

def build_locus_trees(gtf_file):
    transcripts = []
    locus_cluster_trees = collections.defaultdict(lambda: ClusterTree(0,1))
    for locus_transcripts in parse_gtf(open(gtf_file)):
        for t in locus_transcripts:
            is_ref = bool(int(t.attrs[GTFAttr.REF]))
            if not is_ref:
                continue
            i = len(transcripts)
            transcripts.append(t)
            locus_cluster_trees[t.chrom].insert(t.start, t.end, i)
    # build interval trees of loci
    locus_trees = collections.defaultdict(lambda: IntervalTree())
    for chrom, cluster_tree in locus_cluster_trees.iteritems():
        for locus_start, locus_end, indexes in cluster_tree.getregions():
            for i in indexes:
                locus_transcripts = [transcripts[i] for i in indexes]
                locus_trees[chrom].insert_interval(Interval(locus_start, locus_end, value=locus_transcripts))
    return locus_trees

def find_nearest_transcripts(chrom, start, end, strand, locus_trees):
    # first check for overlap
    nearest_features = []
    hits = locus_trees[chrom].find(start, end)
    for hit in hits:
        for t in hit.value:
            if cmp_strand(t.strand, strand):
                c = Category.ENCOMPASSING_SAME_STRAND
            else:
                c = Category.ENCOMPASSING_OPP_STRAND
            nearest_features.append((t, c, 0))
    # look left and right
    left_hits = locus_trees[chrom].before(start, num_intervals=1, max_dist=MAX_LOCUS_DIST)
    right_hits = locus_trees[chrom].after(end, num_intervals=1, max_dist=MAX_LOCUS_DIST)
    # look for nearest hit
    for hits in (left_hits, right_hits):
        nearest_locus_hit = None
        nearest_dist = MAX_LOCUS_DIST
        for hit in hits:
            dist = min(abs(start - hit.end), abs(hit.start - end))
            if dist < nearest_dist:
                nearest_dist = dist
                nearest_locus_hit = hit
        if nearest_locus_hit is not None:
            for t in nearest_locus_hit.value:
                dist = min(abs(start - t.end), abs(t.start - end))
                nearest_features.append((t, Category.INTERGENIC, dist))
    return nearest_features


def _parse_gtf_by_chrom(gtf_file):
    current_chrom = None
    exon_dict = collections.defaultdict(lambda: [])
    transcript_dict = {}
    for feature in GTFFeature.parse(open(gtf_file)):
        if (feature.feature_type != "transcript") and (feature.feature_type != "exon"):
            continue
        if (current_chrom != feature.seqid):
            if len(exon_dict) > 0:
                yield current_chrom, transcript_dict, exon_dict
                exon_dict = collections.defaultdict(lambda: [])
                transcript_dict = {}
            current_chrom = feature.seqid
        t_id = feature.attrs[GTFAttr.TRANSCRIPT_ID]
        if feature.feature_type == "transcript":
            transcript_dict[t_id] = feature
        elif feature.feature_type == "exon":
            exon_dict[t_id].append(feature)
    if len(exon_dict) > 0:
        yield current_chrom, transcript_dict, exon_dict


def add_gtf_file(gtf_file, outfh, is_ref):
    refval = '1' if is_ref else '0'
    for chrom, transcript_dict, exon_dict in _parse_gtf_by_chrom(gtf_file):
        logging.debug("\tfinished chrom %s %d features" % (chrom, len(exon_dict)))
        # output reference transcripts
        for t_id, features in exon_dict.iteritems():
            # sort features (exons) by start position
            features.sort(key=operator.attrgetter('start'))
            # annotate exons as reference features
            for f in features:
                f.attrs[GTFAttr.REF] = refval
                print >>outfh, str(f)
            # transcript feature
            if t_id in transcript_dict:
                f = transcript_dict[t_id]
            else:
                f = GTFFeature()
                f.seqid = features[0].seqid
                f.source = features[0].source
                f.feature_type = 'transcript'
                f.start = features[0].start
                f.end = features[-1].end
                f.score = features[0].score
                f.strand = features[0].strand
                f.phase = '.'
                f.attrs = features[0].attrs.copy()
                if "exon_number" in f.attrs:
                    del f.attrs["exon_number"]
            f.attrs[GTFAttr.REF] = refval
            print >>outfh, str(f)

def impute_transcript_type(catint, length, gene_type, ref_gene_type):
    if (catint == Category.SAME_STRAND or
        catint == Category.READ_THROUGH):
        # impute gene type
        transcript_type = ref_gene_type
    else:
        if gene_type == 'protein_coding':
            # don't change protein coding genes
            transcript_type = gene_type
        elif length < 250:
            # categorize small RNA separately
            transcript_type = 'misc_RNA'
        else:
            if ref_gene_type == 'protein_coding':
                # categorize based on overlap with reference
                transcript_type = PROTEIN_CATEGORY_MAP[catint]
            else:
                # reference is also non-coding
                transcript_type = 'lincRNA'
    return transcript_type

def parse_and_output_gtf_file(input_gtf_file, tmp_dir, is_ref, num_cores):
    parallel_sort_cmd = False
    with open(os.devnull, "w") as fnull:
        cmdline = 'echo "2 1" | %s --parallel=2'
        if subprocess.call(cmdline % 'gsort', stdout=fnull, stderr=fnull, shell=True) == 0:
            parallel_sort_cmd = 'gsort'
        if subprocess.call(cmdline % 'sort', stdout=fnull, stderr=fnull, shell=True) == 0:
            parallel_sort_cmd = 'sort'

    if not parallel_sort_cmd:
        logging.warning('Command line "sort" command does not support '
                        '--parallel flag. For improved performance, consider '
                        'upgrading/installing the latest GNU coreutils to '
                        'enable parallel sort.')
        args = ["sort"]
    else:
        logging.debug('Command line "sort" supports --parallel flag')
        args = [parallel_sort_cmd, '--parallel=%d' % num_cores]

    if tmp_dir is not None:
        args.extend(["-T", tmp_dir])
    args.extend(["-k1,1", input_gtf_file])
    myenv = os.environ.copy()
    myenv["LC_ALL"] = "C"
    subprocess.call(args, stdout=open(os.path.join(tmp_dir, "input." + str(is_ref) + ".srt.gtf"), "w"), env=myenv)

    curr_seqid = None
    outfh = None
    with open(os.path.join(tmp_dir, "input." + str(is_ref) + ".srt.gtf"), "r") as sorted_input_file:
        for input_line in sorted_input_file:
            if (input_line[0] == '#'):
                continue

            seqname, source, feature, start, end, score, strand, frame, attribute = input_line.replace("\n", "").split("\t")

            seqname = seqname.translate(None, string.punctuation).replace(" ", "")

            if (curr_seqid == seqname):
                outfh.write(input_line)
            else:
                try:
                    outfh.close()
                except:
                    pass
                seq_id_folder = os.path.join(tmp_dir, str(seqname))
                if not os.path.exists(seq_id_folder):
                    logging.debug("Creating tmp seqid directory '%s'" % (seq_id_folder))
                    os.makedirs(seq_id_folder)

                outfh = open(os.path.join(seq_id_folder, "input." + str(is_ref) + ".srt.gtf"), "w")
                outfh.write(input_line)
                curr_seqid = seqname


def compare_assemblies_worker(input_queue):
    while True:
        output_dir = input_queue.get()

        if (output_dir == None):
            input_queue.task_done()
            break

        ref_file = os.path.join(output_dir, "input.1.srt.gtf")
        test_file = os.path.join(output_dir, "input.0.srt.gtf")

        if not (os.path.isfile(ref_file) and os.path.isfile(test_file)):
            logging.info("Skipping: " + os.path.basename(output_dir) +
                " because reference and test have 0 overlap")
            input_queue.task_done()
            continue

        # merge step
        merged_gtf_file = os.path.join(output_dir, "merged.gtf")
        merged_sorted_gtf_file = os.path.splitext(merged_gtf_file)[0] + ".srt.gtf"
        merge_done_file = os.path.join(output_dir, 'merged.done')
        sort_done_file = os.path.join(output_dir, 'sort.done')
        if not os.path.exists(merge_done_file):
            # merge and sort ref/test gtf files
            logging.info("Merging reference and test GTF files")
            # make temporary file to store merged ref/test gtf files
            with open(merged_gtf_file, "w") as fileh:
                logging.info("Adding reference GTF file")
                add_gtf_file(ref_file, fileh, True)
                logging.info("Adding test GTF file")
                add_gtf_file(test_file, fileh, False)
            open(merge_done_file, 'w').close()

        if not os.path.exists(sort_done_file):
            logging.info("Sorting merged GTF file")
            # create temp directory
            tmp_dir = os.path.join(output_dir, 'tmp')
            if not os.path.exists(tmp_dir):
                logging.debug("Creating tmp directory '%s'" % (tmp_dir))
                os.makedirs(tmp_dir)
            sort_gtf(merged_gtf_file, merged_sorted_gtf_file, tmp_dir)
            # cleanup
            shutil.rmtree(tmp_dir)
            open(sort_done_file, 'w').close()

        # generate transcript_id to source dict
        transcript_id_to_source_dict = {}
        for feature in GTFFeature.parse(open(merged_sorted_gtf_file, "r")):
            transcript_id_to_source_dict[feature.attrs['transcript_id']] = feature.source

        # compare assemblies
        overlapping_gtf_file = os.path.join(output_dir, 'overlapping.gtf')
        intergenic_tmp_gtf_file = os.path.join(output_dir, 'intergenic.tmp.gtf')
        overlapping_file = os.path.join(output_dir, 'overlapping.tsv')
        # overlapping_consensus_file = os.path.join(output_dir, 'overlapping.consensus.tsv')
        overlapping_done_file = os.path.join(output_dir, 'overlapping.done')
        stats_file = os.path.join(output_dir, 'stats.txt')
        stats_obj = GlobalStats()
        num_intergenic = 0
        if not os.path.exists(overlapping_done_file):
            logging.info("Comparing assemblies")
            gtf_fileh = open(overlapping_gtf_file, 'w')
            tmp_gtf_fileh = open(intergenic_tmp_gtf_file, 'w')
            # overlapping_fileh = open(overlapping_file, 'w')
            # overlapping_consensus_fileh = open(overlapping_consensus_file, 'w')
            for locus_transcripts in parse_gtf(open(merged_sorted_gtf_file)):
                locus_chrom = locus_transcripts[0].chrom
                locus_start = locus_transcripts[0].start
                locus_end = max(t.end for t in locus_transcripts)
                logging.debug("[LOCUS] %s:%d-%d %d transcripts" %
                              (locus_chrom, locus_start, locus_end,
                               len(locus_transcripts)))
                for t, match_stats in compare_locus(locus_transcripts):
                    if len(match_stats) == 0:
                        # write intergenic transcripts to analyze separately
                        t.attrs['category'] = Category.to_str(Category.INTERGENIC)
                        for f in t.to_gtf_features(source='assembly'):
                            print >>tmp_gtf_fileh, str(f)
                        num_intergenic += 1
                    else:
                        # get consensus match information
                        consensus_match = MatchStats.consensus(match_stats, transcript_id_to_source_dict)
                        assert consensus_match is not None
                        t.attrs['category'] = consensus_match.category
                        # add gtf attributes and write
                        for f in t.to_gtf_features(source='assembly'):
                            if t.attrs['category'] in ['same_strand', 'read_through']:
                                consensus_match.add_gtf_attributes(f)
                            print >>gtf_fileh, str(f)
                        # tab-delimited text output
                        # print >>overlapping_consensus_fileh, str(consensus_match)
                        # for ms in match_stats:
                            # print >>overlapping_fileh, str(ms)
                # compute global statistics
                stats_obj.compute(locus_transcripts)
            logging.debug("Reporting global statistics")
            with open(stats_file, 'w') as f:
                print >>f, stats_obj.report()
            gtf_fileh.close()
            tmp_gtf_fileh.close()
            # overlapping_fileh.close()
            # overlapping_consensus_fileh.close()
            open(overlapping_done_file, 'w').close()

        # resolve intergenic transcripts
        intergenic_gtf_file = os.path.join(output_dir, 'intergenic.gtf')
        intergenic_file = os.path.join(output_dir, 'intergenic.tsv')
        intergenic_best_file = os.path.join(output_dir, 'intergenic.best.tsv')
        intergenic_done_file = os.path.join(output_dir, 'intergenic.done')
        if not os.path.exists(intergenic_done_file):
            logging.debug("Characterizing transcripts with complex overlap")
            locus_trees = build_locus_trees(merged_sorted_gtf_file)
            logging.debug('Finding nearest matches to intergenic transcripts')
            gtf_fileh = open(intergenic_gtf_file, 'w')
            # intergenic_fileh = open(intergenic_file, 'w')
            intergenic_best_fileh = open(intergenic_best_file, 'w')
            def wc(infile):
                p = subprocess.Popen('wc -l %s' % infile, shell=True, stdout=subprocess.PIPE)
                out, err = p.communicate()
                return int(out.split()[0])
            if wc(intergenic_tmp_gtf_file) != 0:
                for locus_transcripts in parse_gtf(open(intergenic_tmp_gtf_file)):
                    for t in locus_transcripts:
                        # find nearest transcripts
                        nearest_transcripts = find_nearest_transcripts(t.chrom, t.start, t.end, t.strand, locus_trees)
                        match_stats = []
                        best_match = None
                        if len(nearest_transcripts) == 0:
                            best_match = MatchStats.from_transcript(t)
                            best_match.category = Category.to_str(Category.INTERGENIC)
                            match_stats.append(best_match)
                        else:
                            for ref,category,dist in nearest_transcripts:
                                # create a match object
                                ms = MatchStats.from_transcript(t, ref)
                                ms.shared_same_strand_bp = 0
                                ms.shared_opp_strand_bp = 0
                                ms.shared_introns = 0
                                ms.shared_splicing = False
                                ms.category = Category.to_str(category)
                                ms.distance = dist
                                match_stats.append(ms)
                            # choose the consensus match
                            best_match = MatchStats.choose_best(match_stats)
                        # add gtf attributes and write
                        for f in t.to_gtf_features(source='assembly'):
                            # best_match.add_gtf_attributes(f)
                            print >>gtf_fileh, str(f)
                        # write tab-delimited data
                        # print >>intergenic_best_fileh, str(best_match)
                        # for ms in match_stats:
                            # print >>intergenic_fileh, str(ms)
                gtf_fileh.close()
                # intergenic_fileh.close()
                intergenic_best_fileh.close()
                open(intergenic_done_file, 'w').close()
        # merge overlapping and intergenic results
        metadata_file = os.path.join(output_dir, 'metadata.txt')
        metadata_consensus_file = os.path.join(output_dir, 'metadata.consensus.txt')
        assembly_gtf_file = os.path.join(output_dir, 'assembly.cmp.gtf')
        combine_done_file = os.path.join(output_dir, 'combine.done')
        if not os.path.exists(combine_done_file):
            logging.debug('Merging results')
            # filenames = [overlapping_file, intergenic_file]
            # with open(metadata_file, 'w') as outfile:
            #     print >>outfile, '\t'.join(MatchStats.header_fields())
            #     for fname in filenames:
            #         with open(fname) as infile:
            #             for line in infile:
            #                 outfile.write(line)
            filenames = [intergenic_best_file]
            with open(metadata_consensus_file, 'w') as outfile:
                print >>outfile, '\t'.join(MatchStats.header_fields())
                for fname in filenames:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            filenames = [intergenic_gtf_file, overlapping_gtf_file]
            with open(assembly_gtf_file, 'w') as outfile:
                for fname in filenames:
                    with open(fname) as infile:
                        for line in infile:
                            outfile.write(line)
            open(combine_done_file, 'w').close()

        input_queue.task_done()

def compare_assemblies(ref_gtf_file, test_gtf_file, output_dir, output_final, cpat, num_cores):
    # output files
    if not os.path.exists(output_dir):
        logging.debug('Creating tmp dir: %s' % (output_dir))
        os.makedirs(output_dir)

    # split input GTFs into chromosome-large GTFs
    parse_and_output_gtf_file(ref_gtf_file, output_dir, 1, num_cores)
    parse_and_output_gtf_file(test_gtf_file, output_dir, 0, num_cores)


    input_queue = multiprocessing.JoinableQueue(maxsize = num_cores)
    procs = []

    for i in xrange(num_cores):
        p = multiprocessing.Process(target=compare_assemblies_worker, args=(input_queue, ))
        p.start()
        procs.append(p)

    for folder in [x[0] for x in os.walk(output_dir) if x[0] != output_dir]:
        input_queue.put(folder)

    for element in [None] * num_cores:
        input_queue.put(element)

    # close input queue
    input_queue.join()
    input_queue.close()

    # join worker processes
    for p in procs:
        p.join()


    # Merge all the chromosomes
    assembly_gtf_file = os.path.join(output_dir, "assembly.cmp.gtf")

    with open(assembly_gtf_file, "w") as outfh:
        for folder in [x[0] for x in os.walk(output_dir)]:
            assembled_chrom_gtf_file = os.path.join(folder, "assembly.cmp.gtf")
            if not os.path.isfile(assembled_chrom_gtf_file):
                continue

            with open(assembled_chrom_gtf_file, "r") as inputfh:
                for line in inputfh:
                    outfh.write(line)


    if wc(assembly_gtf_file) == 0:
        logging.error('Zero overlap for reference and test GTFs. '
        'Ensure they are from the same species / genome.')
        sys.exit()


    # read compared assembly and add annotation status / final category
    if cpat['run']:
        logging.info('Running coding potential prediction')
        #make bed file
        assembly_bed_cpat_file = os.path.abspath(os.path.join(output_dir, 'assembly.cpat.bed'))
        logging.debug('Converting GTF to BED for CPAT')
        with open(assembly_bed_cpat_file, 'w') as f:
            for transcripts in parse_gtf(open(assembly_gtf_file)):
                for t in transcripts:
                    if 'gene_id' in t.attrs.keys():
                        name = '%s|%s' % (t.attrs['gene_id'], t.attrs['transcript_id'])
                    else:
                        name = t.attrs['transcript_id']
                    fields = write_bed(t.chrom, name, t.strand, 1000, t.exons)
                    print >>f, '\t'.join(fields)

        cpat_tsv_file = os.path.abspath(os.path.join(output_dir, 'cpat.tsv'))

        logging.info('Running CPAT')
        cmd = [
            cpat['exec'],
            '-g', assembly_bed_cpat_file,
            '-o', cpat_tsv_file,
            '-x', cpat['hex'],
            '-d', cpat['model'],
            '-r', cpat['genome']
        ]
        subprocess.call(' '.join(cmd), shell=True, cwd=output_dir)

        cpat_dict = {}
        fileh = open(cpat_tsv_file)
        header = fileh.next().strip().split('\t')
        tid_idx = 0
        orf_idx = 2
        prob_idx = 5
        for line in fileh:
            line = line.strip().split('\t')
            tid = line[tid_idx].split('|')[1]
            orf = line[orf_idx]
            prob = line[prob_idx]
            cpat_dict[tid] = (orf, prob)


    #### RECATEGORIZE THE GENCODE ANNOTATION CATEOGRY #####
    # ANNOTATED = ['same_strand', 'read_through']
    # UNANNOTATED = ['opp_strand', 'intronic_same_strand', 'intronic_opp_strand',
    #     'interleaving_same_strand', 'interleaving_opp_strand',
    #     'encompassing_opp_strand', 'encompassing_same_strand',
    #     'intergenic']
    assembly_refcomp_file = os.path.join(output_final, 'assembly.refcomp.gtf')
    final_gtf_done_file = os.path.join(output_dir, 'finalgtf.done')
    if not os.path.exists(final_gtf_done_file):
        logging.info('Generating final GTF file')
        with open(assembly_refcomp_file, 'w') as gtf_final:
            for locus_transcripts in parse_gtf(open(assembly_gtf_file)):
                for t in locus_transcripts:
                    tid = t.attrs['transcript_id']
                    if cpat['run']:
                        orf, prob = cpat_dict[tid]
                        t.attrs['orf_size'] = orf
                        t.attrs['cpat_coding_prob'] = prob
                    catstr = t.attrs['category']
                    catint = Category.to_int(catstr)
                    length = t.length
                    gene_type = t.attrs.get('gene_type', None)
                    cat = t.attrs['category']
                    t.attrs['category_relative_detail'] = cat
                    if cat in ['same_strand', 'read_through']:
                        cat_rel = 'exonic_overlap'
                    elif cat == 'intergenic':
                        cat_rel = 'intergenic'
                    else:
                        cat_rel = 'intragenic'
                    t.attrs['category_relative'] = cat_rel
                    if 'ref_gene_type' in t.attrs.keys():
                        ref_gene_name = ','.join(set(t.attrs['ref_gene_name'].split(',')))
                        t.attrs['ref_gene_name'] = ref_gene_name
                        ref_gene_type = t.attrs['ref_gene_type']
                        ref_gene_types = set(ref_gene_type.split(','))
                        transcript_types = set(impute_transcript_type(catint, length, gene_type, x) for x in ref_gene_types)
                        t.attrs['ref_gene_type'] = ','.join(transcript_types)
                        transcript_categories = set(gencode_category_map(x) for x in transcript_types)
                        # sorted and join unique types/categories to make conglomerated category assignments
                        transcript_type = ','.join(sorted(transcript_types))
                        transcript_category = ','.join(sorted(transcript_categories))
                        # ref_cat_transform = GENCODE_CATEGORY_MAP[ref_gene]
                        transcript_category_final = transcript_category
                        if transcript_category in ['lncRNA', 'ncRNA', 'lncRNA,ncRNA']:
                            transcript_category_final = 'lncrna'
                        elif (',' in transcript_category) & (('protein' in transcript_category) | ('pseudo' in transcript_category)):
                            transcript_category_final = 'mixed_read_through'
                        t.attrs['annotation'] = 'annotated'
                    else:
                        if cpat['run']:
                            if float(prob) > cpat['cutoff']:
                                    transcript_category_final = 'tucp'
                            else:
                                transcript_category_final = 'lncrna'
                        else:
                            transcript_category_final = 'lncrna'
                        t.attrs['annotation'] = 'unannotated'

                    t.attrs['category'] = transcript_category_final
                    if 'ref' in t.attrs.keys():
                        del t.attrs['ref']
                    if 'ref_num_introns' in t.attrs.keys():
                        del t.attrs['ref_num_introns']
                    if 'ref_source' in t.attrs.keys():
                        del t.attrs['ref_source']
                    for f in t.to_gtf_features(source='assembly'):
                        print >>gtf_final, str(f)
        open(final_gtf_done_file, 'w').close()

    logging.info('Generating metadata TSV file')
    gtf_attrs = FULL_GTF_ATTRS
    logging.debug("Reading GTF attributes to make metadata file")
    metadata_dict = get_gtf_metadata(assembly_refcomp_file, gtf_attrs)
    header_fields = ['transcript_id', 'chrom', 'start', 'end', 'strand',
                     'num_exons', 'transcript_length'] + gtf_attrs
    assembly_metadata_file = os.path.join(output_final, 'assembly.metadata.tsv')
    with open(assembly_metadata_file, 'w') as meta_final:
        print >>meta_final, '\t'.join(header_fields)
        for t_id in sorted(metadata_dict):
            m = metadata_dict[t_id]
            fields = [t_id, m.chrom, m.start, m.end, m.strand, m.num_exons, m.length]
            for attr in gtf_attrs:
                fields.append(getattr(m, attr))
            print >>meta_final, '\t'.join(map(str,fields))
    logging.info("Done")


def cpat_init(cpat_bool, genome, species):
    #determine whether or not script is being run from a packaged executable
    #this is used to identify where the data files are
    if getattr( sys, 'frozen', False ) :
        datadir = os.path.join(sys._MEIPASS, 'share','data')
    else:
        datadir = os.path.join(os.path.dirname(share.__file__), 'data')

    #select OS to identify which CPAT executable to use
    if sys.platform == 'darwin':
        cpat_exec = os.path.join(datadir, 'cpat_execs', 'cpat_exec_mac')
    elif (sys.platform == 'win32'):
        raise OSError("Error: Windows is not supported")
        exit(1)
    else:
        cpat_exec = os.path.join(datadir, 'cpat_execs', 'cpat_exec_linux')

    hexamer = os.path.join(datadir, 'cpat_refs', '%s/hexamer.tab' % species)
    model = os.path.join(datadir, 'cpat_refs', '%s/logitmodel.RData' % species)
    if species == 'human':
        cutoff = 0.5
    else:
        cutoff_file = get_data(os.path.join('cpat_refs', '%s/logitmodel.RData' % species))
        cutoff = float(open(cutoff_file).next().strip().split(' :')[1])

    cpat_dict = {
        'run': cpat_bool,
        'exec': cpat_exec,
        'hex': hexamer,
        'model': model,
        'cutoff': cutoff,
        'genome': genome
    }

    return cpat_dict




def main():
    # parse command line
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--verbose", action="store_true",
                        dest="verbose", default=False)
    parser.add_argument("-o", "--output-dir", dest="output_dir",
                        default="taco_compare",
                        help='Directory for reference comparison output')
    parser.add_argument("-p", "--num-processes", type=int, default=1,
                        dest="num_cores", help='Run tool in parallel with N processes. '
                        '(note: each core processes 1 chromosome) ')
    # parser.add_argument("--cpat", action='store_true', default=False,
    #                     help='Run CPAT tool to for coding potential scoring. '
    #                     '(CPAT function currently only supports '
    #                     'Human, Mouse, and Zebrafish) '
    #                     '(WARNING: The CPAT tool can take over an hour) ')
    # parser.add_argument("--cpat-species", dest='cpat_spec', default='human',
    #                     help='Select either: human, mouse, zebrafish')
    # parser.add_argument("--cpat-genome", dest='cpat_gen',
    #                     help='Provide a genome fasta for the genome used to '
    #                     'produce assemblies being compared. Required '
    #                     'if \"--cpat\" used. CPAT uses this '
    #                     'to obtain sequence for the provided transcripts')
    parser.add_argument("-r", "--ref-gtf", dest='ref_gtf_file',
                        help='Reference GTF file to compare against')
    parser.add_argument("-t", "--test-gtf", dest='test_gtf_file',
                        help='GTF file used for comparison')
    args = parser.parse_args()
    # set logging level
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO


    #establish CPAT associated files
    cpat_dict = cpat_init(args.cpat, args.cpat_gen, args.cpat_spec)


    if not args.ref_gtf_file or not args.test_gtf_file:
        logging.error('Please provide both reference and test GTF files')
        return 1

    if args.cpat:
        if not args.cpat_gen:
            logging.error('A genome FASTA must be provided (with "--cpat-genome" flag) when using "--cpat" flag')
            return 1

    logging.basicConfig(level=level,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    logging.info("TACO Reference Comparison Tool" )
    logging.info("----------------------------------")
    # show parameters
    logging.info("Parameters:")
    logging.info("verbose logging:       %s" % (args.verbose))
    logging.info("number of cores:       %s" % (args.num_cores))
    logging.info("reference gtf file:    %s" % (args.ref_gtf_file))
    logging.info("test gtf file:         %s" % (args.test_gtf_file))
    logging.info("output directory:      %s" % (args.output_dir))
    logging.info("run cpat:              %s" % (args.cpat))
    logging.info("----------------------------------")

    # check command line parameters
    if not os.path.exists(args.ref_gtf_file):
        parser.error("Reference GTF file %s not found" % (args.ref_gtf_file))
    if not os.path.exists(args.test_gtf_file):
        parser.error("Test GTF file %s not found" % (args.test_gtf_file))
    if not os.path.exists(args.output_dir):
        logging.info("Creating output directory '%s'" % (args.output_dir))
        os.makedirs(args.output_dir)

    tmp_dir = os.path.join(args.output_dir, 'TMP')
    if not os.path.exists(args.output_dir):
        os.makedirs(tmp_dir)
    compare_assemblies(args.ref_gtf_file, args.test_gtf_file,
                       tmp_dir, args.output_dir, cpat_dict, args.num_cores)
    shutil.rmtree(tmp_dir)
    return 0

if __name__ == "__main__":
    sys.exit(main())
