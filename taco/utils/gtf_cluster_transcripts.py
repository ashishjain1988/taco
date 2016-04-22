import os
import sys
import argparse
import logging
import itertools
import operator
import collections

from taco.lib.gtf import GTF
from taco.lib.bx.cluster import ClusterTree


class Transcript(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', 'exons', 'attrs')

    def __init__(self):
        self.chrom = None
        self.start = None
        self.end = None
        self.strand = None
        self.attrs = {}
        self.exons = []

    def to_gtf_features(self):
        # transcript feature
        f = GTF.Feature()
        f.seqid = self.chrom
        f.source = 'taco'
        f.feature = 'transcript'
        f.start = self.start
        f.end = self.end
        f.score = 0
        f.strand = self.strand
        f.phase = '.'
        f.attrs = self.attrs.copy()
        features = [f]
        # exon features
        for e in self.exons:
            start,end = e
            f = GTF.Feature()
            f.seqid = self.chrom
            f.source = 'taco'
            f.feature = 'exon'
            f.start = start
            f.end = end
            f.score = 0
            f.strand = self.strand
            f.phase = '.'
            f.attrs = {}
            f.attrs['transcript_id'] = self.attrs['transcript_id']
            f.attrs['gene_id'] = self.attrs['gene_id']
            features.append(f)
        return features

    @staticmethod
    def create(t, exons):
        self = Transcript()
        self.attrs.update(t.attrs)
        self.chrom = t.seqid
        self.strand = t.strand
        self.start = min(f.start for f in exons)
        self.end = max(f.end for f in exons)
        self.exons = []
        for f in exons:
            self.exons.append((f.start, f.end))
        self.exons.sort()
        return self

    @staticmethod
    def parse_gtf(filename):
        # read all transcripts
        t_dict = {}
        e_dict = collections.defaultdict(lambda: [])
        i = 0
        for f in GTF.parse(open(filename)):
            i += 1
            if (i % 100000) == 0:
                logging.debug('\tread %d lines' % (i))
            if f.feature == 'transcript':
                t_id = f.attrs['transcript_id']
                t_dict[t_id] = f
            elif f.feature == 'exon':
                t_id = f.attrs['transcript_id']
                e_dict[t_id].append(f)
        i = 0
        for t_id, t_feature in t_dict.iteritems():
            exon_features = e_dict[t_id]
            yield Transcript.create(t_feature, exon_features)
            i += 1
        logging.debug('Parsed %d transcripts' % (i))


def cluster_transcripts(gtf_file):
    # read all features
    chrom_feature_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    logging.debug('Parsing gtf file: %s' % (gtf_file))
    for f in Transcript.parse_gtf(gtf_file):
        # bin by chromosome and strand
        chrom_feature_dict[f.chrom][f.strand].append(f)
    # cluster transcripts into genes
    logging.debug('Clustering transcripts into genes')
    cur_gene_id = 1
    for strand_feature_dict in chrom_feature_dict.itervalues():
        for strand_features in strand_feature_dict.itervalues():
            # initialize each transcript to be in a 'gene' by itself
            cluster_map = {}
            cluster_tree = ClusterTree(0,1)
            for i, f in enumerate(strand_features):
                cluster_map[i] = set((i,))
                for start,end in f.exons:
                    cluster_tree.insert(start, end, i)
            for start, end, indexes in cluster_tree.getregions():
                # group transcripts into larger clusters
                new_cluster = set()
                for i in indexes:
                    new_cluster.update(cluster_map[i])
                # reassign transcript clusters to new cluster
                for i in new_cluster:
                    cluster_map[i] = new_cluster
            del cluster_tree
            # now all transcripts are assigned to a gene cluster
            # enumerate all gene clusters
            clusters = set()
            for clust in cluster_map.itervalues():
                clusters.add(frozenset(clust))
            del cluster_map
            # now assign gene ids to each cluster
            for clust in clusters:
                new_gene_id = 'G%011d' % (cur_gene_id)
                for i in clust:
                    f = strand_features[i]
                    f.attrs['orig_gene_id'] = f.attrs['gene_id']
                    f.attrs['gene_id'] = new_gene_id
                cur_gene_id += 1
    # output genes
    logging.debug('Writing transcripts')
    for chrom in sorted(chrom_feature_dict):
        strand_feature_dict = chrom_feature_dict[chrom]
        features = []
        for strand_features in strand_feature_dict.itervalues():
            features.extend(strand_features)
        features.sort(key=operator.attrgetter('start'))
        for f in features:
            for gtf_feature in f.to_gtf_features():
                yield str(gtf_feature)


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    gtf_file = args.gtf_file
    gtf_files = []
    logging.info('gtf file: %s' % (gtf_file))
    for f in cluster_transcripts(gtf_file):
        print f
    logging.info("Done")
    return 0


if __name__ == '__main__':
    sys.exit(main())
