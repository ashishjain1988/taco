import os
import sys
import argparse
import logging
import itertools
import operator
from collections import defaultdict

from taco.lib.gtf import GTF
from taco.lib.base import Strand
from taco.lib.transfrag import Transfrag
from taco.lib.bx.intersection import Interval, IntervalTree
from taco.lib.bx.cluster import ClusterTree


def get_gtf_features(t, gene_id_str):
    strand_str = Strand.to_gtf(t.strand)
    f = GTF.Feature()
    f.seqid = t.chrom
    f.source = 'taco'
    f.feature = 'transcript'
    f.start = t.start
    f.end = t.end
    f.score = 0.0
    f.strand = strand_str
    f.phase = '.'
    f.attrs = {'transcript_id': t._id,
               'gene_id': gene_id_str}

    yield f
    for e in t.exons:
        f = GTF.Feature()
        f.seqid = t.chrom
        f.source = 'taco'
        f.feature = 'exon'
        f.start = e.start
        f.end = e.end
        f.score = 0.0
        f.strand = strand_str
        f.phase = '.'
        f.attrs = {'transcript_id': t._id,
                   'gene_id': gene_id_str}
        yield f


def gtf_cluster_transcripts(gtf_file):
    # read all features and create transfrags
    logging.debug('Parsing GTF file')
    with open(gtf_file) as fh:
        transcripts = Transfrag.parse_gtf(fh).values()
    def sort_key_transfrag(t):
        return (t.chrom, t.start)
    transcripts.sort(key=sort_key_transfrag)
    # cluster by chromosome and strand
    logging.debug('Clustering transcripts on each chromosome strand')
    chrom_strand_cluster_trees = defaultdict(lambda: defaultdict(lambda: ClusterTree(0,1)))
    for i, t in enumerate(transcripts):
        cluster_tree = chrom_strand_cluster_trees[t.chrom][t.strand]
        for start, end in t.exons:
            cluster_tree.insert(start, end, i)
    # assign gene ids
    logging.debug('Assigning gene IDs')
    cur_gene_id = 1
    t_gene_map = {}
    for chrom in sorted(chrom_strand_cluster_trees):
        strand_cluster_trees = chrom_strand_cluster_trees[chrom]
        for strand in sorted(strand_cluster_trees):
            cluster_tree = strand_cluster_trees[strand]
            for start, end, indexes in cluster_tree.getregions():
                indexes = set(indexes)
                for i in indexes:
                    t = transcripts[i]
                    t_gene_map[t._id] = cur_gene_id
                cur_gene_id += 1
    for t in transcripts:
        gene_id_str = 'G%011d' % (t_gene_map[t._id])
        for f in get_gtf_features(t, gene_id_str):
            yield f


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    if not os.path.exists(args.gtf_file):
        parser.error('GTF file "%s" not found' % (args.gtf_file))

    logging.info("Clustering transcripts and assigning gene IDs")
    for f in gtf_cluster_transcripts(args.gtf_file):
        print str(f)
    logging.info("Done")
    return 0


if __name__ == '__main__':
    sys.exit(main())
