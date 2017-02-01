'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
'''
import os
import shutil
import logging
from operator import itemgetter
from multiprocessing import Process, JoinableQueue, Value, Lock
from collections import namedtuple

from taco.lib.gtf import GTF
from taco.lib.base import Strand, Results
from taco.lib.batch_sort import batch_sort, batch_merge, sort_key_gtf, \
    sort_key_bed, SORT_BUFFER_SIZE
from taco.lib.transfrag import Transfrag
from taco.lib.locus import Locus
from taco.lib.cpathfinder import find_paths
from taco.lib.path_graph import PathGraphFactory

__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2017"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "GPL"
__version__ = "0.6.2"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


LocusIndex = namedtuple('LocusIndex', ['name', 'chrom', 'start', 'end',
                        'filepos', 'num_lines'])


class Isoform(object):
    __slots__ = ('expr', 'exons', 'rel_frac', 'abs_frac',
                 'gene_id', 'tss_id')

    def __init__(self, exons=None, expr=0.0, rel_frac=1.0, abs_frac=1.0):
        self.exons = exons
        self.expr = expr
        self.rel_frac = rel_frac
        self.abs_frac = abs_frac
        self.gene_id = -1
        self.tss_id = -1


class Cluster(object):
    def __init__(self):
        self._id = 0
        self.expr = 1e-10
        self.nodes = set()
        self.paths = []

    def __len__(self):
        return len(self.paths)

    def _add(self, path, expr):
        self.expr += expr
        self.nodes.update(path)
        self.paths.append((path, expr))

    def iterpaths(self):
        if len(self.paths) == 0:
            return
        best_expr = self.paths[0][1]
        for path, expr in self.paths:
            rel_frac = 0.0 if best_expr == 0.0 else expr / best_expr
            abs_frac = 0.0 if self.expr == 0.0 else expr / self.expr
            yield path, expr, rel_frac, abs_frac

    def merge(self, *others):
        for other in others:
            self.expr += other.expr
            self.nodes.update(other.nodes)
            self.paths.extend(other.paths)
        self.paths.sort(key=itemgetter(1), reverse=True)

    @staticmethod
    def build(paths, min_frac=0.0):
        filtered = []
        clusters = {}
        _id = 0
        for i in xrange(len(paths)):
            path, expr = paths[i]
            # does path overlap existing clusters?
            matches = [c for c in clusters.itervalues()
                       if not c.nodes.isdisjoint(path)]
            if len(matches) == 0:
                # make new cluster
                c = Cluster()
                c._id = _id
                _id += 1
                clusters[c._id] = c
            else:
                # check frac in all clusters
                discard = False
                for c in clusters.itervalues():
                    best_expr = c.paths[0][1]
                    rel_frac = 0.0 if best_expr == 0.0 else expr / best_expr
                    if rel_frac < min_frac:
                        discard = True
                        break
                if discard:
                    filtered.append((path, expr))
                    continue
                # merge clusters
                c = matches[0]
                c.merge(*matches[1:])
                for c2 in matches[1:]:
                    del clusters[c2._id]
            # add to cluster
            c._add(path, expr)
        return clusters.values(), filtered


def get_gtf_features(chrom, strand, exons, locus_id, gene_id, tss_id,
                     transcript_id, expr, rel_frac, abs_frac):
    tx_start = exons[0].start
    tx_end = exons[-1].end
    strand_str = Strand.to_gtf(strand)
    attr_dict = {'locus_id': locus_id,
                 'gene_id': gene_id,
                 'tss_id': tss_id,
                 'transcript_id': transcript_id}
    f = GTF.Feature()
    f.seqid = chrom
    f.source = 'taco'
    f.feature = 'transcript'
    f.start = tx_start
    f.end = tx_end
    f.score = int(round(1000.0 * rel_frac))
    f.strand = strand_str
    f.phase = '.'
    f.attrs = {'expr': '%.3f' % expr,
               'rel_frac': '%.5f' % rel_frac,
               'abs_frac': '%.5f' % abs_frac}
    f.attrs.update(attr_dict)
    yield f
    for e in exons:
        f = GTF.Feature()
        f.seqid = chrom
        f.source = 'taco'
        f.feature = 'exon'
        f.start = e.start
        f.end = e.end
        f.score = int(round(1000.0 * rel_frac))
        f.strand = strand_str
        f.phase = '.'
        f.attrs = {}
        f.attrs.update(attr_dict)
        yield f


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
              Strand.to_gtf(strand),
              str(tx_start),
              str(tx_start),
              '0',
              str(len(exons)),
              ','.join(map(str, block_sizes)) + ',',
              ','.join(map(str, block_starts)) + ',']
    return fields


def assign_ids(isoforms, strand, gene_id_iter, tss_id_iter):
    # map tss positions to unique ids
    tss_pos_id_map = {}
    gene_id = gene_id_iter.next()
    for isoform in isoforms:
        start = isoform.exons[0].start
        end = isoform.exons[-1].end
        # map TSS positions to IDs
        tss_pos = end if strand == Strand.NEG else start
        if tss_pos not in tss_pos_id_map:
            tss_id = tss_id_iter.next()
            tss_pos_id_map[tss_pos] = tss_id
        else:
            tss_id = tss_pos_id_map[tss_pos]
        isoform.gene_id = gene_id
        isoform.tss_id = tss_id


def assemble_isoforms(sgraph, config):
    # read in transfrag paths
    pgf = PathGraphFactory(sgraph)
    K, k = pgf.create_optimal(kmax=config.path_graph_kmax,
                              stats_fh=config.path_graph_stats_fh)
    if K is None or len(K) == 0:
        return []
    # smooth kmer graph
    K.apply_smoothing()
    # find isoforms
    logging.debug('%s finding isoforms in k=%d graph (%d kmers) '
                  'source_expr=%f' %
                  (sgraph, k, len(K), K.exprs[K.SOURCE_ID]))
    paths = []
    for path_kmers, expr in find_paths(K, config.path_frac, config.max_paths):
        # convert path of kmers back to path of nodes in splice graph
        path = K.reconstruct(path_kmers)
        paths.append((path, expr))
    logging.debug('%s isoforms: %d' % (sgraph, len(paths)))

    # build gene clusters
    clusters, filtered = Cluster.build(paths, min_frac=config.isoform_frac)
    logging.debug('%s gene clusters: %d filtered transfrags: %d' %
                  (sgraph, len(clusters), len(filtered)))
    gene_isoforms = []
    for cluster in clusters:
        isoforms = []
        for path, expr, rel_frac, abs_frac in cluster.iterpaths():
            # convert from splice graph node ids to exons
            exons = sgraph.reconstruct_exons(path)
            isoforms.append(Isoform(exons=exons,
                                    expr=expr,
                                    rel_frac=rel_frac,
                                    abs_frac=abs_frac))
        # apply max isoforms limit (per cluster)
        if config.max_isoforms > 0:
            isoforms = isoforms[:config.max_isoforms]
        gene_isoforms.append(isoforms)
    return gene_isoforms


def assemble_gene(sgraph, locus_id_str, config):
    logging.debug('%s locus: %s nodes: %d' %
                  (sgraph, locus_id_str, len(sgraph.G)))
    if config.change_point:
        # detect change points
        changepts = sgraph.detect_change_points(
            pval=config.change_point_pvalue,
            fc_cutoff=config.change_point_fold_change)
        logging.debug('%s locus %s change points: %d' %
                      (sgraph, locus_id_str, len(changepts)))
        for cp in changepts:
            sgraph.apply_change_point(cp, config.change_point_trim)
            # output splice graph change points
            for f in sgraph.get_change_point_gtf(cp):
                print >>config.change_point_gtf_fh, str(f)
        # must recreate splice graph after finding change points
        if len(changepts) > 0:
            sgraph.recreate()

    # output splice graph node data
    for f in sgraph.get_node_gtf():
        print >>config.splice_graph_gtf_fh, str(f)

    # run isoform path finding algorithm, filter and group into genes
    for gene_isoforms in assemble_isoforms(sgraph, config):
        # assign gene_id and tss_id
        assign_ids(gene_isoforms, sgraph.strand, config.gene_id_iter,
                   config.tss_id_iter)
        # write output
        for isoform in gene_isoforms:
            # assign transcript id
            t_id = config.t_id_iter.next()
            # get strings for each id
            t_id_str = "TU%d" % t_id
            tss_id_str = "TSS%d" % (isoform.tss_id)
            gene_id_str = "G%d" % (isoform.gene_id)
            # normalize expression by number of samples
            expr_norm = isoform.expr / config.num_samples
            # write to GTF
            for f in get_gtf_features(chrom=sgraph.chrom,
                                      strand=sgraph.strand,
                                      exons=isoform.exons,
                                      locus_id=locus_id_str,
                                      gene_id=gene_id_str,
                                      tss_id=tss_id_str,
                                      transcript_id=t_id_str,
                                      expr=expr_norm,
                                      rel_frac=isoform.rel_frac,
                                      abs_frac=isoform.abs_frac):
                print >>config.assembly_gtf_fh, str(f)
            # write to BED
            name = "%s|%s(%.1f)" % (gene_id_str, t_id_str, expr_norm)
            fields = write_bed(sgraph.chrom, name, sgraph.strand,
                               int(round(1000.0 * isoform.rel_frac)),
                               isoform.exons)
            print >>config.assembly_bed_fh, '\t'.join(fields)


def assemble_locus(locus_index, transfrags, config):
    if len(transfrags) == 0:
        return
    # create locus
    locus = Locus.create(transfrags,
                         config.guided_strand,
                         config.guided_ends,
                         config.guided_assembly)
    logging.debug('%s locus: %s' % (locus, locus_index.name))
    # resolve unstranded transcripts
    num_resolved = locus.impute_unknown_strands()
    if num_resolved > 0:
        logging.debug('%s locus: %s resolved: %d' %
                      (locus, locus_index.name, num_resolved))
    # write bedgraph files after strand resolved
    locus.write_bedgraph(config.bedgraph_fhs)
    # write splice junctions
    locus.write_splice_bed(config.splice_bed_fh)
    # convert to stranded locus objects
    for sgraph in locus.create_splice_graphs():
        if (sgraph.strand == Strand.NA) and not config.assemble_unstranded:
            continue
        assemble_gene(sgraph, locus_index.name, config)


def parse_locus(locus, fh):
    genome_id_str = '%s:%d-%d' % (locus.chrom, locus.start, locus.end)
    logging.debug('LocusIndex: %s coords: %s transfrags: %d' %
                  (locus.name, genome_id_str, locus.num_lines))
    # fast-forward to 'filepos'
    fh.seek(locus.filepos)
    # parse 'num_lines' from file into Transfrag objects
    transfrags = []
    for i in xrange(locus.num_lines):
        transfrags.append(Transfrag.from_bed(fh.next()))
    return transfrags


def parse_locus_index(filename):
    with open(filename) as fileh:
        for line in fileh:
            fields = line.strip().split('\t')
            name = fields[0]
            chrom = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            filepos = int(fields[4])
            num_lines = int(fields[5])
            yield LocusIndex(name, chrom, start, end, filepos, num_lines)


class LockValue(object):
    def __init__(self, initval=0):
        self.val = Value('L', initval)
        self.lock = Lock()

    def next(self):
        with self.lock:
            cur_val = self.val.value
            self.val.value += 1
            return cur_val


class GlobalIds(object):
    def __init__(self):
        # shared memory values (for parallelism)
        self.gene_id_iter = LockValue(1)
        self.tss_id_iter = LockValue(1)
        self.t_id_iter = LockValue(1)


class WorkerState(object):
    def __init__(self, bed_file, input_queue, global_ids,
                 runtime_args, num_samples, output_dir):
        self.bed_file = bed_file
        self.input_queue = input_queue
        self.num_samples = num_samples
        for k, v in vars(global_ids).iteritems():
            setattr(self, k, v)
        for k, v in vars(runtime_args).iteritems():
            setattr(self, k, v)
        # set locations for worker results
        r = Results(output_dir)
        self.results = r
        # open files
        self.bedgraph_fhs = []
        for filename in r.bedgraph_files:
            self.bedgraph_fhs.append(open(filename, 'w'))
        self.splice_bed_fh = open(r.splice_bed_file, 'w')
        self.splice_graph_gtf_fh = open(r.splice_graph_gtf_file, 'w')
        self.change_point_gtf_fh = open(r.change_point_gtf_file, 'w')
        self.path_graph_stats_fh = open(r.path_graph_stats_file, 'w')
        self.assembly_gtf_fh = open(r.assembly_gtf_file, 'w')
        self.assembly_bed_fh = open(r.assembly_bed_file, 'w')

    def close(self):
        # close files
        for fh in self.bedgraph_fhs:
            fh.close()
        self.splice_bed_fh.close()
        self.splice_graph_gtf_fh.close()
        self.change_point_gtf_fh.close()
        self.path_graph_stats_fh.close()
        self.assembly_gtf_fh.close()
        self.assembly_bed_fh.close()

    def sort_output_files(self):
        # create output directories
        results = self.results
        sort_output_dir = os.path.join(results.output_dir, 'sort')
        sort_tmp_dir = os.path.join(results.output_dir, 'tmp')
        if not os.path.exists(sort_output_dir):
            logging.debug('\tcreating sort dir %s' % (sort_output_dir))
            os.makedirs(sort_output_dir)
        if not os.path.exists(sort_tmp_dir):
            logging.debug('\tcreating sort tmp dir %s' % (sort_tmp_dir))
            os.makedirs(sort_tmp_dir)

        # create new set of sorted results
        sorted_results = Results(sort_output_dir)

        # bedgraph files
        logging.debug('\t%s bedgraph files' % (results.output_dir))
        for filename, sorted_filename in zip(results.bedgraph_files,
                                             sorted_results.bedgraph_files):
            batch_sort(input=filename,
                       output=sorted_filename,
                       key=sort_key_bed,
                       buffer_size=SORT_BUFFER_SIZE,
                       tempdirs=[sort_tmp_dir])
            os.rename(sorted_filename, filename)
        # splice bed file
        logging.debug('\t%s splice bed file' % (results.output_dir))
        batch_sort(input=results.splice_bed_file,
                   output=sorted_results.splice_bed_file,
                   key=sort_key_bed,
                   buffer_size=SORT_BUFFER_SIZE,
                   tempdirs=[sort_tmp_dir])
        os.rename(sorted_results.splice_bed_file, results.splice_bed_file)
        # splice graph gtf
        logging.debug('\t%s splice graph gtf file' % (results.output_dir))
        batch_sort(input=results.splice_graph_gtf_file,
                   output=sorted_results.splice_graph_gtf_file,
                   key=sort_key_gtf,
                   buffer_size=SORT_BUFFER_SIZE,
                   tempdirs=[sort_tmp_dir])
        os.rename(sorted_results.splice_graph_gtf_file,
                  results.splice_graph_gtf_file)
        # change point gtf
        logging.debug('\t%s change point gtf file' % (results.output_dir))
        batch_sort(input=results.change_point_gtf_file,
                   output=sorted_results.change_point_gtf_file,
                   key=sort_key_gtf,
                   buffer_size=SORT_BUFFER_SIZE,
                   tempdirs=[sort_tmp_dir])
        os.rename(sorted_results.change_point_gtf_file,
                  results.change_point_gtf_file)
        # path graph stats
        logging.debug('\t%s path graph stats file' % (results.output_dir))
        batch_sort(input=results.path_graph_stats_file,
                   output=sorted_results.path_graph_stats_file,
                   key=sort_key_bed,
                   buffer_size=SORT_BUFFER_SIZE,
                   tempdirs=[sort_tmp_dir])
        os.rename(sorted_results.path_graph_stats_file,
                  results.path_graph_stats_file)
        # assembly bed
        logging.debug('\t%s assembly bed file' % (results.output_dir))
        batch_sort(input=results.assembly_bed_file,
                   output=sorted_results.assembly_bed_file,
                   key=sort_key_bed,
                   buffer_size=SORT_BUFFER_SIZE,
                   tempdirs=[sort_tmp_dir])
        os.rename(sorted_results.assembly_bed_file,
                  results.assembly_bed_file)
        # assembly gtf
        logging.debug('\t%s assembly gtf file' % (results.output_dir))
        batch_sort(input=results.assembly_gtf_file,
                   output=sorted_results.assembly_gtf_file,
                   key=sort_key_gtf,
                   buffer_size=SORT_BUFFER_SIZE,
                   tempdirs=[sort_tmp_dir])
        os.rename(sorted_results.assembly_gtf_file,
                  results.assembly_gtf_file)
        # remove temporary directories
        logging.debug('\t%s cleaning up' % (results.output_dir))
        os.rmdir(sort_output_dir)
        os.rmdir(sort_tmp_dir)


def assemble_worker(state):
    # process loci via input queue
    bed_fh = open(state.bed_file)
    while True:
        locus = state.input_queue.get()
        if locus is None:
            break
        assemble_locus(locus, parse_locus(locus, bed_fh), state)
        state.input_queue.task_done()
    state.input_queue.task_done()
    # cleanup and close files
    state.close()
    # sort output files
    logging.debug('Sorting worker output files: "%s"' %
                  (state.results.output_dir))
    state.sort_output_files()


def assemble_parallel(args, results, num_samples):
    '''
    args: object containing parameters to configure the assembly process
    results: Results object containing input and output filenames
    num_samples: number of samples in assembly

    Args
    ====
    - guided_strand
    - guided_ends
    - guided_assembly
    - change_point
    - change_point_pvalue
    - change_point_fold_change
    - change_point_trim
    - path_graph_kmax
    - path_frac
    - max_paths
    - isoform_frac
    - max_isoforms
    - assemble_unstranded

    Results
    =======
    Input file attributes:
    - locus_index_file
    - transfrags_bed_file

    Output file attributes:
    - bedgraph_files
    - splice_bed_file
    - splice_graph_gtf_file
    - change_point_gtf_file
    - path_graph_stats_file
    - assembly_gtf_file
    - assembly_bed_file
    '''
    logging.info('Assembling in parallel using %d processes' %
                 (args.num_processes))
    # create queue
    input_queue = JoinableQueue(maxsize=args.num_processes * 2)
    bed_file = results.transfrags_bed_file
    global_ids = GlobalIds()
    # start worker processes
    procs = []
    worker_results = []
    for i in xrange(args.num_processes):
        worker_id = 'worker%03d' % i
        worker_dir = os.path.join(results.tmp_dir, worker_id)
        if not os.path.exists(worker_dir):
            logging.debug("\tcreating worker directory '%s'" % (worker_dir))
            os.makedirs(worker_dir)
        worker_results.append(Results(worker_dir))
        worker_state = WorkerState(bed_file, input_queue, global_ids,
                                   args, num_samples, worker_dir)
        p = Process(target=assemble_worker, args=(worker_state,))
        p.start()
        procs.append(p)
    # parse locus file
    for locus in parse_locus_index(results.locus_index_file):
        input_queue.put(locus)
    for p in procs:
        input_queue.put(None)
    # close input queue
    input_queue.join()
    input_queue.close()
    # join worker processes
    for p in procs:
        p.join()

    # merge output files
    def merge(input_files, output_file, key, header=None):
        fhs = [open(f, 'rb', 64*1024) for f in input_files]
        with open(output_file, 'wb', 64*1024) as output:
            if header is not None:
                output.write(header)
            iterator = batch_merge(key, *fhs)
            output.writelines(iterator)
        for fh in fhs:
            fh.close()

    logging.info('Merging output files')
    logging.debug('\tmerging bedgraph files')
    for i, output_file in enumerate(results.bedgraph_files):
        input_files = [r.bedgraph_files[i] for r in worker_results]
        merge(input_files, output_file, sort_key_bed)
    logging.debug('\tmerging splice bed file')
    header = ('track name=junctions description="Splice Junctions" '
              'graphType=junctions\n')
    merge(input_files=[r.splice_bed_file for r in worker_results],
          output_file=results.splice_bed_file,
          key=sort_key_bed,
          header=header)
    logging.debug('\tmerging splice graph gtf file')
    merge(input_files=[r.splice_graph_gtf_file for r in worker_results],
          output_file=results.splice_graph_gtf_file,
          key=sort_key_gtf)
    logging.debug('\tmerging change point gtf file')
    merge(input_files=[r.change_point_gtf_file for r in worker_results],
          output_file=results.change_point_gtf_file,
          key=sort_key_gtf)
    logging.debug('\tmerging path graph stats file')
    header = ['chrom', 'start', 'end', 'strand', 'k', 'kmax', 'transfrags',
              'short_transfrags', 'short_expr', 'lost_short', 'lost_short_expr',
              'kmers', 'lost_kmers', 'tot_expr', 'graph_expr', 'expr_frac',
              'valid', 'opt', 'is_opt\n']
    header = '\t'.join(header)
    merge(input_files=[r.path_graph_stats_file for r in worker_results],
          output_file=results.path_graph_stats_file,
          key=sort_key_bed,
          header=header)
    logging.debug('\tmerging assembly bed file')
    merge(input_files=[r.assembly_bed_file for r in worker_results],
          output_file=results.assembly_bed_file,
          key=sort_key_bed)
    logging.debug('\tmerging assembly gtf file')
    merge(input_files=[r.assembly_gtf_file for r in worker_results],
          output_file=results.assembly_gtf_file,
          key=sort_key_gtf)
    # cleanup worker data
    logging.info('Removing temporary files')
    def shutil_error_callback(func, path, excinfo):
        logging.error('Error removing tmp files path=%s message=%s' %
                      (path, excinfo))
    for r in worker_results:
        shutil.rmtree(r.output_dir, onerror=shutil_error_callback)
    logging.info('Done')
    return 0
