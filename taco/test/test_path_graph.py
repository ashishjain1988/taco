'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
'''
import networkx as nx

from taco.lib.base import Exon
from taco.lib.splice_graph import SpliceGraph
from taco.lib.path_graph import PathGraph, PathGraphFactory

from taco.test.base import read_single_locus


def test_path_graph_factory():
    t_dict, locus = read_single_locus('path1.gtf')
    sgraph = SpliceGraph.create(t_dict.values())
    pgraphfactory = PathGraphFactory(sgraph)
    g1 = pgraphfactory.create(k=1)
    assert len(g1) == 5
    g2 = pgraphfactory.create(k=2)
    assert len(g2) == 6
    gopt, k = pgraphfactory.create_optimal()
    assert k == 2
    return


def test_reachability():
    return
    source = -1
    sink = -2
    G = nx.DiGraph()
    nx.path_graph(10, create_using=G)
    G.add_edge(source, 0)
    G.add_edge(9, sink)
    assert len(get_unreachable_kmers(G, source, sink)) == 0
    # unreachable from sink
    G.add_edges_from([(2, 20), (20, 21)])
    unreachable = get_unreachable_kmers(G, source, sink)
    assert unreachable == set([20, 21])
    # unreachable from source
    G.add_edges_from([(33, 32), (32, 31), (31, 9)])
    unreachable = get_unreachable_kmers(G, source, sink)
    assert unreachable == set([20, 21, 33, 32, 31])


def test_unreachable_kmers():
    return
    t_dict, locus = read_single_locus('path_graph_k2.gtf')
    sgraph = SpliceGraph.create(t_dict.values())
    K = create_path_graph(sgraph, k=2)
    assert not K.graph['valid']
    assert len(K) == 0

    K = create_path_graph(sgraph, k=1)
    assert K.graph['valid']
    assert K.graph['num_lost_kmers'] == 0
    assert len(K) == 8

    K, k = create_optimal_path_graph(sgraph, kmax=0, loss_threshold=1.0)
    assert k == 1
    assert len(K) == 8


def test_path_graph1():
    return
    # read transcripts
    t_dict, locus = read_single_locus('path1.gtf')
    SG = SpliceGraph.create(t_dict.values())
    # paths
    ABCDE = (SOURCE, Exon(0, 100), Exon(200, 300), Exon(400, 500),
             Exon(600, 700), Exon(800, 900), SINK)
    ACE = (SOURCE, Exon(0, 100), Exon(400, 500), Exon(800, 900), SINK)
    ABCE = (SOURCE, Exon(0, 100), Exon(200, 300), Exon(400, 500),
            Exon(800, 900), SINK)
    ACDE = (SOURCE, Exon(0, 100), Exon(400, 500), Exon(600, 700),
            Exon(800, 900), SINK)
    paths = [ABCDE, ACE, ABCE, ACDE]
    # create path graph k = 2
    k = 2
    G1 = create_path_graph(SG, k)
    G2 = nx.DiGraph()
    for path in paths:
        kmers = list(get_kmers(path, k))
        add_path(G2, kmers, 1.0)
    assert nx.is_isomorphic(G1, G2)


def test_path_graph2():
    return
    t_dict, locus = read_single_locus('change_point2.gtf')
    sgraph = SpliceGraph.create(t_dict.values())

    # trivial case without additional stops or starts
    k = 1
    K = create_path_graph(sgraph, k)
    kmer_id_map = K.graph['kmer_id_map']
    n_id = sgraph.get_node_id(Exon(0, 100))
    kmer_id = kmer_id_map[(n_id,)]
    assert K.node[kmer_id]['expr'] == 12.0
    assert K.node[SOURCE]['expr'] == 12.0
    assert K.node[SINK]['expr'] == 12.0

    # add a stop site
    sgraph.stop_sites.add(50)
    sgraph.recreate()
    K = create_path_graph(sgraph, k=2)
    kmer_id_map = K.graph['kmer_id_map']
    n1 = (sgraph.get_node_id((0, 50)), sgraph.get_node_id((50, 100)))
    kmer1 = kmer_id_map[n1]
    n2 = (sgraph.get_node_id((0, 50)),)
    kmer2 = kmer_id_map[n2]

    assert K.node[kmer1]['expr'] == 1.0
    assert K.node[kmer2]['expr'] == 10.0
    assert K.node[SOURCE]['expr'] == 11.0
    assert K.node[SINK]['expr'] == 11.0
    # smooth kmer graph
    smooth_graph(K)
    assert K.node[kmer1]['expr'] == 1.0
    assert K.node[kmer2]['expr'] == 10.0
    assert K.node[SOURCE]['expr'] == 11.0
    assert K.node[SINK]['expr'] == 11.0

    # TODO: test after rescuing short transfrags

    # add both a start and a stop site
    sgraph.start_sites.add(50)
    sgraph.stop_sites.add(50)
    sgraph.recreate()
    K = create_path_graph(sgraph, k=2)
    smooth_graph(K)
    kmer_id_map = K.graph['kmer_id_map']
    n1 = (sgraph.get_node_id((0, 50)), sgraph.get_node_id((50, 100)))
    n2 = (sgraph.get_node_id((0, 50)),)
    n3 = (sgraph.get_node_id((50, 100)),)
    kmer1 = kmer_id_map[n1]
    kmer2 = kmer_id_map[n2]
    kmer3 = kmer_id_map[n3]
    assert K.node[kmer1]['expr'] == 1.0
    assert K.node[kmer2]['expr'] == 10.0
    assert K.node[kmer3]['expr'] == 1.0
    assert K.node[SOURCE]['expr'] == 12.0
    assert K.node[SINK]['expr'] == 12.0
