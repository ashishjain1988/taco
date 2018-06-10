'''
TACO: Transcriptome meta-assembly for RNA-Seq
'''
from cpython cimport array
import array


__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2018"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "MIT"
__version__ = "0.7.3"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


# constant minimum path score
DEF MIN_SCORE = 1.0e-10


cdef find_path(int* nodes, int nnodes,
               float *exprs, list succs,
               int n, int source, int sink):
    cdef array.array int_array_template = array.array('i')
    cdef array.array float_array_template = array.array('f')
    cdef array.array min_exprs_arr, prevs_arr
    cdef float *min_exprs
    cdef int *prevs
    cdef float min_expr, new_min_expr
    cdef int x, i, j, prev
    cdef list path

    # allocate data structures
    min_exprs_arr = array.clone(float_array_template, n, zero=False)
    prevs_arr = array.clone(int_array_template, n, zero=False)

    # initialize data structures
    min_exprs = min_exprs_arr.data.as_floats
    prevs = prevs_arr.data.as_ints
    for i in xrange(n):
        min_exprs[i] = MIN_SCORE
        prevs[i] = sink
    min_exprs[source] = exprs[source]

    # traverse nodes in topological sort order
    for x in xrange(nnodes):
        i = nodes[x]
        min_expr = min_exprs[i]

        for j in succs[i]:
            new_min_expr = min_expr if min_expr < exprs[j] else exprs[j]
            if ((prevs[j] == sink) or (new_min_expr > min_exprs[j])):
                min_exprs[j] = new_min_expr
                prevs[j] = i

    # traceback to get path
    expr = min_exprs[sink]
    prev = sink
    path = [sink]
    while True:
        prev = prevs[prev]
        path.append(prev)
        if prev == source:
            break
    path.reverse()

    # subtract path
    for i in xrange(len(path)):
        j = path[i]
        new_expr = exprs[j] - expr
        exprs[j] = MIN_SCORE if MIN_SCORE >= new_expr else new_expr
    return path, expr


def find_paths(object G, float path_frac=0, int max_paths=0):
    cdef array.array nodes_arr, exprs_arr
    cdef int *nodes
    cdef float *exprs
    cdef int nnodes, n
    cdef int source, sink, iterations
    cdef float expr, lowest_expr
    cdef list path
    cdef list results

    # don't run if all nodes are zero
    if G.exprs[G.SOURCE_ID] < MIN_SCORE:
        return []

    # initialize data structures
    nodes_arr = array.array('i', G.topological_sort())
    exprs_arr = array.array('f', G.exprs)
    nodes = nodes_arr.data.as_ints
    exprs = exprs_arr.data.as_floats
    nnodes = len(nodes_arr)
    n = len(exprs_arr)
    source = G.SOURCE_ID
    sink = G.SINK_ID

    # find highest scoring path
    path, expr = find_path(nodes, nnodes, exprs, G.succs, n, source, sink)
    results = [(path, expr)]

    # define threshold score to stop producing paths
    lowest_expr = expr * path_frac
    if MIN_SCORE > lowest_expr:
        lowest_expr = MIN_SCORE

    # iterate to find paths
    iterations = 1
    while True:
        if max_paths > 0 and iterations >= max_paths:
            break
        # find path
        path, expr = find_path(nodes, nnodes, exprs, G.succs, n, source, sink)
        if expr <= lowest_expr:
            break
        # store path
        results.append((path, expr))
        iterations += 1
    return results
