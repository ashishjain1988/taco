'''
TACO: Transcriptome meta-assembly from RNA-Seq
'''
from cpython cimport array
import array


__author__ = "Matthew Iyer and Yashar Niknafs"
__copyright__ = "Copyright 2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs"]
__license__ = "GPL"
__version__ = "0.4.0"
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
    cdef array.array min_exprs_arr, sum_exprs_arr, lengths_arr, prevs_arr
    cdef float *min_exprs, *sum_exprs
    cdef int *lengths, *prevs
    cdef float min_expr, sum_expr, new_min_expr, new_sum_expr
    cdef float new_avg_expr, cur_avg_expr
    cdef int x, i, j, length, new_length, prev
    cdef list path

    # allocate data structures
    min_exprs_arr = array.clone(float_array_template, n, zero=False)
    sum_exprs_arr = array.clone(float_array_template, n, zero=False)
    lengths_arr = array.clone(int_array_template, n, zero=False)
    prevs_arr = array.clone(int_array_template, n, zero=False)

    # initialize data structures
    min_exprs = min_exprs_arr.data.as_floats
    sum_exprs = sum_exprs_arr.data.as_floats
    lengths = lengths_arr.data.as_ints
    prevs = prevs_arr.data.as_ints
    for i in xrange(n):
        min_exprs[i] = MIN_SCORE
        sum_exprs[i] = MIN_SCORE
        lengths[i] = 1
        prevs[i] = sink
    min_exprs[source] = exprs[source]
    sum_exprs[source] = exprs[source]

    # traverse nodes in topological sort order
    for x in xrange(nnodes):
        i = nodes[x]
        min_expr = min_exprs[i]
        sum_expr = sum_exprs[i]
        length = lengths[i]

        for j in succs[i]:
            new_min_expr = min_expr if min_expr < exprs[j] else exprs[j]
            new_sum_expr = sum_expr + exprs[j]
            new_length = length + 1
            new_avg_expr = new_sum_expr / new_length
            cur_avg_expr = sum_exprs[j] / lengths[j]

            # update node if 1) not yet visited, 2) can reach this node
            # with a higher min expr, or 3) can reach with an equal
            # min expr but a higher overall average expr
#            if ((prevs[j] == sink) or (new_min_expr > min_exprs[j]) or
#                (new_min_expr == min_exprs[j] and
#                 new_avg_expr > cur_avg_expr)):
            if ((prevs[j] == sink) or (new_min_expr > min_exprs[j])):
                min_exprs[j] = new_min_expr
                sum_exprs[j] = new_sum_expr
                lengths[j] = new_length
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
    return tuple(path), expr


def find_paths(object G, float path_frac=0, int max_paths=0):
    cdef array.array nodes_arr, exprs_arr
    cdef int *nodes
    cdef float *exprs
    cdef int nnodes, n
    cdef int source, sink, iterations
    cdef float expr, lowest_expr
    cdef tuple path
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
