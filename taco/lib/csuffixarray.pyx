'''
TACO: Transcriptome meta-assembly for RNA-Seq
'''
from cpython cimport array
import array
from bisect import bisect_right

# local imports
from sais cimport sais_int, sais_int_bwt

DEF TERMINATOR = 0
DEF SPACER = 1

cdef int suffix_cmp(tuple p, int* t, int* sa, int n, int c):
    cdef int i = 0
    cdef int np = len(p)
    cdef int retval = 0

    while i < len(p) and (sa[c]+i) < n:
        if p[i] > t[sa[c]+i]:
            retval = 1
            break
        elif p[i] < t[sa[c]+i]:
            retval = -1
            break
        i += 1
    return retval


cdef tuple suffix_range(tuple p, int* t, int* sa, int n):
    cdef int l, r, s, c, retval

    l = 0
    r = n
    while l < r:
        c = (l+r) // 2
        retval = suffix_cmp(p, t, sa, n, c)
        if retval > 0:
            l = c + 1
        else:
            r = c
    s = l
    r = n
    while l < r:
        c = (l+r) // 2
        retval = suffix_cmp(p, t, sa, n, c)
        if retval < 0:
            r = c
        else:
            l = c + 1
    return (s, r)


cdef suffix_array(int* t, int* sa, int n, int k):
    sais_int(t, sa, n, k)


cdef class SuffixArrayIndex(object):

    cdef array.array t_arr
    cdef array.array s_arr
    cdef int* t
    cdef int* sa
    cdef int n
    cdef list starts
    cdef list lengths
    cdef int alphabet_size

    def __init__(self, seqs):
        cdef object seq
        cdef list t, starts, lengths
        cdef int i, alphabet_size

        t = []
        starts = []
        lengths = []
        i = 0
        alphabet_size = 0
        for seq in seqs:
            t.extend(seq)
            t.append(SPACER)
            starts.append(i)
            lengths.append(len(seq))
            i += len(seq) + 1
        t.append(TERMINATOR)

        self.alphabet_size = max(t) + 1
        self.starts = starts
        self.lengths = lengths
        self.t_arr = array.array('i', t)
        self.s_arr = array.array('i', t)
        self.t = self.t_arr.data.as_ints
        self.sa = self.s_arr.data.as_ints
        self.n = len(t)
        # create suffix array
        sais_int(self.t, self.sa, self.n, self.alphabet_size)

    def search(self, tuple p):
        cdef int i, l, r, si, istart, start, end

        l, r = suffix_range(p, self.t, self.sa, self.n)
        for i in xrange(l, r):
            si = self.sa[i]
            istart = bisect_right(self.starts, si) - 1
            start = self.starts[istart]
            end = start + self.lengths[istart]
            yield tuple(self.t_arr[start:end])
