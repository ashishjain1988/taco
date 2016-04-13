###############################################################################
# Adapted from Pysam (cfaidx.pyx)
###############################################################################
#
# The MIT License
#
# Copyright (c) 2015 Andreas Heger
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
###############################################################################
import sys
import os
import re
from cpython cimport array
from cpython cimport PyBytes_Check, PyUnicode_Check

from libc.stdlib cimport free

from chtslib cimport \
    faidx_t, \
    faidx_nseq, fai_load, fai_destroy, fai_fetch, \
    faidx_seq_len, \
    faidx_fetch_seq

# hard-coded constants
cdef int MAX_POS = 2 << 29

cdef charptr_to_str(const char* s):
    if s == NULL:
        return None
    return s

cdef force_str(object s):
    if s is None:
        return None
    return s

cdef bytes force_bytes(object s):
    return s

# filename encoding (copied from lxml.etree.pyx)
cdef str _FILENAME_ENCODING
_FILENAME_ENCODING = sys.getfilesystemencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = sys.getdefaultencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = 'ascii'

cdef bytes encode_filename(object filename):
    """Make sure a filename is 8-bit encoded (or None)."""
    if filename is None:
        return None
    elif PyBytes_Check(filename):
        return filename
    elif PyUnicode_Check(filename):
        return filename.encode(_FILENAME_ENCODING)
    else:
        raise TypeError(u"Argument must be string or unicode.")


cpdef parse_region(reference=None,
                   start=None,
                   end=None,
                   region=None):
    cdef int rtid
    cdef long long rstart
    cdef long long rend

    rtid = -1
    rstart = 0
    rend = MAX_POS
    if start != None:
        try:
            rstart = start
        except OverflowError:
            raise ValueError('start out of range (%i)' % start)

    if end != None:
        try:
            rend = end
        except OverflowError:
            raise ValueError('end out of range (%i)' % end)

    if region:
        region = force_str(region)
        parts = re.split("[:-]", region)
        reference = parts[0]
        if len(parts) >= 2:
            rstart = int(parts[1]) - 1
        if len(parts) >= 3:
            rend = int(parts[2])

    if not reference:
        return None, 0, 0

    if not 0 <= rstart < MAX_POS:
        raise ValueError('start out of range (%i)' % rstart)
    if not 0 <= rend <= MAX_POS:
        raise ValueError('end out of range (%i)' % rend)
    if rstart > rend:
        raise ValueError(
            'invalid region: start (%i) > end (%i)' % (rstart, rend))

    return force_bytes(reference), rstart, rend


cdef class FastaFile:
    cdef object _filename, _references, _lengths, reference2length
    cdef faidx_t* fastafile
    cdef char* _fetch(self, char* reference,
                      int start, int end, int* length)

    def __cinit__(self, *args, **kwargs):
        self.fastafile = NULL
        self._filename = None
        self._references = None
        self._lengths = None
        self.reference2length = None
        self._open(*args, **kwargs)

    def is_open(self):
        '''return true if samfile has been opened.'''
        return self.fastafile != NULL

    def __len__(self):
        if self.fastafile == NULL:
            raise ValueError("calling len() on closed file")

        return faidx_nseq(self.fastafile)

    def _open(self, filename, filepath_index=None):
        '''open an indexed fasta file.

        This method expects an indexed fasta file.
        '''

        # close a previously opened file
        if self.fastafile != NULL:
            self.close()

        self._filename = encode_filename(filename)
        cdef char *cfilename = self._filename

        # open file for reading
        if (self._filename != b"-"
            and not os.path.exists(filename)):
            raise IOError("file `%s` not found" % filename)

        with nogil:
            self.fastafile = fai_load(cfilename)

        if self.fastafile == NULL:
            raise IOError("could not open file `%s`" % filename)

        if filepath_index is None:
            filepath_index = filename + ".fai"

        if not os.path.exists(filepath_index):
            raise ValueError("could not locate index file {}".format(
                filepath_index))

        with open(filepath_index) as inf:
            data = [x.split("\t") for x in inf]
            self._references = tuple(x[0] for x in data)
            self._lengths = tuple(int(x[1]) for x in data)
            self.reference2length = dict(zip(self._references, self._lengths))

    def close(self):
        """close the file."""
        if self.fastafile != NULL:
            fai_destroy(self.fastafile)
            self.fastafile = NULL

    def __dealloc__(self):
        self.close()

    # context manager interface
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    property closed:
        """"bool indicating the current state of the file object.
        This is a read-only attribute; the close() method changes the value.
        """
        def __get__(self):
            return not self.is_open()

    property filename:
        """filename associated with this object. This is a read-only attribute."""
        def __get__(self):
            return self._filename

    property references:
        '''tuple with the names of :term:`reference` sequences.'''
        def __get__(self):
            return self._references

    property nreferences:
        """"int with the number of :term:`reference` sequences in the file.
        This is a read-only attribute."""
        def __get__(self):
            return len(self._references) if self.references else None

    property lengths:
        """tuple with the lengths of :term:`reference` sequences."""
        def __get__(self):
            return self._lengths

    def fetch(self,
              reference=None,
              start=None,
              end=None,
              region=None):
        """fetch sequences in a :term:`region`.

        A region can
        either be specified by :term:`reference`, `start` and
        `end`. `start` and `end` denote 0-based, half-open
        intervals.

        Alternatively, a samtools :term:`region` string can be
        supplied.

        If any of the coordinates are missing they will be replaced by the
        minimum (`start`) or maximum (`end`) coordinate.

        Note that region strings are 1-based, while `start` and `end` denote
        an interval in python coordinates.
        The region is specified by :term:`reference`, `start` and `end`.

        Returns
        -------

        string : a string with the sequence specified by the region.

        Raises
        ------

        IndexError
            if the coordinates are out of range

        ValueError
            if the region is invalid

        """

        if not self.is_open():
            raise ValueError("I/O operation on closed file" )

        cdef int length
        cdef char *seq
        cdef char *ref
        cdef int rstart, rend

        reference, rstart, rend = parse_region(reference, start, end, region)

        if reference is None:
            raise ValueError("no sequence/region supplied.")

        if rstart == rend:
            return ""

        ref = reference
        with nogil:
            length = faidx_seq_len(self.fastafile, ref)
        if length == -1:
            raise KeyError("sequence '%s' not present" % reference)
        if rstart >= length:
            return ""

        # fai_fetch adds a '\0' at the end
        with nogil:
            seq = faidx_fetch_seq(self.fastafile,
                                  ref,
                                  rstart,
                                  rend-1,
                                  &length)

        if seq == NULL:
            raise ValueError(
                "failure when retrieving sequence on '%s'" % reference)

        try:
            return charptr_to_str(seq)
        finally:
            free(seq)

    cdef char * _fetch(self, char * reference, int start, int end, int * length):
        '''fetch sequence for reference, start and end'''

        with nogil:
            return faidx_fetch_seq(self.fastafile,
                                   reference,
                                   start,
                                   end-1,
                                   length)

    def get_reference_length(self, reference):
        '''return the length of reference.'''
        return self.reference2length[reference]

    def __getitem__(self, reference):
        return self.fetch(reference)

    def __contains__(self, reference):
        '''return true if reference in fasta file.'''
        return reference in self.reference2length


__all__ = ["FastaFile"]
