'''
Created on Dec 1, 2013

@author: mkiyer
'''
# based on Recipe 466302: Sorting big files the Python 2.4 way
# by Nicolas Lehuen
# http://code.activestate.com/recipes/576755-sorting-big-files-the-python-26-way/

import os
from tempfile import gettempdir
from itertools import islice, cycle
from collections import namedtuple
import heapq

Keyed = namedtuple("Keyed", ["key", "obj"])

# batch sort configuration
SORT_BUFFER_SIZE = 32000


def sort_key_bed(line):
    fields = line.split('\t', 2)
    return (fields[0], int(fields[1]))


def sort_key_gtf(line):
    fields = line.split('\t', 4)
    feature_key = 0 if fields[2] == 'transcript' else 1
    return (fields[0], int(fields[3]), feature_key)


def merge_files(input_files, output_file, key, header=None):
    fhs = [open(f, 'rb', 64*1024) for f in input_files]
    with open(output_file, 'wb', 64*1024) as output:
        if header is not None:
            output.write(header)
        iterator = batch_merge(key, *fhs)
        output.writelines(iterator)
    for fh in fhs:
        fh.close()


def batch_merge(key=None, *iterables):
    # based on code posted by Scott David Daniels in c.l.p.
    # http://groups.google.com/group/comp.lang.python/msg/484f01f1ea3c832d

    if key is None:
        keyed_iterables = iterables
    else:
        keyed_iterables = [(Keyed(key(obj), obj) for obj in iterable)
                            for iterable in iterables]
    for element in heapq.merge(*keyed_iterables):
        yield element.obj


def batch_sort(input, output, key=None, buffer_size=32000, tempdirs=None):
    if tempdirs is None:
        tempdirs = []
    if not tempdirs:
        tempdirs.append(gettempdir())

    chunks = []
    try:
        with open(input,'rb',64*1024) as input_file:
            input_iterator = iter(input_file)
            for tempdir in cycle(tempdirs):
                current_chunk = list(islice(input_iterator,buffer_size))
                if not current_chunk:
                    break
                current_chunk.sort(key=key)
                output_chunk = open(os.path.join(tempdir,'%06i'%len(chunks)),'w+b',64*1024)
                chunks.append(output_chunk)
                output_chunk.writelines(current_chunk)
                output_chunk.flush()
                output_chunk.seek(0)
        with open(output,'wb',64*1024) as output_file:
            output_file.writelines(batch_merge(key, *chunks))
    finally:
        for chunk in chunks:
            try:
                chunk.close()
                os.remove(chunk.name)
            except Exception:
                pass

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser()
    parser.add_option(
        '-b','--buffer',
        dest='buffer_size',
        type='int',default=32000,
        help='''Size of the line buffer. The file to sort is
            divided into chunks of that many lines. Default : 32,000 lines.'''
    )
    parser.add_option(
        '-k','--key',
        dest='key',
        help='''Python expression used to compute the key for each
            line, "lambda line:" is prepended.\n
            Example : -k "line[5:10]". By default, the whole line is the key.'''
    )
    parser.add_option(
        '-t','--tempdir',
        dest='tempdirs',
        action='append',
        default=[],
        help='''Temporary directory to use. You might get performance
            improvements if the temporary directory is not on the same physical
            disk than the input and output directories. You can even try
            providing multiples directories on differents physical disks.
            Use multiple -t options to do that.'''
    )
    options,args = parser.parse_args()

    if options.key:
        options.key = eval('lambda line : (%s)'%options.key)

    batch_sort(args[0],args[1],options.key,options.buffer_size,options.tempdirs)
