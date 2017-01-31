'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
'''
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy
numpy_inc = numpy.get_include()

__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2016"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "GPL"
__version__ = "0.5.0"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


cython_extensions = [
    Extension('taco.lib.cbedgraph',
              sources=['taco/lib/cbedgraph.pyx'],
              include_dirs=[numpy_inc],
              extra_compile_args=['-w']),
    Extension('taco.lib.cchangepoint',
              sources=['taco/lib/cchangepoint.pyx'],
              include_dirs=[numpy_inc],
              extra_compile_args=['-w']),
    Extension('taco.lib.bx.cluster',
              sources=['taco/lib/bx/cluster.pyx',
                       'taco/lib/bx/intervalcluster.c'],
              include_dirs=['taco/lib/bx'],
              extra_compile_args=['-w']),
    Extension('taco.lib.bx.intersection',
              sources=['taco/lib/bx/intersection.pyx'],
              extra_compile_args=['-w']),
    Extension('taco.lib.csuffixarray',
              sources=['taco/lib/csuffixarray.pyx', 'taco/lib/sais.c'],
              extra_compile_args=['-w']),
    Extension('taco.lib.cpathfinder',
              sources=['taco/lib/cpathfinder.pyx'],
              extra_compile_args=['-w']),
    Extension('taco.lib.cbisect',
              sources=['taco/lib/cbisect.pyx', 'taco/lib/bsearch.c'],
              extra_compile_args=['-w']),
    Extension('taco.lib.scipy.norm_sf',
              sources=['taco/lib/scipy/norm_sf.pyx',
                       'taco/lib/scipy/ndtr.c',
                       'taco/lib/scipy/const.c',
                       'taco/lib/scipy/mtherr.c',
                       'taco/lib/scipy/sf_error.c'],
              include_dirs=[numpy_inc, 'taco/lib/scipy'],
              extra_compile_args=['-w']),
    Extension('taco.lib.pysam.cfaidx',
              sources=['taco/lib/pysam/cfaidx.pyx',
                       'taco/lib/htslib/faidx.c',
                       'taco/lib/htslib/bgzf.c',
                       'taco/lib/htslib/hfile.c'],
              include_dirs=['taco/lib/pysam',
                            'taco/lib/htslib'],
              libraries=['z'])
]

extensions = [
    Extension('taco.lib.clocusindex',
              sources=['taco/lib/clocusindex.c'],
              include_dirs=['taco/lib'])
]


def main():
    datadir = os.path.join('share', 'data')
    datafiles = [(d, [os.path.join(d,f) for f in files])
        for d, folders, files in os.walk(datadir)]
    setup(name='taco',
          version=__version__,
          description='transcriptome meta-assembly for rna-seq',
          author=__author__,
          author_email=__email__,
          requires=['numpy', 'cython'],
          license=__license__,
          platforms='Linux',
          url='https://github.com/tacorna/taco',
          ext_modules=extensions + cythonize(cython_extensions),
          packages=['taco', 'taco.lib', 'taco.lib.bx', 'taco.lib.scipy',
                    'taco.lib.pysam'],
          scripts=['taco/taco_run.py', 'taco/taco_refcomp.py'],
          data_files = datafiles)


if __name__ == '__main__':
    main()
