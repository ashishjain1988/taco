'''
TACO: Multi-sample transcriptome assembly from RNA-Seq
'''
from taco.lib.optimize import maximize_bisect

__author__ = "Matthew Iyer, Yashar Niknafs, and Balaji Pandian"
__copyright__ = "Copyright 2012-2018"
__credits__ = ["Matthew Iyer", "Yashar Niknafs", "Balaji Pandian"]
__license__ = "MIT"
__version__ = "0.7.3"
__maintainer__ = "Yashar Niknafs"
__email__ = "yniknafs@umich.edu"
__status__ = "Development"


def test_maximize_bisect():
    def f1(x):
        return x
    x, y = maximize_bisect(f1, 0, 100, 0)
    assert x == 100
    assert y == 100

    def f2(x):
        return -(x - 50) ** 2
    x, y = maximize_bisect(f2, 1, 100, 0)
    assert x == 50
    assert y == 0

    def f3(x):
        a = [1, 5, 11, 14, 16, 50, 100, 10000, 5, 4, 3, 2, 1]
        return a[x]
    x, y = maximize_bisect(f3, 1, 12, 0)
    assert x == 7
    assert y == 10000

    def f4(x):
        return x
    x, y = maximize_bisect(f4, 1, 1, 0)
    assert x == 1
    assert y == 1
    x, y = maximize_bisect(f4, 1, 2, 0)
    assert x == 2
    assert y == 2
    x, y = maximize_bisect(f4, 1, 3, 0)
    assert x == 3
    assert y == 3
