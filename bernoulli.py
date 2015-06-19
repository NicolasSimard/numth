from rat import Rat

def bernoulli(n):
    """Return the nth Bernoulli number (as a rational number).

    The algorithm is taken from Experimental Number Theory by F. R.
    Villegas (the algorithm itself is apparently due to D. Zagier).

    Examples:
    >>> for n in range(20):
    ...     print(bernoulli(n))
    1
    -1/2
    1/6
    0
    -1/30
    0
    1/42
    0
    -1/30
    0
    5/66
    0
    -691/2730
    0
    7/6
    0
    -3617/510
    0
    43867/798
    0
    """

    x, res, s, c = Rat(0), Rat(0), Rat(0), Rat(-1)
    for k in range(1, n+2):
        c *= 1 - Rat(n + 2)/k
        s += x**n
        x += 1
        res += c*s/k
    return res

if __name__ == "__main__":
    import sys, doctest
    if len(sys.argv) > 1:
        print(list(map(bernoulli,map(int,sys.argv[1:]))))
    doctest.testmod()
