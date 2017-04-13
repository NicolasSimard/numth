def binom(n,k):
	"""This module contains a function that computes the binomial coefficients.

	The function binom(n,k) returns the binomial coefficient n!/(k!(n-k)!). The
	binomial terms are computed recursively.

	Examples:
	>>> binom(4,2)
	6
	>>> binom(10,1)
	10
	>>> sum(binom(10,k) for k in range(11)) == 2**10
	True
	"""
    assert 0 <= k and k <= n
    binoms = [[1]]
    for i in range(1,n+1):
        binoms.append([0]*(i+1))
    if k == 0 or k == n:
        binoms[n][k] = 1
        return 1
    if not binoms[n][k] == 0:
        return binoms[n][k]
    binoms[n][k] = binom(n-1,k-1) + binom(n-1,k)
    return binoms[n][k]

if __name__ == "__main__":
    import sys, doctest
    if len(sys.argv) == 3:
        print(binom(int(sys.argv[1]),int(sys.argv[2])))
    doctest.testmod()