"""This module contains functions that compute with Fibonacci numbers.

One goal of this script is to experiment with the different
algorithms to compute Fibonacci numbers. Here are the algorithms:

- fibo_iter: Naive algorithm (use the formula F_n=F_{n-1}+F_{n-2},
but save the values in a list after each call, using the @cache_calls decorator;
- fibo_list: Return a list of the N first fibonacci numbers, using
the naive method, _without_ saving the list after eaach call;
- fibo_cache_list: Same as fibo_list, but save the list using the @cache_list decorator;
- fibonacci: compute the nth Fibonacci number using the powers of quadratic elements;
- fibo: return fibonacci(n);

The algorithm which uses quadratic numbers is by far the most efficient.
"""

from cachetools import cache_calls, cache_list
from quadratic import QuadraticOrder

@cache_calls
def fibo_iter(n):
    """Compute the nth fibonacci number iteratively."""

    if n <= 0: raise ValueError("Invalid value for Fibonacci numbers.")
    if n <= 2: return 1
    a, b = 1, 1
    for _ in range(n - 2):
        c = a + b
        a = b
        b = c
    return b


def fibo_list(N, known = []):
	"""Return a list of the n first Fibonacci numbers.
	
	If the the first k Fibonacci numbers are known for some k<N, the
	list of them can be passed to the function.
	"""
    if N <= 0: raise ValueError("Invalid value for Fibonacci numbers.")

    # Set initial values
    initial_val = [1, 1]

    # If the initial values are not known
    if len(known) < 2:
        known.extend(initial_val[len(known):])

    #If we already know enough values for the precision asked
    if N < len(known):
        return known[:N]

    # Otherwise, we need to extend the list by computing more values
    start = len(known)
    known.extend([0]*(N - start))
    for n in range(start, N):
        known[n] = known[n - 1] + known[n - 2]
    return known

@cache_list
def fibo_cache_list(N, known = []):
	"""Return a list of the n first Fibonacci numbers, but save the list.
	
	If the the first k Fibonacci numbers are known for some k<N, the
	list of them can be passed to the function.
	"""
    if N <= 0: raise ValueError("Invalid value for Fibonacci numbers.")

    # Set initial values
    initial_val = [1, 1]

    # If the initial values are not known
    if len(known) < 2:
        known.extend(initial_val[len(known):])

    #If we already know enough values for the precision asked
    if N < len(known):
        return known[:N]

    # Otherwise, we need to extend the list by computing more values
    start = len(known)
    known.extend([0]*(N - start))
    for n in range(start, N):
        known[n] = known[n - 1] + known[n - 2]
    return known


def fibonacci(n):
    """Compute the nth fibonacci number."""

    if n <= 0:
        raise ValueError("{}: Invalid value for Fibonacci numbers.".format(n))
    phi = QuadraticOrder(5,0,1)
    return (phi**n).b

def fibo(n):
    """Compute the nth fibonacci number."""

    return fibonacci(n)

if __name__ == "__main__":
    import sys
    # if len(sys.argv) > 1:
    #     print(fibonacci(int(sys.argv[-1])))
    N = int(input("Enter an integer: "))
    from time import time
    print("Starting computations...")
    start = time()
    L1 = fibo_list(N)
    print("Without cache: {}".format(time() - start))
    start = time()
    L2 = fibo_cache_list(N)
    print("With cache: {}".format(time() - start))
    print("The length of the list are: {} and {}".format(len(L1),len(L2)))
