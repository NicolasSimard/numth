from cachetools import cache_calls, cache_list
from quadratic import QuadraticOrder

@cache_calls
def fibo_iter(n):
    """Compute the nth fibonacci number iteratively.
    """

    if n <= 0: raise ValueError("Invalid value for Fibonacci numbers.")
    if n <= 2: return 1
    a, b = 1, 1
    for _ in range(n - 2):
        c = a + b
        a = b
        b = c
    return b


def fibo_list(N, known = []):
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
    """Compute the nth fibonacci number.
    """

    if n <= 0:
        raise ValueError("{}: Invalid value for Fibonacci numbers.".format(n))
    phi = QuadraticOrder(5,0,1)
    return (phi**n).b

def fibo(n):
    """Compute the nth fibonacci number.
    """

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
