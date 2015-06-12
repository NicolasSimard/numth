from cachetools import cache_calls
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
    if len(sys.argv) > 1:
        print(fibonacci(int(sys.argv[-1])))