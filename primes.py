def prime_list(n): #primes < n
    """Return a list of all primes strictly less than n.

    Examples:
    >>> prime_list(5)
    [2, 3]
    >>> prime_list(10)
    [2, 3, 5, 7]
    """
    is_prime = [True] * n
    for i in range(3,int(n**0.5)+1,2):
        if is_prime[i]:
            is_prime[i*i::2*i]=[False]*((n-i*i-1)//(2*i)+1)
    return [2] + [i for i in range(3,n,2) if is_prime[i]]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import sys
    if len(sys.argv) > 1:
        print(prime_list(int(sys.argv[1])))
