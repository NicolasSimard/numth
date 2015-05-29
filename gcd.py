def gcd(*num):
    """Return the gcd of the integers.
    
    Example:
    >>> gcd(1,2)
    1
    >>> gcd(4,2)
    2
    >>> gcd(4,-2)
    2
    >>> gcd(-4,-2)
    2
    >>> gcd(2,4,6,8,10)
    2
    >>> gcd(gcd(10,100),3) == gcd(10, 100, 3)
    True
    """
    
    assert not len(num) == 0, " I need at least one number to compute the gcd!"
    value = num[0]
    for x in num[1:]:
        value = _intern_gcd(value,x)
    return value
    
    
def _intern_gcd(m,n):
    while not n == 0:
        m, n = n, m % n
    return abs(m)
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import sys
    print(gcd(*list(map(int,sys.argv[1:]))))