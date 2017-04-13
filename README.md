# numth
Number Theory Python 3 Package.

## Introduction
The initial objective of this repo was to do computation in relation with quadratic fields (quadratic numbers, quadratic orders,
quadratic integers, ideals in quadratic fields, quadratic forms, etc.). At the time, I was falling in love with python and so 
it made sense to implement everything I needed from scratch! In the meantime, I got familiar with PARI/GP and I switched to it 
and the ENT repo was born. Althought Python is not ther ideal language to do advanced number theory (one should use PARI/GP, Sage or Magma),
this repo can still be useful to students who want to de basic number theory computations in Python. Eventually, the project 
became a good testing ground for slightly more advanced features of the Python language (decorators, docstrings, doctest, etc.)

## Content of the repo
Here is a small description of the scripts contained in the numth repo. For more details, read the docstrings or use the help() function!

### Main modules
- quadratic.py: This script defines four classes: QuadraticRat, QuadraticOrder, QuadraticIntId and QuadraticForm.
- rat.py: This script defines the class Rat, to work with rational numbers.

### Auxilary number theory modules
- bernoulli.py: defines the function `bernoulli` to compute the Bernoulli numbers.
- gcd.py: defines the function `gcd` to compute the gcd of arbitrarily many numbers.
- binom.py: defines the function `binom` to compute the binomial coefficients.
- primes.py: defines a function `prime_list` which returns a list of the first N prime numbers (using sieving).
- fibonacci.py: defines many functions to compute the Fibonacci numbers. The most efficient one is `fibonacci`.
- num_pol_exp.py: defines the function `num_pol_exp` to compute the expansion of a numerical polynomial in terms of the standard basis.

### Auxilary tools
- cachetools.py: defines two decorators, namely `cache_list` and `cache_calls`, which save calls or a list of function calls.

The auxilary number theory modules can also be called directly from the command line. For example, the command `python fibonacci.py 5` returns
the 5th Fibonacci number (assuming that the command `python` runs Python 3 on your system and that you are in the numth repo).

## Testing
Most of these modules are doctest compatible. Actually, calling the Auxiliary number theory modules from the __main__ usually
runs doctest on the module. 