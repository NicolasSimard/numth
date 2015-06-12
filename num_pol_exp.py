"""Compute the coefficients of the numerical polynomial in the standard basis.

The polynomials (x|n)=x(x-1)...(x-n+1)/n! form the standard basis for the
group of numerical polynomials, i.e. the polynomials that take integral values
at integers. An example of numerical polynomial is x(x+1)/2. This code
computes the expansion of a numerical polynomial in this basis.

-input: a0 a1 a2 ... an d, where the numerical polynomial is
(a0 + a1*x + a2*x^2 + ... + an*x^n)/d, with an != 0.

-output: [c0,c1,...,cn], the coefficients of the expansion

Examples:
>>> num_pol_exp([0, -1, 1, 2]) # The polynomial x(x-1)/2
[0, 0, 1]
>>> num_pol_exp([0, 1, 0, -7, 0, 21, 21, 6, 42]) # Sum of 6th powers
[0, 1, 63, 602, 2100, 3360, 2520, 720]
"""
import binom

def num_pol_exp(coeff):
    f=[0]*len(coeff)
    d = coeff.pop() # change to coeff[-1]
    for k in range(len(coeff)):
        power = 1
        for c in coeff: # change to coeff[:-1]
            f[k] += c*power
            power *= k
        f[k] /= d
        f[k] = int(f[k])
    # computing the coefficients in the basis, up to the precision
    basis = []
    for n in range(len(coeff)):
        basis.append(
                sum([(-1)**(n-k)*binom.binom(n,k)*f[k] for k in range(n+1)]))
    return basis

if __name__ == "__main__":
    import sys, doctest
    if len(sys.argv) > 1:
        print(num_pol_exp(list(map(int,sys.argv[1:]))))
    doctest.testmod()