"""Define quadratic fields, quadratic forms and other related objects.

This module defines quadratic fields (QuadraticRat class), quadratic integers
(QuadraticInt class), quadratic integral ideals (QuadraticIntId class) and 
quadratic forms.

To import the module, simply enter the following line in the console:
>>> from numth import quadratic
"""

from numth import rat, gcd
from math import sqrt
from time import time

def is_discriminant(D):
    return D % 4 == 0 or D % 4 == 1
    
def class_number(D):
    """Compute the class number of the quadratic field of discriminant D.
    
    This function uses the theory of binary quadratic forms to compute the
    class number of the quadratic order of discriminant D.
    
    Example:
    >>> for D in [D for D in range(-1,-200,-1) if is_discriminant(D)]:
    ...     if class_number(D) == 1:
    ...         print(D)
    ...
    -3
    -4
    -7
    -8
    -11
    -12
    -16
    -19
    -27
    -28
    -43
    -67
    -163
    
    As expected!
    """
    assert is_discriminant(D)
    if D < 0:
        return len(QuadraticForm.Gaussian_forms(D))
    else:
        return None

def fund_decomp(D):
    """Return a decomposition D = f**2 D_K with D_K a fundamental discriminant.
    
    Example:
    >>> fund_decomp(-4)
    (-4, 1)
    >>> fund_decomp(-8)
    (-8, 1)
    >>> fund_decomp(-12)
    (-3, 2)
    >>> fund_decomp(20)
    (5, 2)
    """
    assert is_discriminant(D), "Not a discriminant!"
    f = int(sqrt(abs(D)))
    while not D % f**2 == 0: f -= 1
    # At this point, D = f**2 sf_D, where sf_D % 4 == 1,2 or 3 is square-free.
    sf_D = D//(f**2)
    if sf_D % 4 == 1:
        return (sf_D, f)
    else:
        return (4*sf_D, f//2)

def is_fundamental(D):
    """Determine if the integer D is a fundamental discriminant.
    
    Example:
    >>> is_fundamental(-4)
    True
    >>> is_fundamental(-8)
    True
    >>> is_fundamental(-12)
    False
    >>> is_fundamental(12)
    True
    >>> is_fundamental(20)
    False
    """
    return is_discriminant(D) and fund_decomp(D)[1] == 1
       
class QuadraticRat:
    """Represent an element in a quadratic field of fixed discriminant.
    
    Represents an element of the unique quadratic extension of Q of 
    discriminant D_K, where D_K is a fundamental discriminant. The computations 
    are done in the basis 
    [1,sqrt(D_K//4)]  if D % 4 == 0 and
    [1,sqrt(D_K)]     if D % 4 == 1
    
    Note that this is the standard basis (the square root of the square-free
    part of D_K) since D_K is fundamental.
    
    Examples:
    >>> a = QuadraticRat(-4,rat.Rat(1,2),rat.Rat(3,2)); a
    1/2+3/2*sqrt(-1)
    >>> b = QuadraticRat(-4,rat.Rat(4,2),rat.Rat(7,2)); b
    2+7/2*sqrt(-1)
    >>> (a + b)*(a - b) == a**2 - b**2
    True
    >>> w = QuadraticRat(-3,rat.Rat(-1,2),rat.Rat(1,2))
    >>> w.norm()
    1
    >>> w**3
    1
    >>> w**(-1) == w.conj()
    True
    >>> a*a.conj()/a.norm()
    1
    >>> a.conj()/a.norm() == a**(-1)
    True
    >>> a/b
    5/13+1/13*sqrt(-1)
    >>> a + b + 1 + rat.Rat(3,5)
    41/10+5*sqrt(-1)
    >>> 2*b
    4+7*sqrt(-1)
    """

    def __init__(self,D_K,a=1,b=0):
        """Construct a quadratic element in the basis [1,sqrt(D_K)].
        
        Note: D_K is square-free if it is congruent to 1 mod 4 and D_K//4 is 
        square-free if D_K is congruent to 0 mod 4.
        
        Example:
        >>> i = QuadraticRat(-4, 0, 1)
        >>> i**2
        -1
        """
        assert is_fundamental(D_K), (str(D_K) + " is not fundamental!")
        self.D_K = D_K
        if self.D_K % 4 == 0:
            self.D0 = self.D_K//4
        else:
            self.D0 = self.D_K
        self.a = rat.Rat(a)
        self.b = rat.Rat(b)

    def __add__(self, other):
        if isinstance(self,other.__class__) and self.D_K == other.D_K:
            return QuadraticRat(self.D_K,self.a + other.a, self.b + other.b)
        elif isinstance(other,rat.Rat) or isinstance(other,int):
            return QuadraticRat(self.D_K,self.a + other,self.b)
        else:
            return NotImplemented

    __radd__ = __add__

    def __neg__(self):
        return QuadraticRat(self.D_K,-self.a,-self.b)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self,other):
        return -self + other

    def __eq__(self,other):
        if isinstance(self,other.__class__) and self.D_K == other.D_K:
            return self.a == other.a and self.b == other.b
        elif isinstance(other,rat.Rat) or isinstance(other,int):
            return self.a == other and self.b == 0
        else:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(self,other.__class__) and self.D_K == other.D_K:
            return QuadraticRat(self.D_K,self.a*other.a 
                                + self.D0*self.b*other.b, 
                                self.b*other.a + self.a*other.b)
        elif isinstance(other,rat.Rat) or isinstance(other,int):
            return QuadraticRat(self.D_K,self.a * other,self.b * other)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def conj(self):
        return QuadraticRat(self.D_K,self.a,-self.b)

    def norm(self):
        return self.a**2 - self.D0*self.b**2

    def __truediv__(self,other):
        if isinstance(self,other.__class__) and self.D_K == other.D_K:
            return QuadraticRat(self.D_K,(self.a*other.a 
                                - self.D0*self.b*other.b)/other.norm(),
                                (other.a*self.b - self.a*other.b)/other.norm())
        elif isinstance(other,rat.Rat) or isinstance(other,int):
            return self*(rat.Rat(other)**-1)
        else:
            return NotImplemented

    def __rtruediv__(self,other):
        if isinstance(other,rat.Rat) or isinstance(other,int):
            return QuadraticRat(self.D_K,(self.a*other)/self.norm(),
                                (-self.b*other)/self.norm())
        else:
            return NotImplemented

    def __pow__(self, n):
        if n==0:
            return QuadraticRat(self.D_K,1,0)
        elif n < 0:
            return 1/self.__pow__(-n)
        else:
            temp = self.__pow__(n//2)
            if n%2==0:
                return temp*temp
            else:
                return self*temp*temp

    def __eq__(self,other):
        return self.a == other.a and self.b == other.b

    def __repr__(self):
        s = ""
        if not self.a == 0:
            s += str(self.a)
            if not self.b == 0: s += "+"
        if not self.b == 0: s += str(self.b) + "*sqrt(" + str(self.D0) + ")"
        if not len(s) == 0: return s.replace('+-','-').replace('1*','')
        return "0"

        
# DEVEL: allow interaction with integers. This doesn't seem to affect the speed
# too much.
class QuadraticOrder:
    """Represent the quadratic order of a given discriminant.
    
    Elements in the quadratic order of discriminant D are represented by their
    integral components in the basis
    [1,sqrt(D//4)]      if D % 4 == 0 and
    [1,(1 + sqrt(D))/2] if D % 4 == 1
    See quadratic_comp3.jpg for remarks about this choice of basis. Note that 
    the order Z[sqrt(d)] has discriminant 4d and so our choice of basis gives
    [1,sqrt(d)], which is the natural choice. Note also that there is no need
    to factor the discriminant into its fundamental and square part to use this
    basis.
    
    Example:
    >>> phi = QuadraticOrder(5,0,1); phi
    (1 + sqrt(5))/2
    >>> (phi**25).b # The 25th Fibonacci number
    75025
    >>> i = QuadraticOrder(-4,0,1)
    >>> i**2
    -1
    >>> one = QuadraticOrder(-4,1,0)
    >>> i*(one - i)**2 # 2 ramifies in Z[i]
    2
    """

    def __init__(self,D,a=1,b=0):
        assert is_discriminant(D),(str(D) + " is not a discriminant!")
        self.D = D
        self.a = a
        self.b = b

    def __add__(self, other):
        return QuadraticOrder(self.D,self.a + other.a,self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return QuadraticOrder(self.D,-self.a,-self.b)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self,other):
        return -self + other

    def __eq__(self,other):
        return self.a == other and self.b == other.b

    def __mul__(self, other):
        if self.D % 4 == 1:
            return QuadraticOrder(self.D,
                              self.a*other.a + (self.D - 1)//4*self.b*other.b,
                              self.b*other.a + self.a*other.b + self.b*other.b)
        else:
            return QuadraticOrder(self.D,
                                  self.a*other.a + self.D//4*self.b*other.b,
                                  self.b*other.a + self.a*other.b)

    __rmul__ = __mul__

    def conj(self):
        if self.D % 4 == 1:
            return QuadraticOrder(self.D,self.a + self.b,-self.b)
        return QuadraticOrder(self.D,self.a,-self.b)

    def norm(self):
        if self.D % 4 == 1:
            return self.a*self.a + self.a*self.b + self.b*self.b*(1 - self.D)//4
        return self.a*self.a - self.D//4*self.b*self.b

    def __pow__(self, n):
        if n==0:
            return QuadraticInt(self.D,1,0)
        elif n == 1: # Adding this line saves a bit of time
            return self
        else: # Creating the variable temp saves A LOT of time in computations.
            temp = self.__pow__(n//2)
            if n%2==0:
                return temp*temp
            else:
                return self*temp*temp
    # NOTE: it is not computationaly faster to use ternary exponentiation
    # instead of binary (i.e. use 3 instead of 2). In fact it is empirically
    # 50% slower!

    def __eq__(self,other):
        return self.a == other.a and self.b == other.b

    def __repr__(self):
        s = ""
        if not self.a == 0:
            s += str(self.a)
            if not self.b == 0: s += "+"
        if not self.b == 0:
            if self.D % 4 == 1:
                s += str(self.b) + "*" + "(1 + sqrt(" + str(self.D) + "))/2"
            else:
                s += str(self.b) + "*" + "sqrt(" + str(self.D//4) + ")"
        if not len(s) == 0: return s.replace('+-','-').replace('1*','')
        return "0"
    

# DEVEL: could add a condition to test if the given Z-module is really
# an integral ideal according to the ideal criterion (see MScThesis...)

# DEVEL: could add a function the returns a simplified basis of the ideal
# i.e. a basis with d = 0.

# DEVEL: add a way to represent the ideal.
class QuadraticIntId:
    """Define an integral ideal in the quadratic order of discriminant D.
    
    Defines an integral ideal I = [a, b + c*w], where a,b,c and d are
    integers and the basis of the order is
    [1,sqrt(D//4)]      if D % 4 == 0 and
    [1,(1 + sqrt(D))/2] if D % 4 == 1
    
    Note that this is not the same basis as in MScThesis. 
    """

    def __init__(self,D,a,b,c):
        # Defines the ideal [a, b + c*w] in the unique order of discriminant D.
        assert is_discriminant(D), (str(D_K) + " is not a discr!")
        self.D = D
        # Those sign changes are made so that I has positive orientation, i.e.
        # Im((b+c*w_f)/a) > 0
        self.a = abs(a)
        if c < 0:
            self.b = -b
            self.c = -c
        else:
            self.b = b
            self.c = c

    def norm(self):
        return a*c

    def get_D_K(self):
        return self.D_K

    def get_conductor(self):
        return self.f

    def get_disc(self):
        return self.f**2*self.D_K

    def integral_basis(self):
        return (QuadraticInt(self.D_K,self.a,0),
                QuadraticInt(self.D_K,self.b,self.c*self.f))

    # See quadratic_comp4 for the computations
    def corresponding_form(self):
        """Return the quadratic form corresponding to the ideal.
        
        The form corresponding to the ideal I = [a, b + c*w] is
        N(a*x + (b + c*w)*y)/N(I).
        
        Examples: 
        >>> I = QuadraticIntId(-4,1,0,1)
        >>> I.corresponding_form()
        x^2+y^2
        >>> J = QuadraticIntId(-3,1,0,1)
        >>> J.corresponding_form()
        x^2+xy+y^2
        """
        if self.D % 4 == 1:
            return QuadraticForm(self.a//self.c, 2*self.b//self.c + 1, 
                                 (self.b**2 + self.b*self.c 
                                 + (1 - self.D)//4*self.c**2)//(self.a*self.c))
        else:
            return QuadraticForm(self.a//self.c, 2*self.b//self.c,
                                 (self.b**2 - self.D//4*self.c**2)//
                                 (self.a*self.c))
                                 
    # Computations in the basis [1,wf]        
    # def corresponding_form(self):
        # if self.D_K % 4 == 0:
            # return QuadraticForm(self.a//self.c, 2*self.b//self.c,
                                 # (self.b**2 - self.get_disc()//4*self.c**2)//
                                 # (self.a*self.c))
        # return     QuadraticForm(self.a//self.c, 2*self.b//self.c + self.f, 
                                 # (self.b**2 + self.b*self.c*self.f + 
                                 # (self.f**2 - self.get_disc())//4*self.c**2)//
                                 # (self.a*self.c))

    def __repr__(self):
        basis = self.integral_basis()
        return "[" + basis[0].__str__() + "," + basis[1].__str__() + "]"

        
# DEVEL: Implement composition
class QuadraticForm:
    """Define an integral binary quadratic form.
    
    That is a polynomial of the form a*x**2 + b*x*y + c*y**2. This class
    defines the discriminant. Quadratic form objects are also callable. See
    examples below.
    
    Examples: 
    >>> f = QuadraticForm(1,0,1)
    >>> f
    x^2+y^2
    >>> f(1,2)
    5
    >>> f.disc()
    -4
    """
    
    def __init__(self,a,b,c):
        """Define an integral binary quadratic form.
        
        Defined by its coefficients a, b and c.
        """
        self.a, self.b, self.c = a, b, c
        self._fund_decomp = fund_decomp(self.disc())

    def disc(self):
        return b**2 - 4*a*c

    def __repr__(self):
        s = ""
        for t in zip([self.a,self.b,self.c],["x^2","xy","y^2"]):
            if t[0] == 0:   pass
            elif t[0] == 1: s += '+' + t[1]
            else:           s += '+' + str(t[0]) + t[1]
        if s[0] == '+': s = s[1:]
        return s.replace('+-','-')

    def disc(self):
        return self.b**2 - 4*self.a*self.c

    def __call__(self,x,y):
        return self.a*x**2 + self.b*x*y + self.c*y**2
        
    def corresponding_ideal(self):
        """Return the integral ideal corresponding to the form.
        
        The group of proper fractional ideals is isomorphic to the group of
        Gaussian forms (which is in turn isomorphic to the narrow class group
        of the corresponding order). See MScThesis.
        """

        
    def is_reduced(self):
        assert self.disc() < 0, "Not a definite form!"
        return not ((self.a == self.c and self.b < 0) or (self.a == -self.b))
    
    def is_primitive(self):
        return gcd.gcd(self.a,self.b,self.c) == 1
    
    @staticmethod
    def almost_reduced_forms(D):
        assert is_discriminant(D), (str(D) + " is not a discriminant.")
        if D < 0: 
            M = int((abs(D)/3)**0.5)
            return [QuadraticForm(a,b,(b**2-D)//(4*a)) for a in range(1,M+1) 
                    for b in range(-a,a+1) if (b**2-D)%(4*a) == 0]
        else:
            return None
            
    @staticmethod
    def representatives(D):
        assert is_discriminant(D), (str(D) + " is not a discriminant.")
        if D < 0:
            return QuadraticForm.reduced_representatives(D)
        else: 
            return None
        
    @staticmethod
    def reduced_representatives(D):
        assert D < 0, ("The discriminant is not negative")
        return [f for f in QuadraticForm.almost_reduced_forms(D) 
                if f.is_reduced()]
            
    @staticmethod
    def reduced_forms(D):
        return QuadraticForm.reduced_representatives(D)
                
    @staticmethod
    def Gaussian_forms(D):
        assert D < 0 and is_discriminant(D)
        return [f for f in QuadraticForm.reduced_representatives(D)
                if f.is_primitive()]
                

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    import sys
    for D in sys.argv[1:]:
        print(D,":",QuadraticForm.class_number(int(D)))
    # print("Testing with Fibonacci:")
    # phi = QuadraticRat(5,Rat(1,2),Rat(1,2))
    # print("The 1000th Fibonacci number is:",2*(phi**25).b)
    # phiRat = QuadraticRat(5,rat.Rat(1,2),rat.Rat(1,2))
    # phiInt = QuadraticInt(5,0,1)
    # N = int(input("Enter an integer: "))
    # start = time()
    # print("Computing the",N,"th Fibonacci number with QuadraticInt: ")
    # x1 = (phiInt**N).b
    # print("Done! ",time() - start)
    # start = time()
    # print("Computing the",N,"th Fibonacci number with QuadraticRat: ")
    # x2 = (phiRat**N).b
    # print("Done! ",time() - start)
    # print("Are they equal: ", x1 == 2*x2)
    
    # print("Testing with random quadratic elements:")
    # xRat = QuadraticRat(-8,rat.Rat(2),rat.Rat(1,2))
    # xInt = QuadraticInt(-8,2,1)
    # print("The first is:", xRat, "and the second is",xInt)
    # N = int(input("Enter an integer: "))
    # start = time()
    # print("Computing the",N,"th power of xInt: ")
    # x1 = xInt**N
    # print("Done! ",time() - start)
    # start = time()
    # print("Computing the",N,"th power of xRat: ")
    # x2 = xRat**N
    # print("Done! ",time() - start)
    # print("Are they equal: ", x1.b == 2*x2.b and x1.a == x2.a)