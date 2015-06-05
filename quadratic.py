"""Define quadratic fields, quadratic forms and other related objects.

This module defines quadratic fields (QuadraticRat class), quadratic integers
(QuadraticOrder class), quadratic integral ideals (QuadraticIntId class) and
quadratic forms.
"""

import rat, gcd
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
        return len(QuadraticForm.Gaussian_classes_reps(D))
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
    >>> i**0
    1
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
        """Compute the norm of the element.

        The norm of an element x in a quadratic field is defined x*sigma(x),
        where sigma is the non-trivial automorphism of the field.

        Examples:
        >>> x = QuadraticOrder(125,3,5)
        >>> x.norm() == (x*x.conj()).a # Recall that x*x.conj() = a + b*w
        True
        """

        if self.D % 4 == 1:
            return (self.a*self.a + self.a*self.b
                    + self.b*self.b*(1 - self.D)//4)
        return self.a*self.a - self.D//4*self.b*self.b

    def __pow__(self, n):
        """Compute the positive powers of the element.

        Examples:
        >>> w = QuadraticOrder(-3,0,1)
        >>> w**6
        1
        """

        if n==0:
            return QuadraticOrder(self.D,1,0)
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

    __str__ = __repr__


# DEVEL: could add a condition to test if the given Z-module is really
# an integral ideal according to the ideal criterion (see MScThesis...)

# DEVEL: could add a function the returns a simplified basis of the ideal
# i.e. a basis with d = 0

# DEVEL: add a way to represent the ideal.
class QuadraticIntId:
    """Define an integral ideal in the quadratic order of discriminant D.

    Defines an integral ideal I = [a, b + c*w], where a,b,c and d are
    integers and the basis of the order is
    [1,sqrt(D//4)]      if D % 4 == 0 and
    [1,(1 + sqrt(D))/2] if D % 4 == 1

    Note that this is not the same basis as in MScThesis.
    """

    def __init__(self, *args, **kwargs):
        """Define an integral ideal in the order containing the generators.

        An alternative way of constructing an ideal is by specifying directly
        the order containing it (by letting D = discriminant) and the integers
        a,b and c such that the ideal has the form [a, b + c*w] in the standard
        basis of the order of discriminant D.

        Examples:
        >>> I = QuadraticIntId(QuadraticOrder(-4,1,0), QuadraticOrder(-4,0,1))
        >>> I
        [1, sqrt(-1)]
        >>> J = QuadraticIntId(D = -4, abc = (1,0,1))
        >>> J
        [1, sqrt(-1)]
        """

        if 'abc' in kwargs and 'D' in kwargs:
            self.gens = (QuadraticOrder(kwargs['D'],kwargs['abc'][0],0),
                         QuadraticOrder(kwargs['D'],kwargs['abc'][1],kwargs['abc'][2]))
            self.D = kwargs['D']
        elif len(args) == 2: # We have a list of two generators
            assert (isinstance(args[0],QuadraticOrder) and
                    isinstance(args[1],QuadraticOrder) and
                    args[0].D == args[1].D)
            self.gens = (args[0], args[1])
            self.D = args[0].D
        self.simplified_basis = QuadraticIntId.simplify_basis(*self.gens)

    @staticmethod
    def simplify_basis(gen1, gen2):
        """Return a simplified basis of the integral ideal.

        Given an integral ideal I = [gen1, gen2] in a quadratic order, this
        method resturns a couple (alpha,beta) of quadratic integers
        representing the simplified basis [alpha, beta], where alpha = a is an
        integer and beta = b + c*w, i.e. where on of the two generators is an
        integer. This is possible by applying the euclidian algorithm the
        irrational components of the generators.

        Examples:
        >>> x = QuadraticOrder(-4,1,-1)
        >>> y = QuadraticOrder(-4,1,1)
        >>> QuadraticIntId.simplify_basis(x,y)
        (2, 1+sqrt(-1))
        >>> z = QuadraticOrder(-4,1,0)
        >>> QuadraticIntId.simplify_basis(y,z)
        (1, sqrt(-1))
        """

        assert (isinstance(gen1,QuadraticOrder) and
                isinstance(gen2,QuadraticOrder) and gen1.D == gen2.D)
        a, b = gen1.a, gen1.b # gen1 = a + b*w
        c, d = gen2.a, gen2.b # gen2 = c + d*w
        if b < 0: a, b = -a, -b # So b >= 0
        if d < 0: c, d = -c, -d # So d >= 0
        while not b*d == 0:
            if b >= d:
                q = b//d
                a, b = a - q*c, b - q*d
            else:
                q = d//b
                c, d = c - q*a, d - q*b
        if b == 0:
            return (QuadraticOrder(gen1.D, abs(a), 0),
                    QuadraticOrder(gen1.D, c % abs(a), d))
        else:
            return (QuadraticOrder(gen1.D, abs(c), 0),
                    QuadraticOrder(gen1.D, a % abs(c), b))

    def norm(self):
        """Return the norm of the integral ideal.

        That is the index of this ideal in the order containing it. The norm
        of the ideal [a + b*w, c + d*w] is then |a*d - b*c|.

        Examples:
        >>> I = QuadraticIntId(QuadraticOrder(-4,1,1),QuadraticOrder(-4,1,-1))
        >>> I.norm()
        2
        >>> I = QuadraticIntId(D = -4, abc = (1, 0, 1)) # The maximal order
        >>> I.norm()
        1
        """

        return self.simplified_basis[0].a*self.simplified_basis[1].b

    def get_D_K(self):
        return fund_decomp(self.D)[1]

    def get_conductor(self):
        return fund_decomp(self.D)[0]

    def get_disc(self):
        return self.D

    def  get_abc(self):
        """Return (a, b, c) such that the ideal is [a, b + c*w].

        Examples:
        >>> I = QuadraticIntId(D = 20, abc = (1,2,3))
        >>> I.get_abc()
        (1, 0, 3)
        >>> J = QuadraticIntId(QuadraticOrder(-4,1,1),QuadraticOrder(-4,1,-1))
        >>> J.get_abc()
        (2, 1, 1)
        """

        return (self.simplified_basis[0].a, self.simplified_basis[1].a,
                self.simplified_basis[1].b)

    def integral_basis(self):
        return self.gens

    # See quadratic_comp4 for the computations, but replace y by -y...
    # This is the same as sending the form to its inverse under Gauss
    # composition. By doing so, going from forms, to ideal and then back to
    # forms gives the identity, instead of the inversion map...
    def corresponding_form(self):
        """Return the quadratic form corresponding to the ideal.

        The form corresponding to the ideal I = [a, b + c*w] is
        N(a*x -(b + c*50w)*y)/N(I). The minus sign is there to simplify the
        bijection. This is the same as sending the corresponding form to its
        inverse under Gauss composition. By doing so, the map going from forms
        to ideals and the back to forms gives the identity right away (on the
        level of classes, of course...)

        Examples:
        >>> I = QuadraticIntId(QuadraticOrder(-4,1,0),QuadraticOrder(-4,0,1))
        >>> I.corresponding_form()
        x^2+y^2
        >>> J = QuadraticIntId(QuadraticOrder(-3,1,0),QuadraticOrder(-3,0,1))
        >>> J.corresponding_form()
        x^2-xy+y^2

        Note that the last form is equivalent to the form x^2+xy+y^2.
        """

        a, b, c = self.get_abc()
        if self.D % 4 == 1:
            return QuadraticForm(a//c, -(2*b//c + 1),
                                 (b**2 + b*c + (1 - self.D)//4*c**2)//(a*c))
        else:
            return QuadraticForm(a//c, -2*b//c,
                                 (b**2 - self.D//4*c**2)//(a*c))

    # Computations in the basis [1,wf]
    # def corresponding_form(self):
        # if self.D_K % 4 == 0:
            # return QuadraticForm(self.a//self.c, 2*self.b//self.c,
                                 # (self.b**2 - self.D//4*self.c**2)//
                                 # (self.a*self.c))
        # return     QuadraticForm(self.a//self.c, 2*self.b//self.c + self.f,
                                 # (self.b**2 + self.b*self.c*self.f +
                                 #(self.f**2 - self.D//4*self.c**2)//
                                 # (self.a*self.c))

    def __repr__(self):
        return "[{0[0]!s}, {0[1]!s}]".format(self.simplified_basis)

    __str__ = __repr__


# DEVEL: Implement composition.

# TODO: implement the corresponding ideal function for indefinite forms.

# TODO: implement reduction theory for indefinite forms.

# TODO: wrap functions in decorator to test for discriminant.
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
    >>> f.D
    -4
    """

    def __init__(self,a,b,c):
        """Define an integral binary quadratic form.

        Defined by its coefficients a, b and c.
        """

        self._a, self._b, self._c = a, b, c
        self._D = b**2 - 4*a*c

    def __call__(self,x,y):
        """Evaluate the quadratic form at the given point.

        Examples:
        >>> f = QuadraticForm(1,1,1); f
        x^2+xy+y^2
        >>> f(1,1)
        3
        >>> g = QuadraticForm(5,-6,2)
        >>> g(1,1)
        1

        So g is equivalent to the principal form. This also follows from the
        theory since g.D = -4 and there is only one equivalence class of
        forms of discriminant -4.
        """

        return self.a*x**2 + self.b*x*y + self.c*y**2

    def __getitem__(self,n):
        """Return the nth coefficient of self in the basis x^2, xy, y^2.

        Examples:
        >>> f = QuadraticForm(2,3,5); f
        2x^2+3xy+5y^2
        >>> for n in range(3):
        ...    print(f[n])
        ...
        2
        3
        5

        If the key is not an integer, a TypeError is raised.
        If the key is not in the range, a IndexError is raised.
        """

        if not isinstance(n,int):
            raise TypeError("key must be an integer.")
        elif not (0 <= n and n <= 2):
            raise IndexError("Index out of range.")
        else:
            if n == 0:   return self.a
            elif n == 1: return self.b
            else:        return self.c

    def __eq__(self,other):
        """Determine if two quadratic forms are equal (not only equivalent).

        Verifies that the coefficients of the two forms match.

        Examples:
        >>> f = QuadraticForm(1,0,1)
        >>> g = QuadraticForm(5,-6,2)
        >>> f == g # f is not equal to g...
        False
        >>> f == g.class_rep() # But they are equivalent!
        True
        """

        assert isinstance(other,QuadraticForm), "Not two QuadraticForm object."
        return self.a == other.a and self.b == other.b and self.c == other.c

    def __hash__(self):
        return hash((self.a, self.b, self.c))

    @property
    def D(self):
        """Return the discriminant of the binary quadratic form.

        Given a form f = ax^2+bxy+cy^2, this method returns b^2-4*a*c.

        Examples:
        >>> f = QuadraticForm(1,0,1)
        >>> f.D
        -4
        >>> g = QuadraticForm(1, 0, -10)
        >>> g.D
        40
        """

        return self._D

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def c(self):
        return self._c

    @property
    def disc(self):
        return self._D

    def corresponding_ideal(self):
        """Return the integral ideal corresponding to the form.

        The group of proper fractional ideals is isomorphic to the group of
        Gaussian forms (which is in turn isomorphic to the narrow class group
        of the corresponding order). See MScThesis.

        Examples:
        >>> f = QuadraticForm(1,0,1)
        >>> I = f.corresponding_ideal(); I
        [1, sqrt(-1)]
        >>> g = QuadraticForm(1,1,1)
        >>> J = g.corresponding_ideal(); J
        [1, (1 + sqrt(-3))/2]
        >>> h = QuadraticForm(5, -6, 2); h
        5x^2-6xy+2y^2
        >>> h.D
        -4
        >>> H = h.corresponding_ideal(); H
        [5, 3+sqrt(-1)]

        One can also test the bijection between ideals and forms.

        Examples:
        >>> f = QuadraticForm(1,0,1)
        >>> ff = f.corresponding_ideal().corresponding_form()
        >>> f.is_equiv_to(ff)
        True
        >>> classes = QuadraticForm.Gaussian_classes_reps(-71)
        >>> len(classes) == class_number(-71) # Just to check!
        True
        >>> for g in classes:
        ...     gg = g.corresponding_ideal().corresponding_form()
        ...     if not g.is_equiv_to(gg):
        ...         print("There is a problem!!!")
        ...

        """

        if self.is_pos_def():
            if self.D % 4 == 0:
                return QuadraticIntId(
                            QuadraticOrder(self.D,self.a,0),
                            QuadraticOrder(self.D,-self.b//2, 1)
                            )
            else:
                return QuadraticIntId(
                            QuadraticOrder(self.D, self.a, 0),
                            QuadraticOrder(self.D, -(self.b + 1)//2, 1)
                            )
        else:
            pass

    def is_pos_def(self):
        return self.D < 0 and self.a > 0

    def is_almost_reduced(self):
        return abs(self.b) <= self.a and self.a <= self.c

    def is_reduced(self):
        return (self.is_almost_reduced() and
                not ((self.a == self.c and self.b < 0) or (self.a == -self.b)))

    def reduced_class_rep(self):
        """Return a reduced form that is quivalent to the calling form.

        Examples:
        >>> f = QuadraticForm(1,0,1)
        >>> f.reduced_class_rep()
        x^2+y^2
        >>> g = QuadraticForm(5,-6,2)
        >>> g.D
        -4
        >>> g.reduced_class_rep()
        x^2+y^2
        >>> h = QuadraticForm(1,1,1)
        >>> h.reduced_class_rep()
        x^2+xy+y^2
        """

        a, b, c = self.a, self.b, self.c
        D = self.D
        while a > c or abs(b) > a or (abs(b) == a and b < 0) or (a == c and b < 0):
            if a > c:
                a, c = c, a
                b = -b
            elif abs(b) > a:
                # solve |b+2*a*n| <= a for n
                if b < 0: n = int((a - b)/(2*a))
                else:     n = int((-a - b)/(2*a))
                b = b + 2*a*n
                c = (b**2 - D)//(4*a)
            # At this point he form is almost reduced
            elif abs(b) == a:
                b = abs(b)
            else: # a == c
                b = abs(b)
        return QuadraticForm(a, b, c)

    def class_rep(self):
        """Return the reduced representative in the class of self.

        Given a positive definite form f, return the unique reduced form in the
        SL2(Z) equivalence class of f. This method does not modify f.

        Examples:
        >>> f = QuadraticForm(5,-6,2)
        >>> f.class_rep()
        x^2+y^2
        """

        return self.reduced_class_rep()

    def is_equiv_to(self,other):
        """Determine if self is SL2(Z)-equivalent to other.

        Examples:
        >>> f = QuadraticForm(1,0,1)
        >>> g = QuadraticForm(5,-6,2)
        >>> f.is_equiv_to(g)
        True
        >>> g.is_equiv_to(f)
        True
        >>> f == g
        False
        """

        assert (isinstance(other, QuadraticForm) and self.D == other.D), "Not two QF with same disc."
        if self.D < 0:
            return self.class_rep() == other.class_rep()
        else:
            pass

    def right_neighbor(self):
        """Return the right neighbor of the quadratif form.

        The right neighbor (aa, bb, cc) of an indefinite quadratic form
        (a, b, c) is uniquely determined by the three conditions:
        1) aa = c
        2) b + bb % 2aa == 0 and sqrt(D) - 2|aa| < bb < sqrt(D)
        3) bb**2 - 4*aa*cc = D

        Examples: (an example in Flath, Ch. 4, Sec. 6)
        >>> f = QuadraticForm(2, 8, 3)
        >>> for _ in range(8):
        ...     print(f)
        ...     f = f.right_neighbor()
        ...
        2x^2+8xy+3y^2
        3x^2+4xy-2y^2
        -2x^2+4xy+3y^2
        3x^2+2xy-3y^2
        -3x^2+4xy+2y^2
        2x^2+4xy-3y^2
        -3x^2+2xy+3y^2
        3x^2+4xy-2y^2
        >>> g = QuadraticForm(1, 0, -10)
        >>> for _ in range(5):
        ...     print(g)
        ...     g = g.right_neighbor()
        ...
        x^2-10y^2
        -10x^2+y^2
        x^2+6xy-y^2
        -x^2+6xy+y^2
        x^2+6xy-y^2
        """

        a, b, c = self[0], self[1], self[2]
        tau = (b + sqrt(self.D))/(2*c)
        aa = c
        # Solve b + bb = 2cd for bb s.t sqrt(self.D - |2c|<bb<sqrt(...)
        if c > 0:
            if tau < 0: d = int(tau - 1)
            else:       d = int(tau)
        else:
            if tau < 0: d = int(tau)
            else:       d = int(tau + 1)
        bb = 2*c*d - b
        return QuadraticForm(aa, bb, (bb**2 - self.D)//(4*aa))

    def is_primitive(self):
        return gcd.gcd(self.a,self.b,self.c) == 1

    @staticmethod
    def almost_reduced_forms(D):
        assert is_discriminant(D) and D<0, (str(D) + " is not a discriminant.")
        if D < 0:
            M = int((abs(D)/3)**0.5)
            return [QuadraticForm(a,b,(b**2-D)//(4*a)) for a in range(1,M+1)
                    for b in range(-a,a+1) if (b**2-D)%(4*a) == 0]

    @staticmethod
    def reduced_classes_reps(D):
        if D < 0:
            return QuadraticForm.reduced_forms(D)
        else:
            reps = []
            S = set(QuadraticForm.reduced_forms(D))
            while not len(S) == 0:
                f0 = next(iter(S))
                reps.append(f0)
                S.discard(f0)
                right = f0.right_neighbor()
                while not right == f0:
                    S.discard(right)
                    right = right.right_neighbor()
            return reps

    @staticmethod
    def reduced_forms(D):
        if D < 0:
            return [f for f in QuadraticForm.almost_reduced_forms(D)
                    if f.is_reduced()]
        else: # We need 0 < b < sqrt(D) and sqrt(D) - b < 2|a| < sqrt(D) + b
            r = sqrt(D)
            Lpos = [QuadraticForm(a, b, (b**2 - D)//(4*a))
                  for b in range(1,int(r) + 1)
                  for a in range(int((r-b)/2) + 1, int((r+b)/2) + 1)
                  if (b**2 - D) % (4*a) == 0]
            Lneg = [QuadraticForm(a, b, (b**2 - D)//(4*a))
                  for b in range(1,int(r))
                  for a in range(int((-r-b)/2), int((-r+b)/2))
                  if (b**2 - D) % (4*a) == 0]
            return Lpos + Lneg


    @staticmethod
    def classes_reps(D):
        assert is_discriminant(D), (str(D) + " is not a discriminant.")
        return QuadraticForm.reduced_classes_reps(D)

    @staticmethod
    def Gaussian_classes_reps(D):
        assert is_discriminant(D), ("Not a negative discriminant.")
        return [f for f in QuadraticForm.classes_reps(D)
                if f.is_primitive()]

    def __repr__(self):
        """Represent the binary quadratic form as a polynomial.

        Examples:
        >>> g = QuadraticForm(1, 0, -10); g
        x^2-10y^2
        """

        s = ""
        for t in zip([self.a,self.b,self.c],["x^2","xy","y^2"]):
            if t[0] == 0:    pass
            elif t[0] == 1:  s += '+' + t[1]
            elif t[0] == -1: s += '-' + t[1]
            else:            s += '+' + str(t[0]) + t[1]
        if s[0] == '+': s = s[1:]
        return s.replace('+-','-')

if __name__ == "__main__":
    import doctest
    doctest.testmod()
