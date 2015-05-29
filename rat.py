import gcd

def sign(x):
    if x >= 0: return 1
    return -1
    
# DEVEL: implement comparaison operators
    
class Rat:
    """Represents the quotient of two integers.
    
    Rational numbers can be added, added, divided, etc. They can also
    interact with integers.
    
    Examples:
    >>> from numth import rat
    >>> a = rat.Rat(7)
    >>> a == 7
    True
    >>> b = rat.Rat(1,3); b
    1/3
    >>> c = rat.Rat(15,6); c
    5/2
    >>> b + c
    17/6
    >>> (b - c)*(b + c) == b**2 - c**2
    True
    >>> 1/(1/a) == a
    True
    >>> (b**(-2))**(-3) == b**6
    True
    >>> a*a**(-1) == 1
    True
    >>> b += 1
    >>> b
    4/3
    
    Note that 1/1/a would not work, since 1/1 is a float and a float cannot be
    divided by a rational number.
    
    Conversion works as expected.
    >>> int(c)
    2
    >>> float(c)
    2.5
    """
    
    def __init__(self, num, denom = 1):
        if not isinstance(num,Rat) and not isinstance(denom,Rat):
            num, denom = int(num),int(denom)
            d = gcd.gcd(abs(num), abs(denom))
            # The sign is always in the numerator
            self.num = sign(num*denom)*abs(num)//d
            self.denom = abs(denom)//d
        else:
            fraction = num/denom
            self.num = fraction.num
            self.denom = fraction.denom
        
    def __add__(self, other):
        if isinstance(self,other.__class__):
            return Rat(self.num*other.denom + self.denom*other.num,
                       self.denom*other.denom)
        elif isinstance(other,int):
            return Rat(self.num + self.denom*other, self.denom)
        else:
            return NotImplemented
           
    # If an object X calls __add__ with a Rat object y, the method x.__add__(y)
    # will certainly not be implemented. Python will then try to call
    # y.__radd__(x). Since x+y=y+x, we define __radd__ = __add__.
    __radd__ = __add__

    def __mul__(self, other):
        if isinstance(self,other.__class__):
            return Rat(self.num*other.num, self.denom*other.denom)
        elif isinstance(other,int):
            return Rat(self.num*other, self.denom)
        else:
            return NotImplemented
            
    # See the comment above __radd__.
    __rmul__ = __mul__

    def __neg__(self):
        return Rat(-self.num, self.denom)

    def __abs__(self):
        return Rat(abs(self.num),abs(self.denom))
    
    def sgn(self):
        return sign(self.num)
        
    def __sub__(self, other):
        return self + -other
        
    def __rsub__(self,other):
        return -self + other

    def __eq__(self, other):
        if isinstance(self,other.__class__):
            return self.num*other.denom == self.denom*other.num
        elif isinstance(other,int):
            return self.num == self.denom*other
        else:
            return NotImplemented

    def __truediv__(self,other):        
        if isinstance(self,other.__class__):
            return Rat(self.num*other.denom, self.denom*other.num)
        elif isinstance(other,int):
            return Rat(self.num, self.denom*other)
        else:
            return NotImplemented
            
    def __rtruediv__(self,other):
        if isinstance(other,int):
            return Rat(self.denom*other, self.num)
        else:
            return NotImplemented
            
    def __floordiv__(self,other):        
        if (isinstance(self,other.__class__) and self.denom == 1 
           and other.denom == 1):
            return self.num//other.num
        elif isinstance(other,int) and self.denom == 1:
            return self.num//other
        else:
            return NotImplemented
            
    def __rfloordiv__(self,other):
        if isinstance(other,int) and self.denom == 1:
            return other//self.num
        else:
            return NotImplemented
            
    def __pow__(self,n):
        if isinstance(n, int):
            if n > 0:
                return Rat(self.num**n,self.denom**n)
            elif n < 0:
                return Rat(self.denom**(-n),self.num**(-n))
            else:
                return Rat(1)
        else: 
            return NotImplemented
    
    def __int__(self):
        return int(self.num/self.denom)
        
    def __float__(self):
        return self.num/self.denom
        
    # Since __str__ isn't defined, __str__ = __repr__ by default
    def __repr__(self):
        if self.denom == 1:
            return "{0}".format(str(self.num))
        else:
            return "{0}/{1}".format(str(self.num),str(self.denom))
            
            
if __name__ == "__main__":
    a = Rat(1,2)
    b = Rat(15,6)
    print("Testing addition:")
    print(a,"+",b,"is",a+b)
    print("1+1/2 is",1+a)
    print("1/2+1 is",a+1)
    print("Testing multiplication:")
    print(a,"*",b,"is",a*b)
    print("5*1/2 is",5*a)
    print("1/2*5 is",a*5)
    print(-b,"*",4,"is",-b*4)
    print("Testing subtraction:")
    print(a,"-",b,"is",a-b)
    print(b,"-",a,"is",b-a)
    print(a,"-",a,"is",a-a)
    print("1-1/2 is",1-a)
    print("1/2-1 is",a-1)
    print("Testing division:")
    print(a,"/",b,"is",a/b)
    print(b,"/",a,"is",b/a)
    print(a,"/",a,"is",a/a)
    print("6/(1/2) is",6/a)
    print("(1/2)/6 is",a/6)
    n = Rat(7,1)
    print("Testing floor division:")
    print("Fraction",n,"//",2,"is",n//2)
    print(2,"//","Fraction",n,"is",2//n)
    # print(2,"//",a,"is",2//a) ERROR
    # print(a,"//",b,"is",a//b) ERROR
    print("Testing exponentiation:")
    print(b,"**",3,"is",b**3)
    print(b,"**",-3,"is",b**-3)
    print(b,"**",0,"is",b**0)
    print("Testing cast:")
    print("Casting",b,"to int: ",int(b))
    print("Casting",b,"to float: ",float(b))
    print("Testing problem:")
    print(a+1/2)
