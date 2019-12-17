# Fraction class that provides infinite-precision rational numbers
# (making use of Python default arbitrarily large integers).

from math import gcd
class UsageError(Exception): pass

# This class-operator wrapper automatially returns a `float` type if
# the operand is neither the same type as the operator nor an `int`.
def float_fallback(operator):
    def wrapped_operator(primary, other):
        # Call the standard operator if the types match.
        if (type(other) == type(primary)):
            return operator(primary, other)
        # Handle integers in a special way (so that users can put ints
        # in code without breaking fractions).
        elif (type(other) == int):
            return operator(primary, type(primary)(other))
        # Otherwise, cast "other" to the type of "primary", call the
        # operator, then cast the result back to a 'float' type.
        else:
            return float(operator(primary, type(primary)(other)))
    return wrapped_operator

# This class-operator wrapper automatially returns a `float` type if
# the operand is neither the same type as the operator nor an `int`.
def auto_same_type(operator):
    def wrapped_operator(primary, other):
        # Call the standard operator if the types match.
        if (type(other) == type(primary)):
            return operator(primary, other)
        # Otherwise, cast "other" into the same type.
        else: return operator(primary, type(primary)(other))
    return wrapped_operator

# This can be 'float_fallback' or 'auto_same_type'.
TYPE_HANDLER = auto_same_type

# A fraction instance. Initialized as:
# 
#     Fraction(numerator=0, denominator=1)
# 
# The arguments to the Fraction must either be Python integers or
# convertible to Python float via the 'float' built-in function.
class Fraction:
    def __init__(self, numerator=0, denominator=1, reduce=True):
        # If "None" was provided, then set this fraction to have "None" values.
        if ((numerator is None) or (denominator is None)):
            self._numerator = None
            self._denominator = None
            return
        # If a Fraction was provided, shortcut and exit.
        elif (type(numerator) == Fraction):
            self._numerator   = numerator.numerator
            self._denominator = numerator.denominator
            return
        # Convert input arguments to the appropriate type.
        if (type(numerator) != int):
            try:    numerator = float(numerator)
            except: raise(UsageError(f"Unsupported numerator number type '{type(numerator)}', could not convert to float."))
        if (type(denominator) != int):
            try:    denominator = float(denominator)
            except: raise(UsageError(f"Unsupported denominator number type '{type(denominator)}', could not convert to float."))
        # Convert non-integer numerator and denominator appropriately.
        if (type(numerator) == float):
            # Handle 'infinity' in an expected way (without errors).
            if (abs(numerator) == float('inf')):
                if (abs(denominator) != float('inf')):
                    v = int(numerator < 0)
                    numerator, denominator = (-1)**(numerator < 0), 0
                elif (denominator < 0):
                    numerator, denominator = -numerator, -denominator
            # The normal usage case.
            else:
                numerator, _ = numerator.as_integer_ratio()
                denominator *= _
        if (type(denominator) == float):
            # Handle 'infinity' in an expected way (without errors).
            if (abs(denominator) == float('inf')):
                if (abs(numerator) != float('inf')):
                    numerator, denominator = 0, 1
            # The normal usage case.
            else:
                denominator, _ = denominator.as_integer_ratio()
                numerator *= _
        # Check to make sure the "infinity" cases are handled consistently.
        if (numerator == 0):
            if (denominator != 0): denominator //= abs(denominator)
        elif (denominator == 0):   numerator   //= abs(numerator)
        # Transfer the negativity to the numerator if possible.
        if   (numerator < 0) and (denominator < 0):
            numerator, denominator = -numerator, -denominator
        elif (numerator >= 0) and (denominator < 0):
            numerator, denominator = -numerator, -denominator
        # If this is indeterminate form, do not reduce!
        if (numerator == 0 == denominator): reduce = False
        # Simplify the numerator and denominator if appropriate.
        if reduce:
            divisor = gcd(numerator, denominator)
            numerator   //= divisor
            denominator //= divisor
        # Store the numerator and denominator of this fraction.
        self._numerator = numerator
        self._denominator = denominator

    # Provide access to the numerator and denominator of this Fraction.
    @property
    def numerator(self): return self._numerator
    @property
    def denominator(self): return self._denominator

    # Convert this fraction into a float.
    def __float__(self):
        if (self.denominator == 0):
            if   (self.numerator > 0): return float('inf')
            elif (self.numerator < 0): return -float('inf')
            else:                      return float(0)
        else: return (self.numerator / self.denominator)
    # Convert this fraction to an integer.
    def __int__(self):   return (self.numerator // self.denominator)

    # Provide a represenation of this Fraction.
    def __repr__(self): return f"{self.__class__.__name__}({self.numerator}, {self.denominator})"
    # Provide a string representation of this Fraction.
    def __str__(self): return f"{self.numerator} / {self.denominator}"

    # Add two fractions.
    @TYPE_HANDLER
    def __add__(a, b):
        an, bn = a.numerator, b.numerator
        ad, bd = a.denominator, b.denominator
        return Fraction( an * bd + bn * ad,  ad * bd )
    @TYPE_HANDLER
    def __radd__(b, a): return a + b

    #  Subtract two fractions.
    @TYPE_HANDLER
    def __sub__(a, b):
        an, bn = a.numerator, b.numerator
        ad, bd = a.denominator, b.denominator
        return Fraction( an * bd - bn * ad,  ad * bd )
    @TYPE_HANDLER
    def __rsub__(a, b): return -a + b

    # Multiply two fractions.
    @TYPE_HANDLER
    def __mul__(a, b):
        return Fraction(a.numerator * b.numerator, a.denominator * b.denominator)
    @TYPE_HANDLER
    def __rmul__(a, b): return a * b

    # Divide two fractions.
    @TYPE_HANDLER
    def __truediv__(a, b):
        return Fraction(a.numerator * b.denominator, a.denominator * b.numerator)
    @TYPE_HANDLER
    def __rtruediv__(b, a):
        return a / b

    # Divide two fractions (with integer result).
    def __floordiv__(a, b):
        b = Fraction(b)
        result = a / b
        result.numerator -= result.numerator % result.denominator
        return result.numerator // result.denominator
    def __rfloordiv__(b, a):
        a = Fraction(a)
        return a // b

    # Compute a mod b (and b mod a).
    @TYPE_HANDLER
    def __mod__(a, b):
        mult = a // b
        return a - b * mult
    @TYPE_HANDLER
    def __rmod__(b, a):
        return a % b

    # Compute a raised to the power b.
    # 
    # WARNING: Precision is lost if a root is taken (fractional power).
    def __pow__(a, b):
        # Handle the integer power case.
        if (type(b) == int):
            if (b == 0):  return Fraction(1,1,reduce=False)
            elif (b > 0): return Fraction(a.numerator**b, a.denominator**b)
            else:         return Fraction(a.denominator**(abs(b)), a.numerator**(abs(b)))
        else:
            # Handle the fractional power case (break up into two steps).
            b = Fraction(b)
            # Handle infinite values in a reasonable way.
            if (b.denominator == 0):
                if   (b.numerator > 0): return Fraction(a * float('inf'))
                elif (b.numerator < 0): return Fraction(0)
                else: b = Fraction(0,1)
            # If this is essentially an integer, just use that operation.
            elif (b.denominator == 1): return a ** b.numerator
            # First do the integer power (exact operation).
            intermediate = a ** abs(b.numerator)
            # Second do the inexact operation (converts to a float).
            result = float(intermediate) ** (1 / b.denominator)
            if (b.numerator < 0):
                if (result != 0): return 1 / result
                else:             return Fraction(1,0)
            return result

    # When taking a Fraction as a power..
    def __rpow__(b, a):
        return Fraction(a) ** b

    # Negate this Fraction.
    def __neg__(a): return Fraction(-a.numerator, a.denominator, reduce=False)

    # Make this Fraction positive.
    def __abs__(a): return Fraction(abs(a.numerator), a.denominator, reduce=False)

    # a equals b. (special case is False for comparison with None)
    def __eq__(a, b):
        b = Fraction(b)
        return (a.numerator == b.numerator) and (a.denominator == b.denominator)

    # a less than b.
    def __lt__(a, b):
        b = Fraction(b)
        return a.numerator * b.denominator < a.denominator * b.numerator

    # a greater than b
    def __gt__(a, b):
        b = Fraction(b)
        return a.numerator * b.denominator > a.denominator * b.numerator

    # a less than or equal to b
    def __le__(a, b):
        b = Fraction(b)
        return a.numerator * b.denominator <= a.denominator * b.numerator

    # a greater than or equal to b
    def __ge__(a, b):
        b = Fraction(b)
        return a.numerator * b.denominator >= a.denominator * b.numerator

    # a not equal to zero
    def __bool__(a):
        return (a.numerator != 0) and (a.denominator != 0)

    # Return the rounded version of this Fraction.
    def __round__(self, *args, **kwargs):
        return round(float(self), *args, **kwargs)


# Some minor testing and demonstration code. Anything not tested here
# has not been explicitly tested!! Test cases will be built as more
# features are used.
def _test_Fraction(display=False, print=lambda *args, **kwargs: None):
    # If display is desired, re-assign the real print function.
    if display:
        import builtins
        print = builtins.print

    print()
    print("1/3              : ",1/3)
    print("a = Fraction(1,3): ",Fraction(1,3))
    print("b = Fraction(1/3): ",Fraction(1/3))

    a = Fraction(1,3)
    assert(a.numerator == 1)
    assert(a.denominator == 3)
    b = Fraction(1/3)
    assert(a != b)

    print()
    print("a - b:             ",a - b)
    assert((a - b) == Fraction(1, 54043195528445952))
    print("float(a - b):      ",float(a - b))
    print("a_sqrt = a**(1/2): ",a**(1/2))
    a_sqrt = a**(1/2)
    print()

    print("error = a - a_sqrt**2: ",a - a_sqrt**2)
    error = a - a_sqrt**2
    print("float(error):          ",float(error))
    print("error - float(a-b):    ",float(error - float(a-b)))
    # Error should certainly be less than SQRT error for rounded fractions.
    assert(error <= float(a-b)**(1/2))  
    print()

    print("Fraction():       ",Fraction())
    assert(Fraction() == Fraction(0,1))
    print("bool(Fraction()): ",bool(Fraction()))
    assert(not bool(Fraction()))
    print("bool(a):          ",bool(a))
    assert(bool(a))
    print()

    print("c = Fraction(float('inf')): ",Fraction(float('inf')))
    c = Fraction(float('inf'))
    assert(c == Fraction(1,0))
    print("Fraction(-float('inf')):    ",Fraction(-float('inf')))
    assert(Fraction(-float('inf')) == Fraction(-1,0))
    print("10 / c:    ",Fraction(10) / c)
    assert((Fraction(10) / c) == Fraction(0,1.))
    print("c / 10:    ",c / Fraction(10))
    assert((c / Fraction(10)) == Fraction(1.,0))
    print("c * 100.0: ",c * 100.0)
    assert((c * 100.0) == float('inf'))
    print("-c * 10:   ",-c)
    assert((-c*10) == Fraction(-1,0))
    print("c - c:     ",c - c)
    assert((c-c) == Fraction(0,0))
    print("float(c-c) ", float(c-c))
    assert(float(c-c) == 0.0)
    print()
    print("a * (c - c)", a*(c - c))
    assert((a*(c-c)) == Fraction(0,0))
    print("a / (c - c)", a/(c - c))
    assert((a/(c-c)) == Fraction(0,0))
    print("a ** c     ", a**c)
    assert((a**c) == Fraction(1,0))
    print("a ** (-c)  ", a**(-c))
    assert((a**(-c)) == Fraction(0,1))
    print("(-a) ** c  ", (-a)**c)
    assert(((-a)**c) == Fraction(-1,0))
    print()

    print("d = Fraction(1,float('inf')): ",Fraction(1,float('inf')))
    d = Fraction(1,float('inf'))
    assert(d == Fraction(0,1))
    print("Fraction(1,-float('inf')):    ",Fraction(1,-float('inf')))
    assert(Fraction(1,-float('inf')) == Fraction(0,1))
    print("10 / d:  ",10 / d)
    assert((10/d) == Fraction(1,0))
    print("d / 10:  ",d / 10)
    assert((d/10) == Fraction(0,1))
    print("d * 100: ",d * 100)
    assert((d * 100) == Fraction(0,1))
    print("c * d:   ",c * d)
    assert((c * d) == Fraction(0,0))
    print("-d:      ",-d)
    assert((-d) == Fraction(0,1))
    print()

    # Make sure that the power operator behaves in expected ways.
    a = Fraction(7,3)
    assert(a**2 == a**2.)
    assert(a**(-2) == a**(-2.))
    assert(Fraction(2.0, -1/2) == Fraction(-4))


if __name__ == "__main__":
    print()
    print("Testing Fraction..")
    _test_Fraction(display=False)
    print("done.")
    
