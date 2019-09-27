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

# A fraction instance. Initialized as:
# 
#     Fraction(numerator=0, denominator=1)
# 
# The arguments to the Fraction must either be Python integers or
# convertible to Python float via the 'float' built-in function.
class Fraction:
    def __init__(self, numerator=0, denominator=1, reduce=True):
        # If a Fraction was provided, shortcut and exit.
        if (type(numerator) == Fraction):
            self._numerator   = numerator.numerator
            self._denominator = numerator.denominator
            return
        # If "None" was provided, then set this fraction to have "None" values.
        elif (type(None) in (type(numerator), type(denominator))):
            self._numerator = None
            self._denominator = None
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
                _, denominator = numerator.as_integer_ratio()
                numerator /= _
        # Check to make sure the "infinity" cases are handled consistently.
        if (numerator == 0):
            if (denominator != 0): denominator //= abs(denominator)
        elif (denominator == 0):   numerator   //= abs(numerator)
        # Transfer the negativity to the numerator if possible.
        if   (numerator < 0) and (denominator < 0):
            numerator, denominator = -numerator, -denominator
        elif (numerator >= 0) and (denominator < 0):
            numerator, denominator = -numerator, -denominator
        # Simplify the numerator and denomninator if appropriate.
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
    def __float__(self): return (self.numerator / self.denominator)
    # Convert this fraction to an integer.
    def __int__(self):   return (self.numerator // self.denominator)

    # Provide a represenation of this Fraction.
    def __repr__(self): return f"{self.__class__.__name__}({self.numerator}, {self.denominator})"
    # Provide a string representation of this Fraction.
    def __str__(self): return f"{self.numerator} / {self.denominator}"

    # Add two fractions.
    @float_fallback
    def __add__(a, b):
        an, bn = a.numerator, b.numerator
        ad, bd = a.denominator, b.denominator
        return Fraction( an * bd + bn * ad,  ad * bd )
    @float_fallback
    def __radd__(b, a): return a + b

    #  Subtract two fractions.
    @float_fallback
    def __sub__(a, b):
        an, bn = a.numerator, b.numerator
        ad, bd = a.denominator, b.denominator
        return Fraction( an * bd - bn * ad,  ad * bd )
    @float_fallback
    def __rsub__(a, b): return -a + b

    # Multiply two fractions.
    @float_fallback
    def __mul__(a, b):
        return Fraction(a.numerator * b.numerator, a.denominator * b.denominator)
    @float_fallback
    def __rmul__(a, b): return a * b

    # Divide two fractions.
    @float_fallback
    def __truediv__(a, b):
        return Fraction(a.numerator * b.denominator, a.denominator * b.numerator)
    @float_fallback
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
    @float_fallback
    def __mod__(a, b):
        mult = a // b
        return a - b * mult
    @float_fallback
    def __rmod__(b, a):
        return a % b

    # Compute a raised to the power b.
    # 
    # WARNING: Precision is lost if a root is taken (fractional power).
    #          The number is converted to a float then converted back.
    def __pow__(a, b):
        # Handle the integer power case.
        if (type(b) == int):
            if (b == 0):  return Fraction(1,1,reduce=False)
            elif (b > 0): return Fraction(a.numerator**b, a.denominator**b)
            else:         return Fraction(a.denominator**(abs(b)), a.numerator**(abs(b)))
        else:
            b = Fraction(b)
            # Handle the fractional power case (break up into two steps).
            intermediate = a**b.numerator
            return float(intermediate)**(1/b.denominator)
    def __rpow__(b, a):
        return a ** b

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
if __name__ == "__main__":
    print()
    print("1/3              : ",1/3)
    print("a = Fraction(1,3): ",Fraction(1,3))
    print("b = Fraction(1/3): ",Fraction(1/3))
    a = Fraction(1,3)
    b = Fraction(1/3)
    print()

    print("a - b:             ",a - b)
    print("float(a - b):      ",float(a - b))
    print("a_sqrt = a**(1/2): ",a**(1/2))
    a_sqrt = a**(1/2)
    print()

    print("error = a - a_sqrt**2: ",a - a_sqrt**2)
    error = a - a_sqrt**2
    print("float(error):          ",float(error))
    print()

    print("Fraction():       ",Fraction())
    print("bool(Fraction()): ",bool(Fraction()))
    print("bool(a):          ",bool(a))
    print()

    print("c = Fraction(float('inf')): ",Fraction(float('inf')))
    c = Fraction(float('inf'))
    print("Fraction(-float('inf')):    ",Fraction(-float('inf')))
    print("10 / c:  ",Fraction(10) / c)
    print("c / 10:  ",c / Fraction(10))
    print("c * 100: ",c * Fraction(100))
    print("-c:      ",-c)
    print()

    print("d = Fraction(1,float('inf')): ",Fraction(1,float('inf')))
    d = Fraction(1,float('inf'))
    print("Fraction(1,-float('inf')):    ",Fraction(1,-float('inf')))
    print("10 / d:  ",Fraction(10) / d)
    print("d / 10:  ",d / 10)
    print("d * 100: ",d * 100)
    print("-d:      ",-d)
