# This file provides utilities for constructing polynomial interpolants.
# 
# The following objects are provided:
# 
#   Spline           -- A piecewise polynomial interpolant with
#                       evaluation, derivative, integration, negation,
#                       multiplication, addition, and a string method.
#   Polynomial       -- A polynomial with nonrepeating coefficients,
#                       evaluation, derivative, and a string method.
#   NewtonPolynomial -- A Newton polynomial with stored coefficients,
#                       offsets, evaluation, derivative, and a string method.
#   Vector -- A limited-usage pure Python vector class, used for 
#             convenience, but WARNING it does check for usage errors.
# 
# The following functions are provided:
# 
#   fit              -- Use a local polynomial interpolant to estimate
#                       derivatives and construct an interpolating
#                       Spline with a specified level of continuity.
#   polynomial       -- Given x and y values, this produces a minimum
#                       degree Newton polynomial that interpolates the
#                       provided points (this is a *single* polynomial).
#   local_polynomial -- Construct a local polynomial interpolant about
#                       a given point in a list of points.
#   local_quadratic  -- Given a set of points and an index, construct
#                       a weighted quadratic fit about the point at
#                       "index" in the set. Use Shepard weighting.
#   local_facet      -- Given a set of opints and an index, construct
#                       local linear fits near the index and return
#                       the local linear regression with lowest error.
#   inverse          -- Given a function, use the Newton method to
#                       find the nearest inverse to a guess point.
# 


# This general purpose exception will be raised during user errors.
class UsageError(Exception): pass

# This class-method wrapper function ensures that the class method
# recieves a fraction and returns a "float" if a fraction was not
# provided as input (ensuring internal methods only recieve fractions).
def float_fallback(class_method):
    from util.math.fraction import Fraction
    def wrapped_method(obj, x):
        # Check for vector usage.
        try: return [wrapped_method(obj, v) for v in x]
        except:
            # Return perfect precision if that was provided.
            if (type(x) == Fraction): return class_method(obj, x)
            # Otherwise, use exact arithmetic internally and return as float.
            else: return float(class_method(obj, Fraction(x)))
    return wrapped_method

# A piecewise polynomial function that supports evaluation,
# differentiation, and stitches together an arbitrary sequence of
# function values (and any number of derivatives) at data points
# (knots). Provides evaluation, derivation, and string methods.
# It is recommended to use EXACT ARITHMETIC internally (which will
# automatically return floats unless Fraction objects are provided as
# evaluation points). Exact arithmetic is achieved by providing knots
# and values composed of Fraction objects.
# 
# Spline(knots, values):
#   Given a sequence of "knots" [float, ...] and an equal-length
#   sequence of "values" [[float, ...], ...] that define the
#   function value and any number of derivatives at every knot,
#   construct the piecewise polynomial that is this spline.
class Spline:
    # Define private internal variables for holding the knots, values,
    # and the functions. Provide access to "knots" and "values" properties.
    _knots = None
    _values = None
    _functions = None
    _derivative = 0
    @property
    def knots(self): return self._knots
    @knots.setter
    def knots(self, knots): self._knots = list(knots)
    @property
    def values(self): return self._values
    @values.setter
    def values(self, values): self._values = [[v for v in vals] for vals in values]
    @property
    def functions(self): return self._functions
    @functions.setter
    def functions(self, functions): self._functions = list(functions)
    
    # Create a piecewise polynomial that matches the given values at
    # the given knots, alternatively provide a list of functions that
    # exist over a set of intervals defined by knots.
    def __init__(self, knots, values=None, functions=None):
        assert(len(knots) >= 1)
        self.knots = knots
        # Store the 'values' at each knot if they were provided.
        if (values is not None):
            assert(len(knots) == len(values))
            self.values = values
        # Store the 'functions' over each interval if they were provided.
        if (functions is not None):
            # TODO: Verify that the provided functions match the values.
            assert(len(functions) == len(knots)-1)
            self.functions = functions
        # Use the 'values' to generate 'functions' if no functions were provided.
        elif (values is not None):
            # Create the polynomial functions for all the pieces.
            self._functions = []
            for i in range(len(knots)-1):
                v0, v1 = self.values[i], self.values[i+1]
                k0, k1 = self.knots[i], self.knots[i+1]
                x = [k0, k1]
                y = [v0[0], v1[0]]
                derivs = {i*'d'+'x0':v0[i] for i in range(1,len(v0))}
                derivs.update({i*'d'+'x1':v1[i] for i in range(1,len(v1))})
                # Construct a local polynomial that interpolates the 
                # function (and derivative) values.
                f = polynomial(x, y, **derivs)
                # Store the function (assuming it's correct).
                self._functions.append( f )
        # No 'values' nor 'functions' were provided, usage error.
        else: raise(UsageError("Either 'values' or 'functions' must be provided with knots to construct a spline."))

    # Evaluate this Spline at a given x coordinate.
    def function_at(self, x):
        # If "x" was given as an iterable, then iterate over it.
        try:    return [self.function_at(v) for v in x]
        except: pass
        # Deterimine which interval this "x" values lives in.
        if (x <= self.knots[0]):    return self._functions[0]
        elif (x >= self.knots[-1]): return self._functions[-1]
        # Find the interval in which "x" exists.
        for i in range(len(self.knots)-1):
            if (self.knots[i] <= x < self.knots[i+1]): break
        # If no interval was found, then something unexpected must have happened.
        else:
            class UnexpectedError(Exception): pass
            raise(UnexpectedError("This problem exhibited unexpected behavior. Check code."))
        # Return the applicable function.
        return self._functions[i]

    # Compute the integral of this Spline.
    def integral(self, i=1): return self.derivative(-i)

    # Compute the first derivative of this Spline.
    def derivative(self, d=1):
        # For integration, adjust the additive term to reflect the
        # expected lefthand side value of each function.
        if (d < 0):
            deriv_funcs = self._functions
            for i in range(-d):
                deriv_funcs = [f.integral(1) for f in deriv_funcs]
                total = self(self.knots[0])
                for i in range(len(deriv_funcs)):
                    deriv_funcs[i].coefficients[-1] = (
                        total - deriv_funcs[i](self.knots[i]))
                    total = deriv_funcs[i](self.knots[i+1])
        else:
            # Create a spline with the same knot sequence and all the
            # derivative functions and associated values.
            deriv_funcs = self._functions
            for i in range(d):
                deriv_funcs = [f.derivative(1) for f in deriv_funcs]

        # Construct the new spline, pass the "values" even though
        # nothing will be done with them. Assign the new "functions".
        s = Spline(self.knots, self.values, functions=deriv_funcs)
        s._derivative = self._derivative + d
        # Return the new derivative Spline.
        return s

    # Evaluate this Spline at a given x coordinate.
    def __call__(self, x):
        # If "x" was given as an iterable, then iterate over it.
        try:    return [self(v) for v in x]
        except: pass
        # Get the appropriate function and compute the output value.
        return self.function_at(x)(x)
        
    # Add this "Spline" object to another "Spline" object.
    def __add__(self, other):
        # Check for correct usage.
        if (type(other) != type(self)):
            raise(UsageError(f"Only '{type(self)} objects can be added to '{type(self)}' objects, but '{type(other)}' was given."))
        # Generate the new set of knots.
        knots = sorted(set(self._knots + other._knots))
        # Compute the functions over each interval.
        functions = []
        for i in range(len(knots)-1):
            # Get the knot, nearby knots, and order of resulting
            # polynomial at this particular knot.
            left, right = knots[i], knots[i+1]
            k = knots[i]
            my_poly = self.function_at(k)
            other_poly = other.function_at(k)
            order = max(len(my_poly.coefficients), len(other_poly.coefficients))
            # Evaluate the function at equally spaced "x" values, TODO:
            # this should be Chebyshev nodes for numerical stability.
            x = [(step / (order-1)) * (right - left) + left for step in range(order)]
            y = [self(node) + other(node) for node in x]
            # Construct the interpolating polynomial.
            functions.append( polynomial(x, y) )
            f = functions[-1]
        # Return the new added function.
        return Spline(knots, functions=functions)

    # Multiply this "Spline" obect by another "Spline".
    def __mul__(self, other):
        # Check for correct usage.
        if (type(other) != type(self)):
            raise(UsageError(f"Only '{type(self)} objects can be multiplied by '{type(self)}' objects, but '{type(other)}' was given."))
        # Generate the new set of knots.
        knots = sorted(set(self._knots + other._knots))
        # Compute the functions over each interval.
        functions = []
        for i in range(len(knots)-1):
            # Get the knot, nearby knots, and order of resulting
            # polynomial at this particular knot.
            left, right = knots[i], knots[i+1]
            k = knots[i]
            my_poly = self.function_at(k)
            other_poly = other.function_at(k)
            order = max(2,len(my_poly.coefficients) + len(other_poly.coefficients) - 1)
            # Evaluate the function at equally spaced "x" values, TODO:
            # this should be Chebyshev nodes for numerical stability.
            x = [(step / (order-1)) * (right - left) + left for step in range(order)]
            y = [self(node) * other(node) for node in x]
            # Construct the interpolating polynomial.
            functions.append( polynomial(x, y) )
            f = functions[-1]
        # Return the new added function.
        return Spline(knots, functions=functions)

    # Raise an existing spline to a power.
    def __pow__(self, number):
        if (type(number) != int) or (number < 1):
            raise(TypeError(f"Only possible to raise '{type(self)}' to integer powers greater than  or equal to 1."))
        # Start with a copy of "self", multiply in until complete.
        outcome = Spline(self.knots, values=self.values,
                         functions=[Polynomial(f) for f in self.functions])
        for i in range(number-1): outcome = outcome * self
        return outcome

    # Subtract another spline from this spline (add after negation).
    def __sub__(self, other): return self + (-other)

    # Negate this spline, create a new one that is it's negative.
    def __neg__(self):
        # Create a new spline, but negate all internal function coefficients.
        return Spline(self.knots, functions=[-f for f in self.functions])

    # Produce a string description of this spline.
    def __str__(self):
        s = "Spline:\n"
        s += f" [-inf, {self.knots[1]}]  =  "
        s += str(self._functions[0]) + "\n"
        for i in range(1,len(self.knots)-2):
            s += f" ({self.knots[i]}, {self.knots[i+1]}]  =  "
            s += str(self._functions[i]) + "\n"
        s += f" ({self.knots[-2]}, inf)  =  "
        s += str(self._functions[-1])
        return s

# A generic Polynomial class that stores coefficients in a fully
# expanded form. Provides numerically stable evaluation, derivative,
# integration, and string operations for convenience. Coefficients
# should go highest to lowest order.
# 
# Polynomial(coefficients):
#    Given coefficients (or optionally a "NewtonPolynomial")
#    initialize this Polynomial representation of a polynomial function.
# 
# EXAMPLE:
# Polynomial([a,b,c])
#   = a x^2  +  b x  +  c
class Polynomial:
    # Initialize internal storage for this Polynomial.
    _coefficients = None
    # Protect the "coefficients" of this class with a getter and
    # setter to ensure a user does not break them on accident.
    @property
    def coefficients(self): return self._coefficients
    @coefficients.setter
    def coefficients(self, coefs): self._coefficients = list(coefs)
    # Initialize this Polynomial given coefficients.
    def __init__(self, coefficients=(0,)):
        # If the user initialized this Polynomial with a Newton
        # Polynomial, then extract the points and coefficients.
        if (type(coefficients) == NewtonPolynomial):
            # Unpack the NewtonPolynomial into coefficients and points
            # by multiplying everything out (and unrolling it).
            newton_polynomial = coefficients
            c, p = newton_polynomial.coefficients, newton_polynomial.zeros
            coefficients = [c[0]]
            for i in range(1,len(c)):
                # Compute the old coefficients multiplied by a constant and
                # add the lower power coefficients that are shifted up.
                coefficients.append(c[i])
                coefficients = [coefficients[0]] + \
                               [coefficients[j+1]-p[i]*coefficients[j]
                                for j in range(len(coefficients)-1)]
        # If the user initialized this Polynomial, do a copy.
        elif (type(coefficients) == type(self)):
            coefficients = coefficients.coefficients
        # Store the coeficients.
        self.coefficients = coefficients # <- implicit copy
        # Remove all leading 0 coefficients.
        for i in range(len(self._coefficients)):
            if (self._coefficients[i] != 0): break
        else: i = len(self._coefficients)-1 # <- include the last 0, if needed
        self._coefficients = self._coefficients[i:]

    # Evaluate this Polynomial at a point "x" in a numerically stable way.
    def __call__(self, x):
        if (len(self._coefficients) == 0): return 0
        total = self._coefficients[0]
        for d in range(1,len(self._coefficients)):
            total = self._coefficients[d] + x * total
        return total

    # Construct the polynomial that is the integral of this polynomial.
    def integral(self, i=1): return self.derivative(-i)

    # Construct the polynomial that is the derivative (or integral) of this polynomial.
    def derivative(self, d=1):
        if (d == 0): return self
        elif (d > 1):  return self.derivative(1).derivative(d-1)
        elif (d == 1): return Polynomial([c*i for (c,i) in zip(
                self._coefficients, range(len(self._coefficients)-1,0,-1))])
        elif (d < -1):  return self.derivative(-1).derivative(d+1)
        elif (d == -1): return Polynomial([c/(i+1) for (c,i) in zip(
                self._coefficients, range(len(self._coefficients)-1,-1,-1))]+[0])

    # Determines if two polynomials are equivalent.
    def __eq__(self, other):
        if type(other) == NewtonPolynomial: other = Polynomial(other)
        return all(c1 == c2 for (c1, c2) in zip(self._coefficients, other._coefficients))

    # Negate this polynomial.
    def __neg__(self): return Polynomial([-c for c in self._coefficients])

    # Construct a string representation of this Polynomial.
    def __str__(self):
        s = ""
        for i in range(len(self._coefficients)):
            if (self._coefficients[i] == 0): continue
            if   (i == len(self._coefficients)-1): x = ""
            elif (i == len(self._coefficients)-2): x = "x"
            else:   x = f"x^{len(self._coefficients)-1-i}"
            s += f"{self._coefficients[i]} {x}  +  "
        # Remove the trailing 
        s = s.rstrip(" +")
        # Return the final string.
        return s

# Extend the standard Polymomial class to hold Newton polynomials with
# zeros in addition to the coefficients. This is more convenient when
# constructing interpolating polynomials from divided difference tables.
# 
# NewtonPolynomial(coefficients, zeros):
#    Given a set of coefficients and a set of zeros (offsets), of
#    the same length, construct a standard Newton
#    Polynomial. Coefficients are stored from highest order term to
#    lowest order. Earlier zeros are evaluated earlier.
# 
# EXAMPLE:
# NewtonPolynomial([a,b,c], [s1, s2, s3])
#   =  c + (x - s3)(b + (x - s2)(a))
# 
# NOTE:
#   Notice that "s1" is never used in the computation, as it is
#   redundant with the term "a".
class NewtonPolynomial(Polynomial):
    _zeros = None
    @property
    def zeros(self): return self._zeros
    @zeros.setter
    def zeros(self, zeros): self._zeros = list(zeros)

    # Store the coefficients and zeros for this Newton Polynomial.
    def __init__(self, coefficients, zeros):
        self.coefficients = coefficients
        # Strip off the 0-valued coefficients (all these will be 0'd out).
        for i in range(len(self._coefficients)):
            if (self._coefficients[i] != 0): break
        else: i = len(self._coefficients)-1
        # Store locally.
        self._coefficients = self._coefficients[i:]
        self._zeros = [p for (j,p) in enumerate(zeros) if j >= i]

    # Construct the polynomial that is the derivative of this
    # polynomial by converting to polynomial form and differntiating.
    def derivative(self, d=1): return Polynomial(self).derivative(d)

    # Evaluate this Newton Polynomial (in a numerically stable way).
    def __call__(self, x):
        total = self._coefficients[0]
        for d in range(1,len(self._coefficients)):
            total = self._coefficients[d] + (x - self._zeros[d]) * total
        return total

    # Negate this polynomial.
    def __neg__(self): return -Polynomial(self)

    # Construct a string representation of this Newton Polynomial.
    def __str__(self):
        s = f"{self._coefficients[0]}"
        for i in range(1,len(self._coefficients)):
            sign = "-" if (self._zeros[i] >= 0) else "+"
            s = f"{self._coefficients[i]} + (x {sign} {abs(self._zeros[i])})({s})"
        return s


# A "Vector" subclass of a "list" for convenience, allows component-
# wise addition, subtraction, multiplication, and (float) division.
# Also provides negation, "dot" product, and "norm" methods (2-norm).
# 
# WARNING: This not not an extensively supported vector class, it only
#          handles a few expected simple operations. It does not
#          check for usage errors and may provided unexpected results
#          if used incorrectly.
class Vector(list):
    # Define add, subtract, multiply, divide, assume a vector is given
    # and if the object given is not iterable assume it's scalar.
    def __add__(self, value):
        try:              return Vector(v1 + v2 for (v1,v2) in zip(self,value))
        except TypeError: return Vector(v + value for v in self)
    def __sub__(self, value):
        try:              return Vector(v1 - v2 for (v1,v2) in zip(self,value))
        except TypeError: return Vector(v - value for v in self)
    def __mul__(self, value):
        try:              return Vector(v1 * v2 for (v1,v2) in zip(self,value))
        except TypeError: return Vector(v * value for v in self)
    def __truediv__(self, value):
        try:              return Vector(v1 / v2 for (v1,v2) in zip(self,value))
        except TypeError: return Vector(v / value for v in self)
    # Define "absolute value", "dot product" and "norm".
    def __abs__(self):  return Vector(map(abs,self))
    def dot(self, vec): return sum(v1*v2 for (v1,v2) in zip(self, vec))
    def norm(self):     return self.dot(self)**(1/2)


# Given data points "x" and data values "y", construct an
# interpolating spline over the given points with at least the
# specified level of continuity using a sufficiently continuous
# polynomial fit over neighboring points.
#  
# x: A strictly increasing sequences of numbers.
# y: Function values associated with each point.
# 
# continuity:
#   The level of continuity desired in the interpolating function.
# 
def fit(x, y, continuity=0):
    assert( len(x) == len(y) )
    # Sort the "x" values if they were not given in sorted order.
    if not all(x[i] < x[i+1] for i in range(len(x)-1)):
        indices = sorted(range(len(x)), key=lambda i: x[i])
        x = [x[i] for i in indices]
        y = [y[i] for i in indices]
    # Get the knots and initial values for the spline.
    knots = [v for v in x]
    values = [[v] for v in y]
    order = continuity + 1
    for i in range(len(x)):
        lp = local_polynomial(x,y,i,order=order)
        for d in range(1,continuity+1):
            dlp = lp.derivative()
            values[i].append( dlp(x[i]) )
    # Return the interpolating spline.
    return Spline(knots, values)


# Given unique "x" values and associated "y" values (of same length),
# construct an interpolating polynomial with the Newton divided
# difference method. Allows for the incorporation of derivative
# constraints via keyword arguments. I.e. the keyword argument "dx0=0"
# makes the derivative of the polynomial 0 at the first element of x.
# Returns a NewtonPolynomial object.
# 
# EXAMPLE:
# polynomial([1, 2, 3], [-1, 2, -3], dx1=0)
#   => cubic polynomial interpolating the points (0,-1), (1,2), and
#      (2,-3) that has a first derivative of 0 at x=2.
def polynomial(given_x, given_y, **derivs):
    # Sort the data by "x" value.
    indices = sorted(range(len(given_x)), key=lambda i: given_x[i])
    # Construct the initial x and y (with repetitions to match derivatives).
    x = []; y = []; index_ranges = {}
    for i in indices:
        dxs = sorted(key for key in derivs if (i == int(key.split('x')[1])))
        largest_d = ([key.count("d") for key in dxs] + [0])[0]
        index_ranges[(len(x), len(x)+largest_d)] = i
        x += [given_x[i]] * (largest_d+1)
        y += [given_y[i]] * (largest_d+1)
    # Compute the divided difference table.
    divisor = 1
    dd_values = [y]
    for d in range(1, len(x)):
        divisor *= d
        slopes = []
        for i in range(len(dd_values[-1])-1):
            # Identify which "original point index" this spot belongs to.
            idx = ''
            for (start, end) in index_ranges:
                if (start <= i <= end-d): idx = index_ranges[(start,end)]
            # Substitute in derivative value if it was provided.
            key = d*"d"+f"x{idx}"
            if key in derivs: dd = derivs[key] / divisor
            # Otherwise a value wasn't provided, compute the divided difference.
            elif (x[i+d] == x[i]): dd = 0
            else: dd = (dd_values[-1][i+1] - dd_values[-1][i]) / (x[i+d] - x[i])
            slopes.append( dd )
        # Add in the finished next row of the divided difference table.
        dd_values.append( slopes )
    # Get the divided difference (polynomial coefficients) in reverse
    # order so that the most nested value (highest order) is first.
    # Return as an interpolating polynomial in Newton form.
    return NewtonPolynomial(
        (row[0] for row in reversed(dd_values)),  reversed(x))


# Given a function, use the Newton method to find the nearest inverse
# to a guess. If this method gets trapped, it will exit without having
# computed an inverse.
def inverse(f, y, x=0, accuracy=2**(-52), max_steps=100):
    import fraction # To access the "ACCURACY" internal.
    from util.math.fraction import Fraction    
    # Change the accuracy to only be enough to achieve an acceptable answer.
    previous_acc = fraction.ACCURACY
    fraction.ACCURACY = 8 * int(1 / accuracy)
    # Get the derivative of the function.
    df = f.derivative()
    # Convert 'x' and 'y' to Fractions.
    x, y = Fraction(x), Fraction(y)
    # Search for a solution iteratively.
    for i in range(max_steps):
        diff = f(x) - y
        x -= diff / df(x)
        # Stop the loop if we have gotten close to the correct answer.
        if (abs(diff) <= accuracy): break
    # Check for correctness, warn the user if result is bad.
    if (abs(f(x) - y) > accuracy):
        import warnings
        warnings.warn(f"\n\n  The calculated inverse has high error ({float(abs(f(x)-y)):.2e}).\n"+
                      "  Consider providing better initial position.\n"+
                      "  This problem may also not be solvable.\n")
    # Reset the accuracy of the instantiated Fraction classes.
    fraction.ACCURACY = previous_acc
    # Return the final value.
    return Fraction(x)


# Construct a local polynomial interpolant over nearby data points.
# Given "x" values, "y" values, and the index "i" at which a local
# polynomial should be (roughly) centered.
def local_polynomial(x, y, i, order=3, local_indices=None,
                     local_derivs=None, **derivs):
    # Sort indices first by their index nearness to x[i], then
    # secondly by their value nearness to x[i].
    all_indices = list(range(len(x)))
    all_indices.sort(key=lambda j: ( abs(j - i), abs(x[j] - x[i])*abs(y[j]-y[i]) ))
    # Convert all provided derivative information to the relative
    # index in this local polynomial. Only incorporate as much
    # information as the given "order" polynomial can contain.
    if local_indices is None: local_indices = list()
    if local_derivs is None: local_derivs = dict()
    for step,j in enumerate(all_indices):
        if (order <= 0): break
        if j not in local_indices:
            local_indices.append(j)
            order -= 1
        # Append the next point if it is the same distance away,
        # before considering the derivatives of the current point.
        if ((order > 0) and (step+1 < len(all_indices))
            and (abs(x[j]-x[i]) == abs(x[all_indices[step+1]] - x[i]))):
            local_indices.append(all_indices[step+1])
            order -= 1
        # Find all the derivatives that were given for this index and
        # are capable of being approximated, highest derivative first.
        given_derivatives = sorted(key for key in derivs
                                   if int(key.split('x')[1]) == j)
        kept_derivatives = [key for key in given_derivatives
                            if order >= key.count('d')]
        # If derivatives are being approximated, then lower the 
        # remaining approximation power by the derivative number.
        if (len(kept_derivatives) > 0):
            order -= kept_derivatives[0].count('d')
        # Store the transformed derivative assignments.
        for key in kept_derivatives:
            # Update the index of the keys that are being kept.
            d, k = key.split('x')
            local_derivs['x'.join((d,str(step)))] = derivs[key]
    # Return the interpolating polynomial of specified order.
    return polynomial([x[j] for j in local_indices],
                      [y[j] for j in local_indices],
                      **local_derivs)


# Given x (points, list of floats) and y (values list of floats), and
# an idx (integer), build a weighted least squares quadratic fit over
# the point values using inverse squared distances from x[idx] as
# weights.
def local_quadratic(x, y, idx, **derivs):
    # Check for proper usage.
    assert(len(x) == len(y))
    assert(0 <= idx < len(x))
    # Convert "x" and "y" to vectors.
    x, y = Vector(x), Vector(y)
    # Shift the "x" and "y" to interpolate the point at "idx", remove
    # that point from both vectors.
    shift_x, shift_y = x[idx], y[idx]
    x.pop(idx); y.pop(idx)
    # Pop out any repetitions of the center "x" point (they don't affect the fit).
    while (shift_x in x):
        i = x.index(shift_x)
        x.pop(i)
        y.pop(i)
    # Construct the shifted point set (to interpolate the origin).
    sx, sy = x - shift_x, y - shift_y
    # Compute x and y normalized by their weight.
    nx, ny = sx/abs(sx), sy/abs(sx)
    # Compute the linear term in the quadratic fit (interpolating
    # the means of normalized data (on left and right).
    if (max(nx) > 0):
        pos_y = [yv for (xv,yv) in zip(nx,ny) if xv > 0]
        right = (1, sum(pos_y) / len(pos_y))
    else: right = (0, 0)
    if (min(nx) < 0):
        neg_y = [yv for (xv,yv) in zip(nx,ny) if xv < 0]
        left = (-1, sum(neg_y) / len(neg_y))
    else: left = (0, 0)
    # Compute the linear term.
    if (left[0] == right[0]): b = 0
    else: b = (right[1] - left[1]) / (right[0] - left[0])
    # Overwrite the linear term if it was provided.
    b = derivs.get(f"dx{idx}", b)
    # Compute the quadratic term in the fit, considering we want to
    # minimize:
    # 
    #   sum(( a*abs(sx) + [b*sx/abs(sx) - sy/abs(sx)] )**2)
    # 
    # Compute "a" by finding the length of "abs(sx)" that places it
    # closest to the point "ny - b*nx". This can be achieved by
    # projecting "ny - b*nx" onto the normalized line "abs(sx)". Then
    # the solution is (projection length / current length).
    a = (ny - nx*b).dot(abs(sx)) / sx.norm()**2
    # Overwrite the quadratic term if it was provided.
    a = derivs.get(f"ddx{idx}", a)
    # Return the quadratic in NewtonPolynomial form.
    return NewtonPolynomial([a,b,shift_y],[None,shift_x,shift_x])


# Construct local linear fits, pick the one that has the lowest
# residual and return it in a Polynomial form.
def local_facet(x, y, idx, size=3, local_indices=None, **derivs):
    # Check for correct usage.
    assert(len(x) >= size)
    assert(0 <= idx < len(x))
    assert(len(x) == len(y))
    # Initialize the "used indices".
    if (local_indices is None): local_indices = []
    # Look for this derivative being pre-defined, if so return the required function.
    key = f"dx{idx}"
    if key in derivs:
        local_indices += [idx]
        return polynomial([x[idx]], [y[idx]], dx0=derivs[key])
    # Adjust size to be a reasonable value based on available points,
    # should always be an odd number, should be at least 3.
    size = min(max(3, 2*(size//2)+1), len(x))
    # Construct local liner models about each group of points. This
    # function returns the (slope, intercept) pair for the given points.
    def fit_line(x, y):
        n = len(x)
        # Compute the denominator.
        x_sq_sum = sum(x[i]**2 for i in range(n))
        denominator = (n-1) * x_sq_sum - 2 * sum(x[i] * sum(x[i+1:]) for i in range(n))
        # Compute the numerator for slope term.
        y_sum = sum(y)
        x_dot_y = sum(x[i] * y[i] for i in range(n))
        slope_numerator = (n-1) * x_dot_y - sum(x[i] * (y_sum - y[i]) for i in range(n))
        # Compute the numerator for the intercept term.
        x_sum = sum(x)
        inter_numerator = sum(y[i]*(x_sq_sum - x[i]**2 - x[i]*(x_sum - x[i])) for i in range(n))
        # Return the slope and the intercept.
        return slope_numerator / denominator, inter_numerator / denominator
    # Store the ranges, functions, and residuals for the local linear fits.
    ranges = []
    functions = []
    residuals = []
    for center in range(idx-size//2, idx+size//2 +1):
        # Skip "center" indices outside of feasible values.
        if (center < 0) or (center >= len(x)): continue
        lower = max(0,center-size//2)
        upper = min(len(x), center+size//2+1)
        # Skip intervals with not enough points.
        if (upper - lower) < size: continue
        # Get the local points and values (for the local linear fit).
        ranges.append( list(range(lower,upper)) )
        local_x = x[lower:upper]
        local_y = y[lower:upper]
        # The residuals for 0 points are infinity, so these will not be chosen.
        if (len(local_x) == 0):
            residuals.append( float('inf') )
            functions.append( Polynomial() )
            continue
        # Fit a linear model.
        slope, intercept = fit_line(local_x,local_y)
        function = Polynomial((slope, intercept))
        # Store the function and its "error" (squared sum of errors).
        functions.append( function )
        # residual = sum((slope*lx + intercept - ly)**2
        #                for (lx,ly) in zip(local_x, local_y))
        # residuals.append( residual )
        residuals.append( abs(slope*x[idx] + intercept - y[idx]) )

    # Return the linear function that has the best local fit.
    center = residuals.index(min(residuals))
    local_indices += ranges[center]
    return functions[center]


# --------------------------------------------------------------------
#                            TESTING CODE

def _test_local_quadratic():
    # Minimize the error function to verify correctness of solution.
    def opt_min_fit(x,y,idx):
        # Numpy functions.
        from numpy import asarray, sum, dot
        from numpy.linalg import norm
        from util.optimize import minimize

        # Convert "x" and "y" to vectors.
        x, y = list(x), list(y)
        # Shift the "x" and "y" to interpolate the point at "idx", remove
        # that point from both vectors.
        shift_x, shift_y = float(x[idx]), float(y[idx])
        x.pop(idx); y.pop(idx)
        # Pop out any repetitions of the center "x" point (they don't affect the fit).
        while (shift_x in x):
            i = x.index(shift_x);  x.pop(i); y.pop(i)
        # Convert to numpy arrays.
        x, y = asarray(x, dtype=float), asarray(y, dtype=float)
        # Construct the shifted point set (to interpolate the origin).
        sx, sy = x - shift_x, y - shift_y
        # Compute x and y normalized by their weight.
        nx, ny = sx/abs(sx), sy/abs(sx)
        # Compute the linear term in the quadratic fit (interpolating
        # the means of normalized data (on left and right).
        if (max(nx) > 0):
            pos_y = [yv for (xv,yv) in zip(nx,ny) if xv > 0]
            right = (1, sum(pos_y) / len(pos_y))
        else: right = (0, 0)
        if (min(nx) < 0):
            neg_y = [yv for (xv,yv) in zip(nx,ny) if xv < 0]
            left = (-1, sum(neg_y) / len(neg_y))
        else: left = (0, 0)
        # Compute the linear term and quadratic term estimates.
        if (left[0] == right[0]): b = 0
        else: b = (right[1] - left[1]) / (right[0] - left[0])
        a = (ny - nx*b).dot(abs(sx)) / norm(sx)**2
        # Minimize this error function.
        w = 1 / sx**2
        error = lambda ab: sum(w*(a*sx**2 + b*sx - sy)**2)
        # Perform optimization.
        ab = minimize(error, [a,b], max_steps=100, display=False)
        # Check for greater than a "epsilon" error.
        if (error([a,b]) - error(ab)) > 2**(-52):
            print("ERROR: Quadratic fit was not as good as expected.")
            print("  error([a,b]): ",error([a,b]))
            print("  error(ab):    ",error(ab))
            exit()
        elif (error([a,b]) - error(ab)) > 0:
            print(f"WARNING: Deterministic quadratic fit was off by {error([a,b]) - error(ab)}.")
            print("  error([a,b]): ",error([a,b]))
            print("  error(ab):    ",error(ab))
        # Return optimized interpolating polynomial.
        return NewtonPolynomial([a,b,shift_y],[None,shift_x,shift_x])

    # Generate random "x", "y", and "i" to do the fitting.
    import random
    random.seed(0)
    def random_xyi(n=9):
        x = list(range(n))
        y = [random.random() for i in range(n)]
        i = random.randint(0,n-1)
        return (x,y,i)

    from util.math.fraction import Fraction
    # Cycle through some randomly generated tests.
    for size in range(3, 34, 10):
        # print("size:", size)
        for trial in range(10):
            # print("  trial:", trial)
            size = 9
            trial = 0
            x,y,i = random_xyi(n=size)
            x, y = list(map(Fraction,x)), list(map(Fraction,y))
            # Get the exact value from local quadratic function.
            f = Polynomial(local_quadratic(x, y, i))
            # Try to optimize a better solution.
            opt_f = Polynomial(opt_min_fit(x, y, i))

            # Check for consistency between optimized solutoin and calculated solution.
            if any(abs(c1 - c2) > 2**(-51) for c1,c2 in zip(
                    f.coefficients, opt_f.coefficients)):
                print("ERROR: Unexpectedly large error in local quadratic fit.")
                print("f:     ",f)
                print("opt_f: ",opt_f)
                print()
                print("size:", size)
                print("  trial:", trial)
                exit()
            # elif any(abs(float(c1) - float(c2)) > 0 for c1,c2 in zip(
            #         f.coefficients, opt_f.coefficients)):
            #     print("WARNING: Did not find consistent local quadratic.")
            #     print("f:     ",f)
            #     print("opt_f: ",opt_f)


# Test the Polynomial class for basic operation.
def _test_Polynomial():
    f = Polynomial([3,0,1])
    assert(str(f) == "3 x^2  +  1")
    f = Polynomial([3,2,1])
    assert(str(f) == "3 x^2  +  2 x  +  1")
    assert(str(f.derivative()) == "6 x  +  2")
    assert(str(f.derivative(2)) == "6")
    assert(str(f.derivative(3)) == "")
    assert(str(f.derivative(4)) == "")
    assert(f.derivative(3)(10) == 0)
    f = Polynomial(NewtonPolynomial([3,2,1],[0,0,0]))
    assert(str(f) == "3 x^2  +  2 x  +  1")
    assert(str(f.derivative()) == "6 x  +  2")
    assert(str(f.derivative(2)) == "6")
    assert(str(f.derivative(3)) == "")
    assert(str(f.derivative(4)) == "")
    assert(f.derivative(3)(5) == 0)
    f = Polynomial(NewtonPolynomial([-1,10,-16,24,32,-32], [1,1,1,-1,-1,-1]))
    assert(str(f) == "-1 x^5  +  9 x^4  +  6 x^3  +  -22 x^2  +  11 x  +  -3")
    assert(str(f.derivative()) == "-5 x^4  +  36 x^3  +  18 x^2  +  -44 x  +  11")
    # Check that integrals work too.
    assert(str(f.derivative().derivative(-1)) == "-1.0 x^5  +  9.0 x^4  +  6.0 x^3  +  -22.0 x^2  +  11.0 x")
    assert(str(f.derivative().derivative(-1).derivative()) == "-5.0 x^4  +  36.0 x^3  +  18.0 x^2  +  -44.0 x  +  11.0")


# Test the Polynomial class for basic operation.
def _test_NewtonPolynomial():
    f = NewtonPolynomial([-1,2], [1,-1])
    assert(str(f) == "2 + (x + 1)(-1)")
    assert(str(Polynomial(f)) == "-1 x  +  1")
    f = NewtonPolynomial([-1,10,-16,24,32,-32], [1,1,1,-1,-1,-1])
    assert(str(f) == "-32 + (x + 1)(32 + (x + 1)(24 + (x + 1)(-16 + (x - 1)(10 + (x - 1)(-1)))))")


# Test the "polynomial" interpolation routine (uses Newton form).
def _test_polynomial(plot=True):
    SMALL = 1.4901161193847656*10**(-8) 
    # ^^ SQRT(EPSILON(REAL(1.0)))
    x_vals = [0,1,2,3,4,5]
    y_vals = [1,2,1,2,1,10]
    f = polynomial(x_vals, y_vals)
    for (x,y) in zip(x_vals,y_vals):
        try:    assert( abs(y - f(x)) < SMALL )
        except:
            string =  "\n\nFailed test.\n"
            string += f" x:    {x}\n"
            string += f" y:    {y}\n"
            string += f" f({x}): {f(x)}"
            class FailedTest(Exception): pass
            raise(FailedTest(string))


# Test the Spline class for basic operation.
def _test_Spline():
    from util.math.fraction import Fraction
    knots = [0,1,2,3,4]
    values = [[0],[1,-1,0],[0,-1],[1,0,0],[0]]
    # Create the knots and values.
    knots = [Fraction(k) for k in knots]
    values = [[Fraction(v) for v in vals] for vals in values]
    # Create the spline.
    f = Spline(knots, values)
    for (k,v) in zip(knots,values):
        for d in range(len(v)):
            try: assert(f.derivative(d)(k) == v[d])
            except:
                print()
                print('-'*70)
                print("      TEST CASE")
                print("Knot:           ", k)
                print("Derivative:     ", d)
                print("Expected value: ", v[d])
                print("Received value: ", f.derivative(d)(k))
                print()
                print(f)
                print('-'*70)
                raise(Exception("Failed test case."))
    # TODO: Test cases for the addition, multiplication, negation
    #       and integration of a Spline object.
    # 
    # from util.plot import Plot
    # g = Spline(knots, values)
    # s_add = f + g
    # add = lambda x: f(x) + g(x)
    # s_mult = f * g
    # mult = lambda x: f(x) * g(x)
    # p = Plot()
    # p.add_func("f", f, [min(knots), max(knots)])
    # p.add_func("f'", f.derivative(1), [min(knots), max(knots)])
    # p.add_func("f''", f.derivative(2), [min(knots), max(knots)])
    # p.add_func("g", g, [min(knots), max(knots)])
    # p.add_func("s add", s_add, [min(knots), max(knots)])
    # p.add_func("add", add, [min(knots), max(knots)])
    # p.add_func("s mult", s_mult, [min(knots), max(knots)])
    # p.add_func("mult", mult, [min(knots), max(knots)])
    # # Negated
    # p.add_func("-f", -f, [min(knots), max(knots)])
    # p.add_func("-f true", lambda x: -(f(x)), [min(knots), max(knots)])    
    # p.add_func("-mult", -s_mult, [min(knots), max(knots)])
    # p.add_func("-mult true", lambda x: -s_mult(x), [min(knots), max(knots)])
    # p.add_func("int(f)", f.derivative(-1), [min(knots), max(knots)])
    # p.add_func("int(-mult)", (-s_mult).derivative(-1), [min(knots), max(knots)])
    # p.show()
    # exit()

# Test the "fit" function. (there is testing code built in, so this
# test is strictly for generating a visual to verify).
def _test_fit(plot=False):
    from util.math.fraction import Fraction
    x_vals = list(map(Fraction, [0,.5,2,3.5,4,5.3,6]))
    y_vals = list(map(Fraction, [1,2,2.2,3,3.5,4,4]))
    f = fit(x_vals, y_vals, continuity=2)
    # Execute with different operational modes, (tests happen internally).
    f = fit(x_vals, y_vals, continuity=2)
    for i in range(len(f._functions)):
        f._functions[i] = Polynomial(f._functions[i])
    if plot:
        f = fit(x_vals, y_vals, continuity=1)
        print("f: ",f)
        from util.plot import Plot
        plot_range = [min(x_vals)-.1, max(x_vals)+.1]
        p = Plot()
        p.add("Points", list(map(float,x_vals)), list(map(float,y_vals)))
        p.add_func("f", f, plot_range)
        p.add_func("f'", f.derivative(1), plot_range, dash="dash")
        p.add_func("f''", f.derivative(2), plot_range, dash="dot")
        # def L2(f, low, upp):
        #     f2 = f**2
        #     p.add_func("f''<sup>2</sup>", f2, plot_range)
        #     int_f2 = f2.integral()
        #     return int_f2(upp) - int_f2(low)
        # print("L2(f''):", float(L2(f.derivative(2), x_vals[0], x_vals[-1])))
        p.show()


def _test_inverse(plot=False):
    # Construct a polynomial piece to test the inverse operation on.
    f = polynomial([1, 3], [0, 0], dx0=1, ddx0=10, dx1=-1)
    i0 = inverse(f, 0, accuracy=2**(-26))
    assert(abs(f(i0)) < 2**(-2))
    i12 = inverse(f, 1/2, accuracy=2**(-2))
    # print(float(abs(f(i12) - 1/2)), abs(f(i12) - 1/2))
    assert(abs(f(i12) - 1/2) < 2**(-2))
    if plot:
        print(f"inverse(f, 20): {float(inverse(f, 20, 4.2)):.4f}")
        im1 = inverse(f, -1)
        from util.plot import Plot
        p = Plot()
        p.add_func("f", f, [0, 5])
        p.show()

if __name__ == "__main__":
    # Run the tests on this file.
    print()
    print("Running tests..")
    print(" Polynomial")
    _test_Polynomial()
    print(" NewtonPolynomial")
    _test_NewtonPolynomial()
    print(" polynomial")
    _test_polynomial()
    print(" Spline")
    _test_Spline()
    print(" fit")
    _test_fit(plot=True)
    print(" local_quadratic")
    _test_local_quadratic()
    print(" inverse")
    _test_inverse(plot=False)
    print("tests complete.")


