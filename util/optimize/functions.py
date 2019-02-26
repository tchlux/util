# Get numpy data structures (in differentiable form)
from autograd.numpy import zeros, ones, arange, linspace, \
    identity, roll, radians, random, linalg
# Get numpy mathematical expressions (in differentiable form)
from autograd.numpy import pi, sum, exp, sin, cos, log, outer, dot, matmul
# Use an automatic elementwise-gradient function
from autograd import elementwise_grad as grad

# Set up a repeatable random number generator
SEED = 0
RANDOM = random.RandomState(SEED)
# Store all the objective functions in a global variable
FUNCTIONS = []
functions = FUNCTIONS


# Decorator for converting a function into an objective function
def objectify(lower=None, upper=None, max_val=None, solution=None, gradient=None):
    # Define the "normalize" decorator that handles given arguments
    def decorator(func):
        def obj_func(x, *args, **kwargs):
            # Rescale the x out of [-1,1] and into the domain of the
            # function if a normalization was necessary
            if (type(lower) != type(None)) and (type(upper) != type(None)):
                x = ((x + 1) * (upper - lower) / 2) + lower
            # Return the objective value
            value = func(x, *args, **kwargs)
            # Rescale the objective value if necessary
            if (type(max_val) != type(None)): value /= max_val*len(x)
            return value
        # Set a random point generator (in the -1,1 range)
        obj_func.rand = lambda d, gen=RANDOM: gen.rand(d) * 2 - 1
        # Set a solution
        if solution: obj_func.sol  = solution
        else:        obj_func.sol  = lambda d: zeros(d)
        # Set a gradient
        if gradient: obj_func.grad = gradient
        else:        obj_func.grad = grad(obj_func)
        obj_func.lower, obj_func.upper = -1., 1.
        # Store the maximum gradient magnitude
        obj_func.max_grad = max(abs(obj_func.grad(ones(1)*obj_func.lower)),
                                abs(obj_func.grad(ones(1)*obj_func.upper)))
        # Copy the name over to the new function
        obj_func.__name__ = func.__name__ + "_(objectified)"
        obj_func.func = func
        # Add the objective function to the list of functions
        FUNCTIONS.append(obj_func)
        # Return the new objective function
        return obj_func
    # Return the decoroator, which will be passed a single function as
    # an argument by python implicitly.
    return decorator


# # A standard quadratic function.
# @objectify(max_val=1)
# def quadratic(x):
#     return sum(x**2)

# A sub-quadratic with sudden increase in steepness near solution at 0.
@objectify(max_val=1)
def mejr(x, k=2):
    return sum(abs(x)**((2*k)/(2*k-1)))

# A super-quadratic with apparent 'flatness' around solution at 0.
@objectify(max_val=1)
def basin(x, k=2):
    return sum((x)**(2*k))

# A polynomial with saddle points at +/- "s", solution at 0.
@objectify(max_val=(1/6 - .75**2 / 2 + .75**4 / 2))
def saddle(x, s=.75):
    return sum((s**4)*(x**2)/2 - (s**2)*(x**4)/2 + (x**6)/6)

# A chebyshev polynomial with "m" minima, solution at 0.
@objectify(max_val=4)
def multimin(x, m=3, a=2):
    # Compute the chebyshev polynomial with "m" minima
    previous = 1
    current = x
    for n in range(2, m*2 + 1):
        previous, current = (current, 2*x*current - previous)
    # Return the chebyshev polynomial with a quadratic overture
    return len(x) + a*sum(x**2) + sum(current)

# An approximation of a nowhere differentiable function.
@objectify(max_val=3.998046875)
def weirerstrass(x, k=10, a=.5, b=3):
    # Calculate some constants
    ak = a**(arange(k+1))
    bk = b**(arange(k+1))
    # Calculate the function value
    ak_xi = cos( pi * outer(bk,(x+1)) )
    return sum(dot(ak, ak_xi)) - len(x) * sum(ak * cos(pi*bk))

