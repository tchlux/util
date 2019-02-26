from util.optimize.functions import *

# ============================================
#      Objective Function Transformations     
# ============================================

# Errors that may be rasied
class IllegalNoise(Exception): pass
class IllegalSkew(Exception): pass
class IllegalRotation(Exception): pass
class IllegalRecenter(Exception): pass

# Given a noise ratio in [0,1], add noise to the gradient that is
# proportional to "ratio" times the max magnitude of the gradient.
def noise(ratio, dimension):
    # Check for proper usage
    if not (0 <= ratio <= 1): 
        raise(IllegalNoise("Only noise ratios in the range [0,1] are accepted."))
    # Create the decorator for an objective function
    def decorator(func):
        # Calculate the magnitude of the noise that will be added
        noise_magnitude = ratio * func.max_grad
        # Generate a new copy of the objective function (different gradient)
        def obj_func(x, *args, **kwargs):
            return func(x , *args, **kwargs)

        # Store the solution, gradient, lower, and upper bounds for this new function
        obj_func.rand = func.rand
        obj_func.sol = func.sol
        obj_func.grad = lambda x: func.grad(x) + (RANDOM.rand(dimension)-.5)*noise_magnitude/2
        obj_func.max_grad = func.max_grad
        obj_func.lower = ones(dimension) * func.lower
        obj_func.upper = ones(dimension) * func.upper
        obj_func.__name__ = func.__name__ + f"_(noised-{ratio:.2f})"
        obj_func.func = func
        # Return new function
        return obj_func
    # Return decorator function
    return decorator    

# Given a skew in [0,1) go from the original space at 0 to 1
# where the difference in stretched-ness between the least and most
# stretched dimension approaches infinity.
def skew(skew, dimension):
    # Check for proper usage
    if not (0 <= skew < 1): 
        raise(IllegalSkew("Only skews in the range [0,1) are accepted."))
    # Generate a skew vector to apply 
    skew_vec = linspace(1-skew, 1, dimension)
    # Create the decorator for an objective function
    def decorator(func):
        def obj_func(x, *args, **kwargs):
            return func(x * skew_vec, *args, **kwargs)
        # Store the new random point generation function
        obj_func.rand = lambda d, gen=RANDOM: func.rand(d,gen) / skew_vec
        # Store the new solution
        obj_func.sol = lambda d: func.sol(d) / skew_vec
        # Store the new gradient
        obj_func.grad = lambda x: func.grad(x*skew_vec) / skew_vec
        obj_func.max_grad = func.max_grad
        # Store the new bounds
        obj_func.lower = (ones(dimension)*func.lower) / skew_vec
        obj_func.upper = (ones(dimension)*func.upper) / skew_vec
        # Copy over the name
        obj_func.__name__ = func.__name__ + f"_(skewed-{skew:.2f})"
        obj_func.func = func
        # Return new function
        return obj_func
    # Return decorator function
    return decorator

# Given a rotation in [0,1] go from the original space as 0 to 1 when
# the rotation corresponds to all axes being rotated to be exactly
# between the original dimensions.
def rotate(rotation_ratio, dimension):
    # Check for proper usage
    if not (0 <= rotation_ratio <= 1): 
        raise(IllegalRotation("Only rotations in the range [0,1] are accepted."))
    # Convert the rotation from degrees to radians
    rotation = radians(rotation_ratio * 45)
    # Determine the order of the dimension rotations (done in pairs)
    rotation_order = arange(dimension)
    # Generate the rotation matrix by rotating two dimensions at a time.
    rotation_matrix = identity(dimension)
    for (d1, d2) in zip(rotation_order, roll(rotation_order,1)):
        next_rotation = identity(dimension)
        # Establish the rotation
        next_rotation[d1,d1] =  cos(rotation)
        next_rotation[d2,d2] =  cos(rotation)
        next_rotation[d1,d2] =  sin(rotation)
        next_rotation[d2,d1] = -sin(rotation)
        # Compound the paired rotations
        rotation_matrix = matmul(next_rotation, rotation_matrix)
        # When there are two dimenions or fewer, do not keep iterating.
        if (dimension <= 2): break
    inverse_rotation = linalg.inv(rotation_matrix)
    # Create the decorator for an objective function
    def decorator(func):
        def obj_func(x, *args, **kwargs):
            return func(matmul(rotation_matrix, x), *args, **kwargs)
        # Store the new random point generation function
        obj_func.rand = lambda d, gen=RANDOM: matmul(inverse_rotation,func.rand(d,gen))
        # Store the new solution
        obj_func.sol = lambda d: matmul(inverse_rotation,func.sol(d))
        # Store the new gradient function
        obj_func.grad = lambda x: matmul(inverse_rotation,func.grad(matmul(rotation_matrix,x)))
        obj_func.max_grad = func.max_grad
        # Store the new bounds
        obj_func.lower = matmul(rotation_matrix, ones(dimension) * func.lower)
        obj_func.upper = matmul(rotation_matrix, ones(dimension) * func.upper)
        obj_func.lower[:] = min(obj_func.lower)
        obj_func.upper[:] = max(obj_func.upper)
        # Copy over the name
        obj_func.__name__ = func.__name__ + f"_(rotated-{rotation_ratio:.2f})"
        obj_func.func = func
        # Return new function
        return obj_func
    # Return decorator function
    return decorator

# Given a new center, shift all inputs to be centered about that point.
def recenter(new_center, dimension):
    # Create the decorator for an objective function
    def decorator(func):
        # Calculate the shift based on the new center
        shift = func.sol(dimension) - new_center
        # Check for proper usage
        if not (all(func.lower <= new_center) and 
                all(new_center <= func.upper)):
            raise(IllegalRecenter("New center must be within existing bounds."))
        def obj_func(x, *args, **kwargs):
            return func(x + shift, *args, **kwargs)
        # Store the solution, gradient, lower, and upper bounds for this new function
        obj_func.rand = lambda d, gen=RANDOM: func.rand(d,gen) - shift
        obj_func.sol = lambda d: func.sol(d) - shift
        obj_func.grad = lambda x: func.grad(x + shift)
        obj_func.max_grad = func.max_grad
        obj_func.lower = ones(dimension) * func.lower - shift
        obj_func.upper = ones(dimension) * func.upper - shift
        obj_func.__name__ = func.__name__ + f"_(recentered-{new_center})"
        obj_func.func = func
        # Return new function
        return obj_func
    # Return decorator function
    return decorator

