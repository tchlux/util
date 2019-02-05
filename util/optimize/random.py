import numpy as np

# ====================================================================
#                     A PURELY RANDOM MINIMIZER
# 
#   Generates random solutions inside of bounds, keeping track of best
# solution obtaiend.
# 
# Pre:  "objective" a numpy array and returns an object that has a
#         comparison operator 
#       "halt" is a function that takes no parameters and returns
#         "False" when the algorithm should continue and "True" when
#         it should halt
#       "bounds" is a list of tuples (lower,upper) for each parameter
#       "solution" is a valid initial solution for this algorithm (in
#         numpy array format) and can be modified with expected side
#         effects outside the scope of this algorithm 
# Post: Returns an (objective value, solution) tuple, best achieved
def Random(objective, halt, bounds, solution, args=tuple()):
    bound_range = np.array([u-l for u,l in bounds])
    bound_shift = np.array([l for u,l in bounds])    
    random_solution = lambda: ( np.random.random( (len(bounds),) )
                                * bound_range + bound_shift )
    best_solution = random_solution()
    best_obj = objective(best_solution, *args)
    while not halt():
        new_solution = random_solution()
        new_obj = objective(new_solution, *args)
        if new_obj < best_obj:
            best_solution = new_solution
            best_obj = new_obj
    return (best_obj, best_solution)
# ====================================================================
# ====================================================================
