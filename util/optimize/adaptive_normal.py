import numpy as np

# ====================================================================
#              DEFINITION OF ADAPTIVE NORMAL MINIMIZER     
# 
#   This is an optimization algorithm (coincidentally similar to
# annealing) that performs competitively with AMPGO on simple test
# problems. Empirically it has been observed to perform poorly in many
# dimensions.
#   This algorithm works by searching along individual dimensions
# using random normal distributions. Failure causes the standard
# deviations to shrink (and search locally). Success causes the
# standard deviations to grow (and search farther away).
# 
# 
# Adaptive Normal bibtex entry with added URL component
# 
# @inproceedings{lux2016convergence,
#   title={Convergence Rate Evaluation of Derivative-Free Optimization Techniques},
#   author={Lux, Thomas},
#   booktitle={International Workshop on Machine Learning, Optimization and Big Data},
#   pages={246--256},
#   year={2016},
#   organization={Springer},
#   url={https://link.springer.com/chapter/10.1007/978-3-319-51469-7_21},
# }
# 
# ====================================================================
ACCURACY_POW_2 = 10
DIMS_TO_MODIFY = 1
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
def AdaptiveNormal(objective, halt, bounds, solution,
                   exhaustion_limit=ACCURACY_POW_2,
                   num_dim=DIMS_TO_MODIFY, args=tuple()):
    exhaustion = 0
    obj_value = objective(solution, *args)
    # Useful variables for searching
    lower = np.array([l for l,u in bounds])
    upper = np.array([u for l,u in bounds])
    sprawl = upper - lower
    # Main search loop
    while not halt():            
        if exhaustion == 0:
            # Search globally with random solution
            new_sol = np.random.random( (len(bounds),) ) * sprawl + lower
        else:
            # Search locally with small random modifications
            for i in range(num_dim):
                param = np.random.randint( 0, len(bounds) )
                stdev = sprawl[param] / 2**exhaustion
                new_sol = solution.copy()
                new_sol[param] = np.random.normal( solution[param], stdev )
                # Put the new solution back into bounds
                new_sol[param] = max(lower[param], new_sol[param])
                new_sol[param] = min(upper[param], new_sol[param])
        # Calculate the objective value of the new solution
        new_obj_value = objective(new_sol, *args)
        if new_obj_value < obj_value:
            # If it is better, keep it and reset exhaustion
            exhaustion = 0
            solution = new_sol
            obj_value = new_obj_value
        else:
            # Otherwise increment exhaustion
            exhaustion += 1 
        # Check to see if exhaustion has gotten too large
        if exhaustion >= exhaustion_limit:
            exhaustion = 0
    return (obj_value, solution)
# ====================================================================
# ====================================================================

