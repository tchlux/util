import numpy as np
from scipy.optimize import minimize as scipy_minimize

# ====================================================================
#                  DEFINITION OF AMPGO MINIMIZER     
# 
#   AMPGO is a general framework that uses an arbitrary local convex
# optimization technique to identify globally optimal solutions in a
# multimin search space by "tunneling" (modifying the objective
# function to avoid previously identified optimal solutions).
#   This implementation utilizes the minimize function from
# scipy.optimize and the "BFGS" minimization technique that
# approximates the hessian of the objective function (2nd derivative
# with respect to all dimensions) and jumps to the absolute minimum.
# 
# 
# AMPGO bibtex entry with added URL component
# 
# @article{lasdon2010adaptive,
#   title={Adaptive memory programming for constrained global optimization},
#   author={Lasdon, Leon and Duarte, Abraham and Glover, Fred and Laguna, Manuel and Mart{\'\i}, Rafael},
#   journal={Computers \& Operations Research},
#   volume={37},
#   number={8},
#   pages={1500--1509},
#   year={2010},
#   publisher={Elsevier},
#   url={http://leeds-faculty.colorado.edu/glover/Publications/TS%20-%20AMP%20GO%20Constrained%20Global%20Opt%20w%20Lasdon%20et%20al%20.pdf},
# }
# 
# ====================================================================
EPS1 = 0.02
EPS2 = 0.1
MAX_TABU_ITER = 5
MAX_TABULIST_SIZE = 5
LOCAL_TOLERANCE = 1e-8
# Local minimum solver function, returns new solution and obj. value
def LOCAL_SOLVER(objfun, initial, bounds, args=()):
    result = scipy_minimize(objfun, initial, args=args, method="L-BFGS-B",
                            bounds=bounds, tol=LOCAL_TOLERANCE)
    return result["x"], result["fun"]
# Sub-problem used to escape local minima
def TUNNEL(solution, objective, aspiration, tabulist):
    numerator = (objective(solution) - aspiration)**2
    denominator = 1.0
    for tabu in tabulist:
        denominator *= np.sqrt(np.sum((solution - tabu)**2))
    if denominator == 0: return float('inf')
    return numerator / denominator
# Drop the point that is farthest away from the new point
def DROP_TABU_POINTS(solution, tabulist):
    if len(tabulist) >= MAX_TABULIST_SIZE:
        # Drop the farthest points
        distance = np.sqrt(np.sum(
            (tabulist - solution)**2, axis=1))
        index = np.argmax(distance)
        tabulist.pop(index)
    tabulist.append(solution)

# Pre:  "objective" a numpy array and returns an object that has a
#         comparison operator 
#       "halt" is a function that takes no parameters and returns
#         "False" when the algorithm should continue and "True" when
#         it should halt
#       "bounds" is a list of tuples (lower,upper) for each parameter
#       "solution" is a valid initial solution for this algorithm (in
#         numpy array format) and can be modified with expected side
#         effects outside the scope of this algorithm 
#       "args" is a set of arguments to be passed to the objective
#         function after the first parameter (a solution)
# Post: Returns an (objective value, solution) tuple, best achieved
def AMPGO(objective, halt, bounds, solution,
          max_tabu_iter=MAX_TABU_ITER, args=tuple()):
    calc_obj = lambda sol: objective(sol, *args)
    tabulist = [solution]
    lower = np.array([b[0] for b in bounds])
    upper = np.array([b[1] for b in bounds])
    best_obj = calc_obj(solution)
    best_sol = solution.copy()
    new_sol, obj_val = LOCAL_SOLVER(calc_obj, solution, bounds)
    # Record this solution if it is the best
    if obj_val < best_obj:
        best_sol = new_sol.copy()
        best_obj = obj_val
    while not halt():
        # Drop old points and add the new one
        DROP_TABU_POINTS(new_sol, tabulist)
        # Use an altered objective function and randomization to escape
        # current local minimum, attempt this "max_tabu_iter" times
        for i in range(max_tabu_iter):
            if halt(): break
            # Generate a random new solution by picking a random search
            # direction. Walk "EPS2"*length of current solution.
            r = np.random.uniform(-1.0, 1.0, size=(solution.shape[0],))
            beta = EPS2 * np.linalg.norm(new_sol) / np.linalg.norm(r)
            if np.abs(beta) < LOCAL_TOLERANCE: beta = EPS2
            new_sol = new_sol + beta*r
            # Put the new solution back into bounds
            new_sol = np.where(new_sol < lower, lower, new_sol)
            new_sol = np.where(new_sol > upper, upper, new_sol)
            # Generate the aspired objective function value
            aspiration = best_obj - EPS1*(1.0 + np.abs(best_obj))
            tunnel_args = (calc_obj, aspiration, tabulist)
            # Run local solver on sub-problem ("_" is because the obj
            # value is not the actual obj_value for the function)
            new_sol, _ = LOCAL_SOLVER(TUNNEL, new_sol, bounds,
                                      args=tunnel_args)
            obj_val = calc_obj(new_sol)
            if obj_val < best_obj:
                best_sol = new_sol.copy()
                best_obj = obj_val
                break
    return best_obj, best_sol
# ====================================================================
# ====================================================================
