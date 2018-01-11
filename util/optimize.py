import time, random
import numpy as np
from scipy.optimize import minimize as scipy_minimize

CHECKPOINT_FILE = "optimization_checkpoint.py"

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
MAX_ITER = 100
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
# standard deviations to groaw (and search farther away).
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

# ====================================================================
#            Generic Optimization and Minimization Tracking          
# 
DEFAULT_SEARCH_SIZE = 2.0**10
DEFAULT_MAX_TIME_SEC = 0.5
DEFAULT_MIN_STEPS = 10
DEFAULT_MIN_IMPROVEMENT = 0.0
# Class for tracking convergence of optimization algorithm
class Tracker:
    # Initialization of storage attributes
    def __init__(self, objective, max_time=DEFAULT_MAX_TIME_SEC,
                 min_steps=DEFAULT_MIN_STEPS, 
                 min_improvement=DEFAULT_MIN_IMPROVEMENT,
                 display=False, checkpoint=True,
                 checkpoint_file=CHECKPOINT_FILE):
        self.max_time = max_time
        self.min_steps = min_steps
        self.min_change = min_improvement
        self.display = display
        self.obj = objective
        self.max_obj = -float('inf')
        self.min_obj = float('inf')
        self.best_sol = None
        self.record = []
        self.tries = 0
        self.start = None
        self.checkpoint = checkpoint
        self.checkpoint_file = checkpoint_file

    # Wrapper for the objective function that tracks progress
    def check(self, sol, *args, **kwargs):
        # Initialize starting time if it hasn't been done
        if self.start == None: self.start = time.time()
        # Don't execute the objective function if stopping criteria met
        if self.done(): return float('inf')
        # Calculate the new objective value
        self.tries += 1
        new_obj = self.obj(sol, *args, **kwargs)
        # Record relevant statistics
        if new_obj > self.max_obj: self.max_obj = new_obj
        # Default to having displays be inline updates, only changed
        # for when new best solutions are discovered
        display_end = "\r"
        # Only record the new best solution if it is (<= best) and unique
        if (new_obj <= self.min_obj):
            old_min_obj = self.min_obj
            # Recording the new solutions (record ones that are just as good as others)
            self.record.append((new_obj, self.tries, sol))
            self.min_obj = new_obj
            self.best_sol = sol
            # Checkpointing and displaying solution progression should
            # only happen if the new solution is actually better.
            if (new_obj < old_min_obj): 
                # Checkpoints for progression (always save the best solutions)
                if self.checkpoint:
                    formatted_sol = ", ".join(list(map(lambda f:"%e"%f,sol)))
                    with open(self.checkpoint_file, "w") as f:
                        print("solution = [%s]"%(formatted_sol), file=f)
                        print("obj_value = %s"%(self.min_obj), file=f)
                display_end = "\n"
        # Display progress to the user
        if self.display:
            if (self.tries == 1): 
                print()
                print("Formatted convergence printout:")
                print("","[delay] ([attempts]):\t[obj] with numbered solutions saved in %s"%self.checkpoint_file)
                print()
            if (len(self.record) == 1): delay = self.tries
            else: delay = self.tries - self.record[-2][1]
            print("","%i (%i): \t%.1e"%(delay, self.tries, new_obj), end=display_end)
        return new_obj

    # Function that returns "True" when optimization is over
    def done(self):
        # Most importantly, check for time
        time_done = (time.time() - self.start > self.max_time)
        if time_done: return time_done
        # Next, check for convergence
        converged = False
        if len(self.record) > self.min_steps:
            change = [abs(o - self.min_obj) for (o,t,s) in self.record[-self.min_steps:]]
            divisor = max(self.record)[0] - min(self.record)[0]
            if divisor == 0: 
                # This means max_obj == min_obj, clearly not converged
                converged = False
            else:
                # Check if the relative convergence of the most recent
                # window of solutions is slower than a prescribed value
                converged = max(change) / divisor < self.min_change
        return converged
# ====================================================================
# ====================================================================


# Minimize arbitrary objective function, optionally given an initial
# solution, bounds, max computation time, minimum number of steps, or
# a minimum improvement stopping criteria.
# 
# INPUT:
#   objective -- A callable function that takes a numpy array of
#                the same length as "solution" and returns a float.
#   solution  -- A 1D numpy array (or list) that is an initial
#                starting point for the search (also implicilty
#                defines the number of dimenions for the problem)
# 
# OPTIONAL INPUT:
#   bounds          -- [(low_1, upp_1), ..., (low_d, upp_d)]
#                      A list of tuples of lower and upper bounds for
#                      each dimension of valid candidate solutions
#   args            -- A tuple of arguments to pass to the objectivev function
#   max_time        -- The maximum amount of computation time for optimization
#   min_steps       -- The minimum number of evaluations of the objective function
#   min_improvement -- The minimum improvement needed to continue searching
#   display         -- Boolean, if True then try to use plotly to render
#                      a plot of the best solution convergence by
#                      iterations, otherwise print a console summary.
#   method          -- A minimization function that takes as
#                      parameters, (objective function, halting
#                      function , bounds, initial solution, obj args)
#                      and will find a candidate minimum solution.
#      This file provides the methods:
#       AMPGO (default) - BFGS with tunneling away from prior solutions,
#                         best for smooth locally convex multi-min spaces
#       AdaptiveNormal  - Quasi simulated annealing with exhaustion,
#                         best for noisy multi-min spaces
#       Random          - Pure random baseline for comparison, best for
#                         discontinuous and noisy spaces.
#   checkpoint      -- Boolean indicating whether or not a checkpoint
#                      file should be saved with intermediate best solutions
#   checkpoint_file -- String name of the file to save for checkpointing
# 
# OUTPUT:
#   solution        -- The best solution identified by the minimizer provided
#                      the initial restraints + stopping condition(s).
def minimize(objective, solution, bounds=None, args=tuple(),
             max_time=DEFAULT_MAX_TIME_SEC,
             min_steps=DEFAULT_MIN_STEPS,
             min_improvement=DEFAULT_MIN_IMPROVEMENT, display=False,
             method=AdaptiveNormal, checkpoint=True,
             checkpoint_file=CHECKPOINT_FILE):
    # Convert the solution into a float array (if it's not already)
    solution = np.asarray(solution, dtype=float)

    # Generate some arbitrary bounds
    if type(bounds) == type(None):
        upper = solution + np.ones((len(solution),))*DEFAULT_SEARCH_SIZE
        lower = solution - np.ones((len(solution),))*DEFAULT_SEARCH_SIZE
        bounds = list(zip(lower, upper))

    # Initialize a tracker for halting the optimization
    t = Tracker(objective, max_time, min_steps, min_improvement,
                display, checkpoint, checkpoint_file)

    # Call the optimization function and get the best solution

    method(t.check, t.done, bounds, solution, args=args)

    if display: print()
    # Remove the checkpoint file
    if os.path.exists(checkpoint_file): os.remove(checkpoint_file)
    return t.best_sol

