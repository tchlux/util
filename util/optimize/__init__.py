# Global libraries needed.
import time, os
import numpy as np
# Algorithm for one-dimensional optimization.
from util.optimize.line import secant as zero_on_line
from util.optimize.line import golden_section as min_on_line
from util.optimize.line import binary as binary_search
# Algorithms for multi-dimensional optimization.
from util.optimize.gradient_free import AMPGO, AdaptiveNormal, DiRect, Random
from util.optimize.gradient_based import L_BFGS, SGD, ADAGRAD, ADAM

# Optimization default settings.
DEFAULT_SEARCH_SIZE = 2.0**(10)
DEFAULT_MAX_TIME_SEC = 2.0**(0)
DEFAULT_MIN_STEPS = 10
DEFAULT_MAX_STEPS = 10000
DEFAULT_MIN_IMPROVEMENT = 0.0

# Tracker settings.
CHECKPOINT_FILE = "optimization_checkpoint.py"

# ====================================================================
#            Generic Optimization and Minimization Tracking          
# 
# Class for tracking convergence of optimization algorithm
class Tracker:
    # Initialization of storage attributes
    def __init__(self, objective,
                 max_time=DEFAULT_MAX_STEPS,
                 min_steps=DEFAULT_MIN_STEPS,
                 max_steps=DEFAULT_MAX_STEPS,
                 min_improvement=DEFAULT_MIN_IMPROVEMENT,
                 display=False, checkpoint=False,
                 checkpoint_file=CHECKPOINT_FILE):
        self.max_time = max_time
        self.min_steps = min_steps
        self.max_steps = max_steps
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
            self.best_sol = sol.copy()
            # Checkpointing and displaying solution progression should
            # only happen if the new solution is actually better.
            if (new_obj < old_min_obj): 
                # Checkpoints for progression (always save the best solutions)
                if self.checkpoint:
                    formatted_sol = ", ".join(list(map(lambda f:"%e"%f,self.best_sol)))
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
        # Check for convergence (or other stopping criteria)
        converged = False
        if (self.tries >= self.min_steps):
            # Check for max steps requirement.
            if (self.tries >= self.max_steps):
                if self.display:
                    print("Attempt limit exhausted.")
                    self.display = False
                return True
            # Check for time requirement
            if (time.time() - self.start > self.max_time):
                if self.display:
                    print("Time limit exhausted.")
                    self.display = False
                return True
            # Check for change requirement
            change = [abs(o - self.min_obj) for (o,t,s) in self.record[-self.min_steps:]]
            divisor = max(self.record)[0] - min(self.record)[0]
            if divisor == 0: 
                # This means max_obj == min_obj, clearly not converged
                converged = False
            else:
                # Check if the relative convergence of the most recent
                # window of solutions is slower than a prescribed value
                converged = max(change) / divisor < self.min_change
        if (converged and self.display):
            print("Achieved convergence limit.")
            self.display = False
        return converged
# ====================================================================
# ====================================================================


# Minimize arbitrary objective function given an initial solution and
# optionally given bounds, max computation time, minimum number of
# steps, or a minimum improvement stopping criteria.
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
#   max_steps       -- The maximum number of evaluations of the objective function
#   min_improvement -- The minimum improvement needed to continue searching
#   display         -- Boolean, if True then try to use plotly to render
#                      a plot of the best solution convergence by
#                      iterations, otherwise print a console summary.
#   method          -- A minimization function that takes as
#                      parameters, (objective function, halting
#                      function , bounds, initial solution, obj args)
#                      and will find a candidate minimum solution.
#      This file provides the methods:
#       AdaptiveNormal - Quasi simulated annealing with exhaustion,
#            (default)   best for noisy multi-min spaces.
#       DiRect         - The divided rectangles method, provably convergent
#                        at a slow rate, offers "region" of best solution.
#       AMPGO          - BFGS with tunneling away from prior solutions,
#                        best for smooth locally convex multi-min spaces
#       Random         - Pure random baseline for comparison, best for
#                        discontinuous and noisy spaces.
#   checkpoint      -- Boolean indicating whether or not a checkpoint
#                      file should be saved with intermediate best solutions
#   checkpoint_file -- String name of the file to save for checkpointing
# 
# OUTPUT:
#   solution        -- The best solution identified by the minimizer provided
#                      the initial restraints + stopping condition(s).
def minimize(objective, solution, bounds=None, args=tuple(),
             max_time=DEFAULT_MAX_TIME_SEC, min_steps=DEFAULT_MIN_STEPS,
             max_steps=DEFAULT_MAX_STEPS,  min_improvement=DEFAULT_MIN_IMPROVEMENT, 
             display=True, method=AMPGO, checkpoint=True,
             checkpoint_file=CHECKPOINT_FILE):
    # Convert the solution into a float array (if it's not already)
    solution = np.asarray(solution, dtype=float)

    # Generate some arbitrary bounds
    if type(bounds) == type(None):
        upper = solution + np.ones((len(solution),))*DEFAULT_SEARCH_SIZE
        lower = solution - np.ones((len(solution),))*DEFAULT_SEARCH_SIZE
        bounds = list(zip(lower, upper))

    # Initialize a tracker for halting the optimization
    t = Tracker(objective, max_time, min_steps, max_steps, min_improvement,
                display, checkpoint, checkpoint_file)
    # Get the initial objective value of the provided solution.
    t.check(solution, *args)

    # Call the optimization function and get the best solution
    method(t.check, t.done, bounds, solution, args=args)

    if display: 
        print()
        try:
            from util.plot import Plot
            name = objective.__name__.title()
            p = Plot(f"Minimization Performance on Objective '{name}'",
                     "Trial Number", "Objective Value")
            trial_numbers = [n for (o,n,s) in t.record]
            obj_values = [o for (o,n,s) in t.record]
            p.add("", trial_numbers, obj_values, color=p.color(1),
                  mode="lines+markers") 
            p.show(show_legend=False)
        except: pass
            
    # Remove the checkpoint file if checkpoint is turned off.
    if (not checkpoint) and os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)
    return t.best_sol

