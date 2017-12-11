import numpy as np

# =================================================
#           Optimization and Minimization          
# =================================================

DEFAULT_SEARCH_SIZE = 2.0**10
DEFAULT_MAX_TIME_SEC = 1
DEFAULT_MIN_STEPS = 10
DEFAULT_MIN_IMPROVEMENT = 0.0000

import time
# from .DiRect import DiRect as optimize
from .algorithms.AMPGO import AMPGO as optimize

# Class for tracking convergence of optimization algorithm
class Tracker:
    # Initialization of storage attributes
    def __init__(self, objective, max_time=DEFAULT_MAX_TIME_SEC,
                 min_steps=DEFAULT_MIN_STEPS, 
                 min_improvement=DEFAULT_MIN_IMPROVEMENT, display=False):
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

    # Wrapper for the objective function that tracks progress
    def check(self,sol):
        # Initialize starting time if it hasn't been done
        if self.start == None: self.start = time.time()
        # Don't execute the objective function if time has been exceeded
        if time.time() - self.start > DEFAULT_MAX_TIME_SEC: return float('inf')
        # Calculate the new objective value
        self.tries += 1
        new_obj = self.obj(sol)
        # Record relevant statistics
        if new_obj > self.max_obj: self.max_obj = new_obj
        if new_obj <= self.min_obj:
            self.record.append(new_obj)
            self.min_obj = new_obj
            self.best_sol = sol
            if self.display: print(new_obj, sol)
        return new_obj

    # Function that returns "True" when optimization is over
    def done(self):
        # Most importantly, check for time
        time_done = (time.time() - self.start > self.max_time)
        if time_done: return time_done
        # Next, check for convergence
        converged = False
        if len(self.record) > self.min_steps:
            change = [abs(v - self.min_obj) for v in self.record[-self.min_steps:]]
            divisor = max(self.record) - min(self.record)
            if divisor == 0: 
                converged = False
            else:
                converged = max(change) / divisor < self.min_change
        return converged

# Minimization function
def minimize(objective, solution, bounds=None,
             max_time=DEFAULT_MAX_TIME_SEC,
             min_steps=DEFAULT_MIN_STEPS,
             min_improvement=DEFAULT_MIN_IMPROVEMENT, plot=False):
    solution = np.array(list(map(float,solution)))

    # Generate some arbitrary bounds
    if type(bounds) == type(None):
        upper = solution + np.array([DEFAULT_SEARCH_SIZE]*len(solution))
        lower = solution - np.array([DEFAULT_SEARCH_SIZE]*len(solution))
        bounds = list(zip(lower, upper))

    # Initialize a tracker for halting the optimization
    t = Tracker(objective, max_time, min_steps, min_improvement, plot)

    # Call the optimization function and get the best solution
    optimize(t.check, t.done, bounds, solution)

    if plot:
        import pygal
        CUSTOM_CSS_FILE = 'file://'+os.path.expanduser("~")+os.path.join(
            "Git", "pyvis", "pyvis", "Visualize",
            "Value_Distribution", "no_dot.css")
        config = pygal.Config()
        config.css.append(CUSTOM_CSS_FILE)
        plot = pygal.XY(config)
        plot.add("Objective", list(zip(range(len(t.record)),t.record)))
        plot.render_in_browser()
        

    return t.best_sol

