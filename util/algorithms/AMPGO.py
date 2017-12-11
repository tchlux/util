import numpy
from scipy.optimize import minimize

EPS1 = 0.02
EPS2 = 0.1
MAX_ITER = 100
MAX_TABU_ITER = 5
MAX_TABULIST_SIZE = 5
LOCAL_TOLERANCE = 1e-8

# Local minimum solver function, returns new solution and obj. value
def LOCAL_SOLVER(objfun, initial, bounds, maxiter, args=()):
    # options = {'maxiter':max(1, maxiter)}
    result = minimize(objfun, initial, args=args, method="L-BFGS-B",
                      bounds=bounds, tol=LOCAL_TOLERANCE)# ,
                      # options=options)
    return result["x"], result["fun"]

# Sub-problem used to escape local minima
def TUNNEL(solution, objective, aspiration, tabulist):
    numerator = (objective(solution) - aspiration)**2
    denominator = 1.0
    for tabu in tabulist:
        denominator *= numpy.sqrt(numpy.sum((solution - tabu)**2))
    return numerator / denominator

# Drop the point that is farthest away from the new point
def DROP_TABU_POINTS(solution, tabulist):
    if len(tabulist) >= MAX_TABULIST_SIZE:
        # Drop the farthest points
        distance = numpy.sqrt(numpy.sum(
            (tabulist - solution)**2, axis=1))
        index = numpy.argmax(distance)
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
#       "objargs" is a set of arguments to be passed to the objective
#         function after the first parameter (a solution)
# Post: Returns an (objective value, solution) tuple, best achieved
def AMPGO(objective, halt, bounds, solution, maxiter=MAX_ITER,
          max_tabu_iter=MAX_TABU_ITER, *objargs):
    calc_obj = lambda sol: objective(sol, *objargs)
    tabulist = [solution]
    lower = numpy.array([b[0] for b in bounds])
    upper = numpy.array([b[1] for b in bounds])
    best_obj = calc_obj(solution)
    best_sol = solution.copy()
    new_sol, obj_val = LOCAL_SOLVER(calc_obj, solution, 
                                    bounds, maxiter)
    # Record this solution if it is the best
    if obj_val < best_obj:
        best_sol = new_sol.copy()
        best_obj = obj_val

    while not halt():
        # Drop old points and add the new one
        DROP_TABU_POINTS(new_sol, tabulist)

        # Use an altered objective function and randomization to
        # escape current local minimum
        for i in range(max_tabu_iter):
            if halt(): break
            # Generate a random new solution
            r = numpy.random.uniform(-1.0, 1.0, size=(solution.shape[0],))
            beta = EPS2 * numpy.linalg.norm(new_sol) / numpy.linalg.norm(r)
            if numpy.abs(beta) < LOCAL_TOLERANCE:
                beta = EPS2
            new_sol = new_sol + beta*r

            # Put the new solution back into bounds
            new_sol = numpy.where(new_sol < lower, lower, new_sol)
            new_sol = numpy.where(new_sol > upper, upper, new_sol)

            # Generate the aspired objective function value
            aspiration = best_obj - EPS1*(1.0 + numpy.abs(best_obj))
            tunnel_args = (calc_obj, aspiration, tabulist)
            # Run local solver on sub-problem ("_" is because the obj
            # value is not the actual obj_value for the function)
            new_sol, _ = LOCAL_SOLVER(TUNNEL, new_sol, bounds,
                                      max_tabu_iter, args=tunnel_args)
            obj_val = calc_obj(new_sol)
            if obj_val < best_obj:
                best_sol = new_sol.copy()
                best_obj = obj_val
                break

    return best_obj, best_sol
