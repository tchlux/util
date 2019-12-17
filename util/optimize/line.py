
# Perform a binary search on a line over a boolean objective function.
# Find the point at which the the function switches outputs.
# 
# INPUTS:
#   f -- A boolean function of the form f(x) = (x <= z) for (l <= z < u).
#   l -- The lower bound for the range to search, f(l) = True.
#   u -- The upper bound for the range to search, f(u) = False.
# 
# OPTIONAL INPUTS:
#   accuracy -- The maximum allowable distance from the correct solution.
#   round    -- A rounding function to set the values of l and u.
# 
# RETURNS:
#   x -- A point within `accuracy` of the position where f transitions
#        from True to False.
#    
def binary_search(f, l, u, accuracy=2**(-26), round=lambda x: x):
    # Verify the inputs provided are valid.
    assert(f(l) == True)
    assert(f(u) == False)
    assert(l <= u)
    assert(accuracy > 0)
    # Perform the binary search until we reach desired accuracy.
    while (abs(u - l) > accuracy):
        mid = round(l/2 + u/2)
        # If we've reached the limit of achievable accuracy, return.
        if (mid == l) or (mid == u): return mid
        # Transition the search window.
        if f(mid): l = mid
        else:      u = mid
    # Return the midpoint of the upper and lower bounds.
    return round(l/2 + u/2)



# Perform a Golden section search on a line over a convex 1D
# objective function. Find the minimum point of the function.
# 
# INPUTS:
#   f -- A real valued function with a single minimum on (l,u).
#   l -- The lower bound for the range to search.
#   u -- The upper bound for the range to search.
# 
# OPTIONAL INPUTS:
#   accuracy -- The maximum allowable distance from the correct solution.
#   round    -- A rounding function to set the values of l and u.
# 
# RETURNS:
#   x -- A point within `accuracy` of the position that minimizes f.
#    
def golden_section_search(f, l, u, accuracy=2**(-26), round=lambda x: x):
    # Verify the inputs provided are valid.
    assert(l <= u)
    assert(accuracy > 0)
    # Make sure the lower and upper bounds are rounded.
    l, u = round(l), round(u)
    ratio = (5**(1/2) + 1) / 2
    # Initialize the two midpoints.
    d = (u - l) / ratio
    m1 = round(u - d)
    m2 = round(l + d)
    # Compute the four initial function values.
    fl = f(l)
    fm1 = f(m1)
    fm2 = f(m2)
    fu = f(u)
    # Perform the Golden section search.
    while (abs(u - l) > accuracy):
        # Based on the midpoint evaluations, decide which side to search.
        if (fm1 <= fm2):
            # In this case, the minimum must be on the left side.
            m2, u = m1, m2
            fm2, fu = fm1, fm2
            m1 = round(u - (u - l) / ratio)
            # Break if we have reached the limit of achievable accuracy.
            if m1 in (l,u): break
            # Otherwise compute the new value and move forward.
            fm1 = f(m1)
        else:
            # In this case, the minimum must be on the right side.
            l, m1 = m1, m2
            fl, fm1 = fm1, fm2
            m2 = round(l + (u - l) / ratio)
            # Break if we have reached the limit of achievable accuracy.
            if m2 in (l,u): break
            # Otherwise compute the new value and move forward.
            fm2 = f(m2)
    # Return the midpoint of the upper and lower bounds.
    return round(l/2 + u/2)


# ====================================================================
#                          SECANT METHOD
# 
# Given a function f: R -> R, such that f(x) = 0 exactly once for 
#  some (low <= x <= upp), this function will return an x_0 such that
#  abs(x_0 - x) < prec using the Secant method.
# 
def Secant(f, low, upp, prec=2**(-2**1023), max_steps=1000, display=False):
    # Compute the "low" function value and the "upp" function value.
    flow = abs(f(low))
    fupp = abs(f(upp))
    # Compute the new guess for the zero by linearly interpolating.
    low_wt = fupp / (flow + fupp)
    new = low_wt * low + (1 - low_wt) * upp
    if display:
        print("------------------------------------------------")
        print(" step [  low       upp   ] (  new      f(new)  )")
        print("------------------------------------------------")
    # Loop until the movement is less than "prec".
    step = 1
    while (min(new-low, upp-new) > prec) and (step <= max_steps):
        # Compute the new function value.
        fnew = f(new)
        if display: print(" %3i  [%.2e  %.2e] (%.2e  %.2e)"%(step, low, upp, new, fnew))
        # Perform appropriate bound substitution based on function value.
        if  (fnew == 0): return new
        elif (fnew < 0): low, flow = new, abs(fnew)
        elif (fnew > 0): upp, fupp = new, abs(fnew)
        # Compute the new guess for the zero by linearly interpolating.
        low_wt = fupp / (flow + fupp)
        new = low_wt * low + (1 - low_wt) * upp
        # Increment the step (for unexpected failure cases).
        step += 1
    if display: print("------------------------------------------------")
    return new
# ====================================================================




if __name__ == "__main__":
    # Test function for the different searches. Given a search
    # function and an objective function generator (given a solution).
    def _test_serach(search, f_gen, print=print):
        import random
        random.seed(0)

        # Pick which function we want to test.
        print()
        print(f"Testing '{str(search).split()[1]}'..")
        accuracy = 2**(-26)

        # Test `binary_search` without any modification.
        solution = random.random() + 1/2
        f = f_gen(solution)
        lower = 0
        upper = 2
        print("f(lower): ",f(lower))
        print("f(upper): ",f(upper))

        # Store function evaluations.
        function_values = []
        def g(x):
            function_values.append(f(x))
            return function_values[-1]
        # Perform the binary search.
        output = search(g, lower, upper, accuracy=accuracy)
        print()
        print("solution:    ",float(solution))
        print("output:      ",float(output))
        print("evaluations: ",len(function_values))
        print("accuracy:    ",float(accuracy))
        print("error:       ",float(abs(output - solution)))
        assert(abs(output - solution) < accuracy)

        # Test `search` with rounded outputs.
        round_func = lambda x: round(x, 1)
        solution = round_func(100 * (random.random() + 1/2))
        f = f_gen(solution)
        lower = 0
        upper = 200
        # Reset the list of function values.
        function_values = []
        def g(x):
            function_values.append(f(x))
            return function_values[-1]
        # Perform the binary search.
        output = search(g, lower, upper, accuracy=accuracy, round=round_func)
        print()
        print("solution:    ",float(solution))
        print("output:      ",float(output))
        print("evaluations: ", len(function_values))
        print("accuracy:    ", float(accuracy))
        print("error:       ",float(output - solution))
        assert(abs(output - solution) < accuracy)



    # Run the testing code.
    binary_f_gen = lambda solution: lambda x: x <= solution
    golden_f_gen = lambda solution: lambda x: abs(x - solution)
    _test_serach(binary_search, binary_f_gen)
    _test_serach(golden_section_search, golden_f_gen)
