
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


GOLDEN_RATIO = (5**(1/2) + 1) / 2

# Perform a golden section search on a line over a convex 1D objective
# function. Find the minimum point of the function.
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
    # Compute the initial two golden ratio midpoints.
    d = (u - l) / GOLDEN_RATIO
    m1 = u - d
    m2 = l + d
    fm1 = f(m1)
    fm2 = f(m2)
    # Perform the golden section search.
    while (abs(u - l) > accuracy):
        # Based on the midpoint evaluations, decide which side to search.
        if (fm1 <= fm2):
            # In this case, the minimum must be on the left side.
            l, m2, u = l, m1, m2
            d = (u - l) / GOLDEN_RATIO            
            m1 = round(u - d)
            # If we have reached the limit of achievable accuracy, return.
            if m1 in {l,m2,u}: return m1
            # Otherwise compute the new value and move forward.
            fm1 = f(m1)
        else:
            # In this case, the minimum must be on the right side.
            l, m1, u = m1, m2, u
            d = (u - l) / GOLDEN_RATIO            
            m2 = round(u + d)
            # If we have reached the limit of achievable accuracy, return.
            if m2 in {l,m1,u}: return m2
            # Otherwise compute the new value and move forward.
            fm2 = f(m2)
    # Return the midpoint of the upper and lower bounds.
    return round(l/2 + u/2)





if __name__ == "__main__":
    def _binary_serach(print=print):
        import random
        random.seed(0)

        # Pick which function we want to test.
        print()
        print(f"Testing '{str(binary_search).split()[1]}'..")

        # Test `binary_search` without any modification.
        solution = random.random() + 1/2
        f = lambda x: x <= solution
        accuracy = 2**(-26)
        lower = 0
        upper = 2
        # Create a cached version of `f` to save excessive function evaluations.
        function_values = {}
        def g(x):
            if (x not in function_values): function_values[x] = f(x)
            return function_values[x]
        # Perform the binary search.
        output = binary_search(g, lower, upper, accuracy)
        print()
        print("solution:    ",solution)
        print("output:      ",output)
        print("evaluations: ",len(function_values))
        print("error:       ",abs(output - solution))
        assert(abs(output - solution) < accuracy)

        # Test `binary_search` with rounded outputs.
        round_func = lambda x: round(x, 1)
        solution = round_func(100 * (random.random() + 1/2))
        f = lambda x: x <= solution
        accuracy = 2**(-26)
        lower = 0
        upper = 200
        # Create a cached version of `f` to save excessive function evaluations.
        function_values = {}
        def g(x):
            if (x not in function_values): function_values[x] = f(x)
            return function_values[x]
        # Perform the binary search.
        output = binary_search(g, 0, 200, round=round_func)
        print()
        print("solution:    ",solution)
        print("output:      ",output)
        print("evaluations: ", len(function_values))
        print("error:       ",output - solution)
        assert(abs(output - solution) < accuracy)

    # Run the testing code.
    none_function = lambda *args, **kwargs: None
    _binary_serach()
