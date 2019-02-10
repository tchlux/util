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
