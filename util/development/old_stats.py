def modes(data, confidence=.99, tol=1/1000):
    from util.optimize import zero
    num_samples = len(data)
    error = 2*samples(num_samples, confidence=confidence)
    cdf = cdf_fit(data, fit="cubic")
    print("error: ",error)
    # Find all of the zeros of the derivative (mode centers / dividers)
    checks = np.linspace(cdf.min, cdf.max, np.ceil(1/tol))
    second_deriv = cdf.derivative.derivative
    deriv_evals = second_deriv(checks)
    modes = [i for i in range(1, len(deriv_evals)) if
             (deriv_evals[i-1] * deriv_evals[i] <= 0) and
             (deriv_evals[i-1] >= deriv_evals[i])]
    antimodes = [i for i in range(1, len(deriv_evals)) if
                 (deriv_evals[i-1] * deriv_evals[i] <= 0) and
                 (deriv_evals[i-1] < deriv_evals[i])]
    # Compute exact modes and antimodes using a zero-finding function.
    modes = [zero(second_deriv, checks[i-1], checks[i]) for i in modes]
    antimodes = [zero(second_deriv, checks[i-1], checks[i]) for i in antimodes]
    original_antimodes = antimodes[:]
    # Fix the bounds of the antimodes to match the distribution.
    if modes[0] < antimodes[0]: antimodes    = [cdf.min] + antimodes
    else:                       antimodes[0] =  cdf.min
    if modes[-1] > antimodes[-1]: antimodes    += [cdf.max]
    else:                         antimodes[-1] =  cdf.max
    # Make sure that there is an antimode between each mode.
    for i in range(len(modes)):
        if antimodes[i] > modes[i]:
            # Update the next antimode with this one (as long as it's not the max).
            if (i < len(modes)-1):
                antimodes[i+1] = (antimodes[i] + antimodes[i+1]) / 2
            # Always update this antimode to properly be LESS than the mode.
            antimodes[i] = (modes[i] + modes[i-1]) / 2
    print("len(modes):     ",len(modes))
    print("len(antimodes): ",len(antimodes))
    # Define a function that counts the number of modes thta are too small.
    def count_too_small():
        return sum( (cdf(upp) - cdf(low)) < error for (low,upp) in
                    zip(antimodes[:-1],antimodes[1:]) )
    # Show PDF
    from util.plot import Plot
    p = Plot()
    pdf = pdf_fit(cdf.inverse(np.random.random((1000,))))

    # Loop until all modes are big enough to be accepted given error tolerance.
    step = 1
    while count_too_small() > 0:
        print()
        print("step: ",step, (len(modes), len(antimodes)))
        f = len(modes)
        p.add_func("PDF", pdf, cdf(), color=p.color(1), frame=f, show_in_legend=(step==1))
        # Generate the mode lines.
        mode_lines = [[],[]]
        for z in modes:
            mode_lines[0] += [z,z,None]
            mode_lines[1] += [0,.2,None]
        p.add("modes", *mode_lines, color=p.color(0), mode="lines",
              group="modes", show_in_legend=(z==modes[0] and step==1), frame=f)
        # Generate the antimode lines.
        anti_lines = [[],[]]
        for z in antimodes:
            anti_lines[0] += [z,z,None]
            anti_lines[1] += [0,.2,None]
        p.add("seperator", *anti_lines, color=p.color(3,alpha=.3), mode="lines",
              group="seperator", show_in_legend=(z==antimodes[0] and (step==1)), frame=f)
        step += 1


        # Compute the densities and the sizes of each mode.
        sizes = [cdf(antimodes[i+1]) - cdf(antimodes[i])
                 for i in range(len(modes))]
        densities = [(cdf(antimodes[i+1]) - cdf(antimodes[i])) /
                     (antimodes[i+1] - antimodes[i])
                     for i in range(len(modes))]
        # Compute those modes that have neighbors that are too small.
        to_grow = [i for i in range(len(modes))
                   if (i > 0 and sizes[i-1] < error)
                   or (i < len(sizes)-1 and sizes[i+1] < error)]
        if len(to_grow) == 0: break
        print("modes:     ",modes)
        print("antimodes: ",antimodes)
        print("sizes:     ",sizes)
        print("densities: ",densities)
        print("to_grow:   ",to_grow)
        # Sort the modes to be grown by their size, largest first.
        to_grow = sorted(to_grow, key=lambda i: -densities[i])
        # Keep track of the modes that have already been absorbed.
        preference = {}
        taken = set()
        conflicts = set()
        modes_to_remove = []
        anti_to_remove  = []
        while len(to_grow) > 0:
            i = to_grow.pop(0)
            # Pick which of the adjacent nodes to absorb.
            to_absorb = None
            if (i < len(modes)-1) and (sizes[i+1] < error):
                direction = 1
                to_absorb = i + 1
            if (i > 0) and (sizes[i-1] < error):
                # If there wasn't a right mode, take the left by default.
                if (to_absorb == None):
                    direction = -1
                    to_absorb = i - 1
                # Otherwise we have to pick based on the density similarity.
                elif (abs(modes[i-1]-modes[i]) < abs(modes[i+1]-modes[i])):
                    # Take the other one if its density is more similar.
                    direction = -1
                    to_absorb = i - 1
            # If there is no good option to absorb, the skip.
            if (to_absorb in preference): continue
            # Record the preferred pick of this mode.
            preference[i] = (direction, to_absorb)
            # If this mode is already absorbed, then add it to conflict list.
            if to_absorb in taken: conflicts.add( to_absorb )
            # Remove the ability to 'absorb' from modes getting absorbed.
            if to_absorb in to_grow: to_grow.remove(to_absorb)
            # Add the absorbed value to the set of "taken" modes.
            taken.add(to_absorb)

        # Resolve conflicts by giving absorbed modes to closer modes.
        for i in sorted(conflicts, key=lambda i: -densities[i]):
            if (abs(modes[i-1] - modes[i]) < abs(modes[i+1] - modes[i])):
                preference.pop(i+1)
            else:
                preference.pop(i-1)

        # Update the boundaries
        for i in sorted(preference, key=lambda i: -densities[i]):
            direction, to_absorb = preference[i]
            # Update the boundary of this mode.
            antimodes[i+(direction>0)] = antimodes[to_absorb + (direction>0)]
            # Update the "to_remove" lists.
            anti_to_remove.append( antimodes[to_absorb + (direction>0)] )
            modes_to_remove.append( modes[to_absorb] )
        # Remove the modes and antimodes that were merged.
        for m in modes_to_remove: modes.remove(m)
        for a in anti_to_remove:  antimodes.remove(a)
        # Update the remaining antimodes to be nearest to the middle
        # of the remaining modes (making them representative dividers).
        for i in range(len(modes)-1):
            middle = (modes[i] + modes[i+1]) / 2
            closest = np.argmin([abs(oam - middle) for oam in original_antimodes])
            antimodes[i+1] = original_antimodes[closest]



    f = len(modes)
    p.add_func("PDF", pdf, cdf(), color=p.color(1), frame=f, show_in_legend=(step==1))
    # Generate the mode lines.
    mode_lines = [[],[]]
    for z in modes:
        mode_lines[0] += [z,z,None]
        mode_lines[1] += [0,.2,None]
    p.add("modes", *mode_lines, color=p.color(0), mode="lines",
          group="modes", show_in_legend=(z==modes[0] and step==1), frame=f)
    # Generate the antimode lines.
    anti_lines = [[],[]]
    for z in antimodes:
        anti_lines[0] += [z,z,None]
        anti_lines[1] += [0,.2,None]
    p.add("seperator", *anti_lines, color=p.color(3,alpha=.3), mode="lines",
          group="seperator", show_in_legend=(z==antimodes[0] and (step==1)), frame=f)

    p.show(append=True, y_range=[0,.15])

    p = Plot()
    p.add_func("CDF", cdf, cdf(), color=p.color(1))
    for z in modes:
        p.add("modes", [z,z], [0,1], color=p.color(0), mode="lines",
              group="modes", show_in_legend=(z==modes[0]))
    for z in antimodes:
        p.add("seperator", [z,z], [0,1], color=p.color(3), mode="lines",
              group="sep", show_in_legend=(z==antimodes[0]))
    p.show(append=True)




# ../development/testing/test_stats.py 
if __name__ == "__main__":
    from util.random import cdf
    np.random.seed(1)
    data = cdf(nodes=3).inverse(np.random.random(100))
    modes(data)

