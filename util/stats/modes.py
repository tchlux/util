from util.stats import *

# ============================================================
#      Automatic detection of modes (based on confidence)     
# ============================================================

def modes(data, confidence=.99, tol=1/1000):
    from util.optimize import zero
    num_samples = len(data)
    error = 2*samples(num_samples, confidence=confidence)
    print()
    print("Smallest allowed mode: ",error)
    print()
    # Get the CDF points (known to be true based on data).
    x, y = cdf_points(data)
    cdf = cdf_fit((x,y), fit="linear")
    x, y = x[1:], y[1:]
    # Generate the candidate break points based on linear interpolation.
    slopes = [(y[i+1] - y[i]) / (x[i+1] - x[i]) for i in range(len(x)-1)]
    candidates = [i for i in range(1,len(slopes)-1)
                  if (slopes[i] < slopes[i-1]) and (slopes[i] < slopes[i+1])]
    # Sort candidates by their 'width', the widest being the obvious divisors.
    candidates = sorted(candidates, key=lambda i: -(x[i+1] - x[i]))
    print("candidates: ",candidates)
    print("slopes:     ",[slopes[c] for c in candidates])
    # Break the data at candidates as much as can be done with confidence.
    breaks = [min(x), max(x)]
    sizes = [1.]
    chosen = []
    print()
    print("breaks: ",breaks)
    print("sizes:  ",sizes)
    print("chosen: ",chosen)
    print()
    # Loop until there are no candidate break points left.
    while len(candidates) > 0:
        new_break_idx = candidates.pop(0)
        new_break = (x[new_break_idx + 1] + x[new_break_idx]) / 2
        b_idx = np.searchsorted(breaks, new_break, side="right")
        # Compute the CDF values at the upper, break, and lower positions.
        upp = cdf(breaks[b_idx+1]) if (b_idx < len(breaks)-1) else cdf.max
        mid = cdf(new_break)
        low = cdf(breaks[b_idx-1])
        print()
        print("new_break: ", new_break)
        print("b_idx: ",b_idx)
        print("  upp: ",upp, upp - mid)
        print("  mid: ",mid)
        print("  low: ",low, mid - low)
        # Compute the size of the smallest mode resulting from the break.
        smallest_result = min(upp - mid, mid - low)
        # Skip the break if it makes a mode smaller than error.
        if smallest_result < error: continue

        print()
        print("Num_modes: ", len(sizes))
        print("breaks:    ", breaks)
        print("sizes:     ", sizes)
        print()

        # Update the "breaks" and "sizes" lists.
        breaks.insert(b_idx, new_break)
        sizes.insert(b_idx, upp - mid)
        sizes[b_idx-1] = mid - low
        chosen.append(new_break_idx)
        
    # From the "breaks" and "sizes", construct a list of "modes".
    # Consider the most dense point between breaks the "mode".
    modes = []
    for i in range(len(chosen)):
        low = 0 if (i == 0) else chosen[i-1]
        upp = len(slopes)-1 if (i == len(chosen)-1) else chosen[i+1]
        mode_idx = low + np.argmax(slopes[low:upp])
        modes.append( (x[mode_idx+1] + x[mode_idx]) / 2 )


    from util.plot import Plot
    p = Plot()
    pdf = pdf_fit(data)
    p.add_func("PDF", pdf, pdf(), color=p.color(1))
    # Generate the mode lines.
    mode_lines = [[],[]]
    for z in modes:
        mode_lines[0] += [z,z,None]
        mode_lines[1] += [0,.2,None]
    p.add("modes", *mode_lines, color=p.color(0), mode="lines", group="modes")
    # Generate the antimode lines.
    break_lines = [[],[]]
    for z in breaks:
        break_lines[0] += [z,z,None]
        break_lines[1] += [0,.2,None]
    p.add("seperator", *break_lines, color=p.color(3,alpha=.3), mode="lines", group="seperator")
    p.show()
    # Show CDF
    p = Plot()
    pdf = pdf_fit(data)
    p.add_func("CDF", cdf, cdf(), color=p.color(1))
    # Generate the mode lines.
    mode_lines = [[],[]]
    for z in modes:
        mode_lines[0] += [z,z,None]
        mode_lines[1] += [0,1,None]
    p.add("modes", *mode_lines, color=p.color(0), mode="lines", group="modes")
    # Generate the antimode lines.
    break_lines = [[],[]]
    for z in breaks:
        break_lines[0] += [z,z,None]
        break_lines[1] += [0,1,None]
    p.add("seperator", *break_lines, color=p.color(3,alpha=.3), mode="lines", group="seperator")
    p.show(append=True)

    # Try and break the data up at the lowest density regions,
    # continue breaking the data until we cannot break it any more
    # based on confidence.

