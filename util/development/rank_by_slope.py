# Originally was a replacement for "rank_by_variation" in "stats.py".

# Given a set of row-vectors, compute a convex weighting that is
# proportional to the inverse average metric slope  between points
# projected onto each vector (component).
def rank_by_slope(components, points, values, metric, max_pairs=10000, display=True):
    # Compute the magnitudes using average metric slope.
    if display: print(" computing average metric slope per component.. ", end="", flush=True)
    avg_slope = np.zeros(len(components))
    update = end = ""
    for i in range(len(components)):
        update = "\b"*len(end)
        end = f"{i+1} of {len(components)}"
        if display: print(update, end=end, flush=True)
        x = np.matmul(points, components[i])
        for (p1, p2) in gen_random_pairs(len(x), count=max_pairs):
            avg_slope[i] += metric(values[p1], values[p2]) / abs(x[p1] - x[p2])
    print("Results so far:")
    from util.plot import Plot
    p = Plot()
    for (comp, slope) in zip(components, avg_slope):
        print(comp, slope)
        x, y = np.matmul(points, comp), values
        p.add(f"Points: {slope:.2f}", x, y)
    print()
    p.show()
    exit()
    # Invert the total average metric slope.
    if (min(avg_slope) <= 0.0):
        avg_slope = np.where(avg_slope == 0, 1., 0.)
    else:
        avg_slope = 1 / avg_slope
    avg_slope /= np.sum(avg_slope)
    # Re-order the components according to inverse average metric slope.
    order = np.argsort(avg_slope)[::-1]
    # If they are not already ordered correctly, re-order the returns.
    if not all(order[i] < order[i+1] for i in range(len(order)-1)):
        if display: print(" reordering components by average metric slope..", end="\r", flush=True)
        components, avg_slope = components[order], avg_slope[order]
    if display: print("                                                ", end="\r", flush=True)
    return components, avg_slope


