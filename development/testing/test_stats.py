from util.stats import *

# ====================================================================
#                           util.stats     
# ====================================================================

def _test_mpca(display=False):
    if display:
        print()
        print("-"*70)
        print("Begin tests for 'mpca'")

    GENERATE_APPROXIMATIONS = True
    BIG_PLOTS = False
    SHOW_2D_POINTS = False

    import random
    # Generate some points for testing.
    np.random.seed(4) # 4
    random.seed(0) # 0
    rgen = np.random.RandomState(1) # 10
    n = 100
    points = (rgen.rand(n,2) - .5) * 2
    # points *= np.array([.5, 1.])

    # Create some testing functions (for learning different behaviors)
    funcs = [
        lambda x: x[1]               , # Linear on y
        lambda x: abs(x[0] + x[1])   , # "V" function on 1:1 diagonal
        lambda x: abs(2*x[0] + x[1]) , # "V" function on 2:1 diagonal
        lambda x: x[0]**2            , # Quadratic on x
        lambda x: (x[0] + x[1])**2   , # Quadratic on 1:1 diagonal
        lambda x: (2*x[0] + x[1])**3 , # Cubic on 2:1 diagonal
        lambda x: (x[0]**3)          , # Cubic on x
        lambda x: rgen.rand()        , # Random function
    ]
    # Calculate the response values associated with each function.
    responses = np.vstack(tuple(tuple(map(f, points)) for f in funcs)).T

    # Reduce to just the first function
    choice = 3
    func = funcs[choice]
    response = responses[:,choice]

    # Run the princinple response analysis function.
    components, values = mpca(points, response)
    values /= np.sum(values)
    conditioner = np.matmul(components, np.diag(values))

    if display:
        print()
        print("Components")
        print(components)
        print()
        print("Values")
        print(values)
        print()
        print("Conditioner")
        print(conditioner)
        print()

    components = np.array([[1.0, 0.], [0., 1.]])
    values = normalize_error(np.matmul(points, components.T), response, abs_diff)
    values /= np.sum(values)
    if display:
        print()
        print()
        print("True Components")
        print(components)
        print()
        print("True Values")
        print(values)
        print()


    # Generate a plot of the response surfaces.
    from util.plot import Plot, multiplot
    if display: print("Generating plots of source function..")

    # Add function 1
    p1 = Plot()
    p1.add("Points", *(points.T), response, opacity=.8)
    p1.add_func("Surface", func, [-1,1], [-1,1], plot_points=100)

    if GENERATE_APPROXIMATIONS:
        from util.approximate import NearestNeighbor, Delaunay, condition
        p = Plot()
        # Add the source points and a Delaunay fit.
        p.add("Points", *(points.T), response, opacity=.8)
        p.add_func("Truth", func, [-1,1], [-1,1])
        # Add an unconditioned nearest neighbor fit.
        model = NearestNeighbor()
        model.fit(points, response)
        p.add_func("Unconditioned Approximation", model, [-1,1], [-1,1],
                    mode="markers", opacity=.8)
        # Generate a conditioned approximation
        model = condition(NearestNeighbor, method="MPCA")()
        model.fit(points, response)
        p.add_func("Best Approximation", model, [-1,1], [-1,1],
                    mode="markers", opacity=.8)

        if display: p.plot(show=False, height=400, width=650)

    if display: print("Generating metric principle components..")

    # Return the between vectors and the differences between those points.
    def between(x, y, unique=True):
        vecs = []
        diffs = []
        for i1 in range(x.shape[0]):
            start = i1+1 if unique else 0
            for i2 in range(start, x.shape[0]):
                if (i1 == i2): continue
                vecs.append(x[i2] - x[i1])
                diffs.append(y[i2] - y[i1])
        return np.array(vecs), np.array(diffs)

    # Plot the between slopes to verify they are working.
    # Calculate the between slopes
    vecs, diffs = between(points, response)
    vec_lengths = np.sqrt(np.sum(vecs**2, axis=1))
    between_slopes = diffs / vec_lengths
    bs = ((vecs.T / vec_lengths) * between_slopes).T
    # Extrac a random subset for display
    size = 100
    random_subset = np.arange(len(bs))
    rgen.shuffle(random_subset)
    bs = bs[random_subset[:size],:]
    # Normalize the between slopes so they fit on the plot
    max_bs_len = np.max(np.sqrt(np.sum(bs**2, axis=1)))
    bs /= max_bs_len
    # Get a random subset of the between slopes and plot them.
    p2 = Plot("","Metric PCA on Z","")
    p2.add("Between Slopes", *(bs.T), color=p2.color(4, alpha=.4))

    if SHOW_2D_POINTS:
        # Add the points and transformed points for demonstration.
        new_pts = np.matmul(np.matmul(conditioner, points), np.linalg.inv(components))
        p2.add("Original Points", *(points.T))
        p2.add("Transformed Points", *(new_pts.T), color=p2.color(6, alpha=.7))

    # Add the principle response components 
    for i,(vec,m) in enumerate(zip(components, values)):
        vec = vec * m
        p2.add(f"PC {i+1}", [0,vec[0]], [0,vec[1]], mode="lines")
        ax, ay = (vec / sum(vec**2)**.5) * 3
        p2.add_annotation(f"{m:.2f}", vec[0], vec[1])


    p3 = Plot("", "PCA on X", "")
    p3.add("Points", *(points.T), color=p3.color(4, alpha=.4))

    # Add the normal principle components
    components, values = pca(points)
    values /= np.sum(values)
    for i,(vec,m) in enumerate(zip(components, values)):
        vec = vec * m
        p3.add(f"PC {i+1}", [0,vec[0]], [0,vec[1]], mode="lines")
        ax, ay = (vec / sum(vec**2)**.5) * 3
        p3.add_annotation(f"{m:.2f}", vec[0], vec[1])


    if BIG_PLOTS:
        if display: p1.plot(file_name="source_func.html", show=False)
        if display: p2.plot(append=True, x_range=[-8,8], y_range=[-5,5])
    else:
        # Make the plots (with manual ranges)
        p1 = p1.plot(html=False, show_legend=False)
        p2 = p2.plot(html=False, x_range=[-1,1], y_range=[-1,1], show_legend=False)
        p3 = p3.plot(html=False, x_range=[-1,1], y_range=[-1,1], show_legend=False)
        # Generate the multiplot of the two side-by-side figures
        if display: multiplot([p1,p2,p3], height=126, width=650, append=True)
 
    if display: print("-"*70)


def _test_effect(display=False):
    if display:
        print()
        print("-"*70)
        print("Begin tests for 'effect'")

    a = {"a":.4, "b":.1, "c":.5, "d":.0}
    b = {"a":.1, "b":.3, "c":.5, "d":.1}
    c = {"a":.0, "b":.0, "c":.0, "d":1.}
    assert(.3 == categorical_diff(a, b))
    assert(1. == categorical_diff(a, c))
    assert(.9 == categorical_diff(b, c))
    a = ['a','a','a','a','b','b','b','b']
    b = ['a','a','b','b','b','a','a','b']
    c = ['a','a','a','b','b','a','a','a']
    assert(0. == effect(a,b))
    assert(0. == effect(a,c))
    # assert((4/6 + 3/6 + 0/6 + 0/6)/4 == effect(b, c))
    a = list(range(1000))
    b = list(range(1000))
    assert(effect(a,b) == 1.0)
    b = np.random.random(size=(1000,))
    assert(effect(a,b) < .1)
    a = ['a', 'a', 'a', 'b', 'a', 'b', 'a', 'b', 'a', 'b', 'a', 'a', 'a', 'c', 'a', 'c', 'a', 'c', 'a', 'c']
    b = [890.79, 1048.97, 658.43, 659.39, 722.0, 723.34, 1040.76, 1058.02, 1177.0, 1170.94, 415.56, 462.03, 389.09, 676.82, 688.49, 735.56, 552.58, 1127.29, 1146.42, 1251.12]
    assert(0.0 == effect(a,b) - 0.5264043528589452)
    if display: print("-"*70)


def _test_epdf_diff(display=False):
    if display:
        print()
        print("-"*70)
        print("Begin tests for 'epdf_diff'")

    # ----------------------------------------------------------------
    def demo(seq):
        if display:
            print('~'*70)
            print(len(seq), seq)
            print()
        total = 0
        for vals in edf_pair_gen(seq):
            total += vals[-1]
            if display: print("[% 4s, % 3s] (%.2f) --"%vals, round(total,3))
        if display:
            print('~'*70)
            print()
    demo( [0] + list(range(9)) )
    demo( sorted(np.random.random(size=(10,))) )
    demo( list(range(9)) + [8] )
    # ----------------------------------------------------------------
    # a = [1, 1, 3, 3, 5, 6]
    # b = [0, 1, 2, 3, 4]
    # 
    n = 100
    if display:
        print(f"a = [0,100] ({n} samples)")
        print(f"b = [v + d for v in a]")
        print()
    for d in (.0, .01, .1, .5, .55, .9, 1., 1.5):
        a = [v / n for v in list(range(n+1))]
        b = [v+d for v in a]
        if display: print(f"d = {d:.2f}   a~b = {epdf_diff(a,b):.2f}   b~a = {epdf_diff(b,a):.2f}   a~a = {epdf_diff(a,a):.2f}   b~b = {epdf_diff(b,b):.2f}")
    if display: print()

    for d in (.0, .01, .1, .5, .55, .9, 1., 1.5):
        # Generate a random sequence.
        a = sorted((np.random.random(size=(10,))))
        b = sorted((np.random.random(size=(1000,)) + d))
        diff = epdf_diff(a, b)
        if display:
            print(f"d = {d:.2f}","",
                  "[%.2f, %.2f]"%(min(a),max(a)), "[%.2f, %.2f]"%(min(b),max(b)),"",
                  diff
            )
    if display: print()
    # ----------------------------------------------------------------
    from util.plot import Plot

    # Generate a random sequence.
    a = sorted((np.random.random(size=(2000,))))
    b = sorted((np.random.random(size=(2000,))))

    p = Plot("Empirical PDF Diff Test")
    p1 = pdf_fit(a, smooth=0.00001)
    p2 = pdf_fit(b, smooth=0.00001)
    p.add_func("a", p1, p1()) #(-.5,2))
    p.add_func("b", p2, p2()) #(-.5,2))
    if display: p.show(show=False, y_range=[-.5,1.5])

    p = Plot("Empirical CDF Diff Test")
    p1 = cdf_fit(a)
    p2 = cdf_fit(b)
    p.add_func("a", p1, p1()) #(-.5,2))
    p.add_func("b", p2, p2()) #(-.5,2))
    if display: p.show(append=True)
    # ----------------------------------------------------------------
    if display: print("-"*70)


def _test_fit_funcs(display=False):
    if display:
        print()
        print("-"*70)
        print("Begin tests for 'fit_funcs'")

    from util.plot import Plot

    # ==============================================
    #      Test the fit functions and smoothing     
    # ==============================================
    # Make data
    smooth = .1
    data = np.random.normal(size=(1000,))
    # data[:len(data)//2] += 2
    min_max = (min(data) - .1, max(data) + .1)
    if display:
        print()
        print("(min, max) : (%.2f, %.2f)"%(min_max))
        print("Normal confidence: %.2f%%"%(100*normal_confidence(data)))
        print()
    # Make PDF fits
    pfit = pdf_fit(data)
    smooth_pfit = pdf_fit(data, smooth=smooth)
    # Make CDF fits
    cfit = cdf_fit(data)
    stdev = .05 * (cfit.max - cfit.min)
    smooth_cfit = gauss_smooth(cfit, stdev)
    stdev = smooth * (cfit.max - cfit.min)
    smoother_cfit = gauss_smooth(cfit, stdev)

    # Make PDF plots
    p = Plot()
    p.add_func("PDF", pfit, min_max)
    # Make smooth PDF
    p.add_func("Smooth PDF", smooth_pfit, min_max)
    if display: p.show(show=False)
    # Make CDF plots
    p = Plot()
    p.add_func("CDF", cfit, min_max)
    # Make CDF whose derivative is the default PDF.
    p.add_func("CDF for default PDF", smooth_cfit, min_max)
    # Make smoother cdf.
    p.add_func("Smooth CDF", smoother_cfit, min_max)
    if display: p.show(append=True)
    # Make an animation transitioning between two normal distributions.
    np.random.seed(0)
    d1 = np.random.normal(0, .5, size=(500,))
    d2 = np.random.normal(3, 1, size=(500,))
    f1 = cdf_fit(d1, smooth=.1)
    f2 = cdf_fit(d2, smooth=.1)
    p = Plot()
    for w in np.linspace(0,1,21):
        w = round(w,2)
        f3 = w*f1 + (1-w)*f2
        p.add_func("0, 1/2", f1, f1(), frame=w)
        p.add_func("3, 1", f2, f2(), frame=w)
        p.add_func("weighted sum", f3, f3(), frame=w)
    if display: p.show(bounce=True, append=True)
    if display: print()
    if display: print("-"*70)


def _test_samples(display=True, test_correctness=False):
    from util.math import Fraction
    from util.plot import Plot
    for size in tuple(range(2,42))+(128, 129, 256, 257):
        p = Plot(f"Error at x with {size} samples")
        for confidence in (Fraction(9,10), Fraction(185,200),
                           Fraction(95,100), Fraction(97,100),
                           Fraction(99,100)):
            f = lambda x: samples(size=size, confidence=confidence, at=x[0])
            p.add_func(f"{confidence} confidence", f, [0, 1])
        p.show(append=True, show=(size==2))
    exit()



    if display:
        print()
        print("-"*70)
        print("Begin tests for 'samples'")
        print()
    key_values = [
        (11, Fraction(3,10), Fraction(95,100)),
        (25, Fraction(2,10), Fraction(95,100)),
        (97, Fraction(1,10), Fraction(95,100)),
        (385, Fraction(5,100), Fraction(95,100)),
        (9604, Fraction(1,100), Fraction(95,100)),

        (19, Fraction(3,10), Fraction(99,100)),
        (42, Fraction(2,10), Fraction(99,100)),
        (166, Fraction(1,10), Fraction(99,100)),
        (664, Fraction(5,100), Fraction(99,100)),
        (16588, Fraction(1,100), Fraction(99,100)),

        (33733, Fraction(1,10), Fraction(999,1000)),
        (134930, Fraction(5,100), Fraction(999,1000)),
        (3373242, Fraction(1,100), Fraction(999,1000)),
    ]
    if display: print("samples (max error, confidence)")
    for (s, e,c) in key_values[:-3]:
        needed = samples(error=e, confidence=c)
        print("needed: ",needed)
        if display: print("%6d  (%2.0f%%, %2.0f%%)"%(needed,100*e,100*c))
        # if (s != None): assert(needed == s)

    if display:
        print()
        for n in (10, 25, 90, 350, 9000):
            print(f"With {n:4d} samples we are 95% confident in max CDF error <=",
                  round(samples(n, confidence=.95), 1 if n < 350 else 2))
        print()
        for n in (20, 40, 160, 660, 16000):
            print(f"With {n:5d} samples we are 99% confident in max CDF error <=",
                  round(samples(n, confidence=.99), 1 if n < 600 else 2))
        print("-"*70)

    if test_correctness:
        # Generate a random CDF
        from util import random
        from util.plot import Plot

        TESTS = 10000
        DIFFS = 100
        N = 100

        fit = "linear"
        truth = random.cdf(nodes=5, fit=fit)
        max_error = samples(N, confidence=.99)
        print("Largest expected error:", max_error)
        mid_error = []
        errors = {(1/100):[], (1/4):[], (1/3):[], (1/2):[]}
        max_errors = []
        mean_failed = []
        for i in range(TESTS):
            sample = truth.inverse(np.random.random((N,)))
            guess = cdf_fit(sample, fit=fit)
            diff = truth - guess
            diff_func = lambda x: abs(truth(x) - guess(x))
            diffs = diff_func(np.linspace(0,1,DIFFS))
            mean_failed += [sum(diffs > max_error)]
            print(f"Failed: {mean_failed[-1]:4d}   {sum(mean_failed)/len(mean_failed):.0f}")
            max_errors.append(diff)
            for v in errors: errors[v].append(truth(v) - guess(v))
            # if (diff > max_error):
            #     print(i, N, diff, max_error)
            #     p = Plot()
            #     p.add_func("Truth", truth, truth())
            #     p.add_func("Guess", guess, guess())
            #     p.add_func("Error", diff_func, truth())
            #     p.show(show=False)
            #     p = Plot()
            #     p.add_func("Truth Inverse", truth.inverse, (0.,1.))
            #     p.add_func("Guess Inverse", guess.inverse, (0.,1.))
            #     p.show(show=False, append=True)
            #     break

        total_failed = sum(e > max_error for e in max_errors)
        print(f"Failed {total_failed} out of {TESTS}, or {100*total_failed/TESTS:.1f}%.")


        p = Plot()
        # Add the distribution of maximum errors.
        f = cdf_fit(max_errors)
        p.add_func(f"{len(max_errors)} max errors", f, f())

        # Add the distribution of errors at different values.
        for v in sorted(errors)[::-1]:
            mean = np.mean(errors[v])
            std = np.std(errors[v])
            f = cdf_fit(errors[v])
            p.add_func(f"{v:.1f} errors ({mean:.1e}, {std:.1e}) {samples(N,confidence=.99):.2f}", f, f())

        p.show(append=True)

        p = Plot()
        p.add_func("Truth", truth, truth())
        p.show(append=True)



def _test_Distribution(display=False):
    # Verify that the distribution works under a weighted sum.
    import numpy as np
    d = []
    count = 3
    scale = 2**(-20)
    scale = 2**(30)
    for i in range(count):
        pts = np.random.random(20) * np.random.random() * scale
        d.append( cdf_fit(pts, fit="cubic") )
    wts = np.random.random((count,))
    wts /= sum(wts)
    min_max = (-scale / 3, scale + scale/3)
    if display:
        print(sum(dist*w for (dist,w) in zip(d, wts)))
        from util.plot import Plot
        p = Plot("Weighted cubic fit")
        out = sum(dist*w for (dist,w) in zip(d, wts))
        p.add_func("Weighted sum", out, min_max)
        for i,(dist,w) in enumerate(zip(d, wts)):
            p.add_func(f"Dist {i+1} -- {round(w,3)}", dist, min_max, opacity=.3)
        p.show()

def _test_cdf_fit():
    import numpy as np
    from util.math import SMALL
    from util.random import cdf
    from util.plot import Plot


    n = 10000
    f = cdf(nodes=5, fit="linear")

    sample = f.inverse(np.random.random((n,)))
    g = cdf_fit(sample, fit=None)
    print("Fit first CDF point: ", g.nodes[0])
    print("Fit second CDF point:", g.nodes[1])
    print("Fit last CDF point:  ", g.nodes[-1])

    print("Expected max error:", samples(n, confidence=.99))
    print("Actual max error:  ", f - g)

    min_max = (f.min - SMALL, f.max+SMALL)

    p = Plot()
    p.add_func("Truth", f, min_max)
    p.add_func("EDF", g, min_max)
    p.show(show=False, height=700, width=800)

    p = Plot()
    p.add_func("Truth", f.inverse, min_max)
    p.add_func("EDF", g.inverse, min_max)
    p.show(append=True, height=700, width=800)

    


def test():
    print(f"Testing 'util.stats'..")
    # _test_mpca()
    # _test_effect()
    # _test_epdf_diff()
    # _test_fit_funcs()
    # _test_Distribution()
    _test_samples(True)
    print("done.")
    
if __name__ == "__main__":
    test()
