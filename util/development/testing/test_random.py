from util.random import *

# ====================================================================
#                           util.random     
# ====================================================================

def test_latin(display=False):
    from util.random import latin
    print("Testing latin..", end=" ")
    if display: print()

    if display:
        from util.plot import Plot
        D = 2
        N = 400

        p = Plot("Latin hyper cube design")
        # Add grid lines for the cube.
        for i in range(N+1):
            p.add(f"gly {i}", [0,1], [i / N]*2, mode="lines", group="grid",
                  color="rgba(.6,.1,.1,.2)", show_in_legend=(i==0))
            p.add(f"glx {i}", [i / N]*2, [0,1], mode="lines", group="grid",
                  color="rgba(.6,.1,.1,.2)", show_in_legend=False)
        # Add the random points.
        p.add("Points", *(latin(N,D).T))
        p.show(file_name="/Users/thomaslux/Desktop/lhc.html")

    from numpy import sort
    N = 10000
    D = 1000
    pts = latin(N, D)
    for i in range(D):
        values = sort(pts[:,i])
        max_adj_diff = max(abs(values[:-1] - values[1:]))
        if max_adj_diff >= (2/N):
            print("BAD!", "   ", i, max_adj_diff, 1/N, 2/N)
        assert(max_adj_diff < (2/N))

    print("passed.")


def test_random_range(display=False):
    # print("Testing random_range..", end=" ")
    # if display: print()
    from util.random import random_range

    if display:
        print()
        print("-"*70)
        print("Begin tests for 'random_range'")

    import random
    # Check unique values witnessed with sliced range.
    seen = {}
    for i in range(1000):
        random.seed(i)
        v = tuple(random_range(4))
        seen[v] = seen.get(v,0) + 1
    if display:
        for row in sorted(seen):
            print("%5d"%seen[row], row)
    # vals = {}
    # mults = {}
    # for i in range(10000):
    #     results = tuple(random_range(20))
    #     key, (val, mult) = results[:-1], results[-1]
    #     diff = set(key[i-1] - key[i] for i in range(len(key)))
    #     if len(diff) <= 2:
    #         vals[key] = vals.get(key, set())   .union({val})
    #         mults[key] = mults.get(key, set()) .union({mult})
    # for v in sorted(vals): print(v, "% 10s"%sorted(vals[v]), sorted(mults[v]))
    # print(len(vals))
    if display: print()
    if display: print("-"*70)


def test_random_cdf(display=True):
    if not display: return

    from util.plot import Plot
    p = Plot("")
    for nodes in range(1, 100):
        f = cdf()
        p.add_func(f"Random PDF {nodes}", f.derivative, f(),
                   color=p.color(nodes-1), group=nodes)
        # p.add(f"Points {nodes}", *list(zip(*f.nodes)),
        #       color=p.color(nodes-1,alpha=.3), group=nodes)
    p.show(show=False)

    print()

    p = Plot("")
    for nodes in range(1, 30):
        f = cdf(nodes=nodes)
        p.add_func(f"Random {nodes} node CDF", f, f(),
                   color=p.color(nodes-1), group=nodes)
        p.add(f"Points {nodes}", *list(zip(*f.nodes)),
              color=p.color(nodes-1,alpha=.3), group=nodes)
    p.show(append=True)


def test():
    test_latin()
    test_random_range()
    test_random_cdf()


if __name__ == "__main__":
    test()
