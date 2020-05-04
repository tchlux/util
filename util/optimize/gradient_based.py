import numpy as np

STEPS = 100000 // 2
SGD_ALPHA = .1
ADA_ALPHA = .01

# Implement ReLU activation network in Fortran, this is a waste.


# L-BFGS minimization.
def L_BFGS(grad, start, budget=1000, m=0, alpha=.99, eps=2.**(-56)):
    from numpy import zeros, dot, sqrt, roll, asarray, array
    points = []
    dim = len(start)
    if (m == 0): m = max(10, int(sqrt(dim)))
    # Initialize storage for internal variables.
    x = start
    g = zeros(dim)
    s = zeros((m, dim))
    y = zeros((m, dim))
    rho = zeros(m)
    a = zeros(m)
    b = zeros(m)
    old_x = zeros(dim)
    old_g = zeros(dim)
    # Loop until the budget is exhausted.
    for i in range(budget):
        points.append(x.copy())
        g = grad(x)
        # Very first step does different stuff for initialization.
        if (i == 0):
            old_x[:] = x[:]
            old_g[:] = g[:]
            x -= alpha * old_g
            continue
        # Normal case.
        # 
        # Roll the history arrays (free the first index) and also
        # check for termination (numerically unstable small step).
        # 
        s = roll(s, 1, axis=0)
        y = roll(y, 1, axis=0)
        rho = roll(rho, 1, axis=0)
        s[0,:] = x[:] - old_x[:]
        y[0,:] = g[:] - old_g[:]
        ys = (dot(y[0], s[0])) # <- Original L-BFGS doesn't "abs",
        #                              but for noisy functions this is
        #                              a necessary change to continue.
        if (sqrt(abs(ys)) < eps): break
        rho[0] = 1 / ys
        # Copy current iterates for next iteraction.
        old_x[:] = x[:]
        old_g[:] = g[:]
        # Reconstruct the BFGS update.
        for j in range(min(m, i)):
            a[j] = rho[j] * dot(s[j], g)
            g -= a[j] * y[j]
        g *= (ys / dot(y[0],y[0]))
        for j in range(min(m, i)-1, -1, -1):
            b[j] = rho[j] * dot(y[j], g)
            g += s[j] * (a[j] - b[j])
        # Compute the rescaled update.
        x -= alpha * g
    else:
        # If the for loop didn't break, then the last point needs to be recorded.
        points.append( x.copy() )
    # Return the set of points.
    return points


# Standard SGD optimization algorithm.
#    grad() - Function that returns the gradient at any point.
#    start  - Contains the starting iterate.
#    budget - Optional argument that sets the number of gradient evaluations.
#    alpha  - Optional argument containing the step size.
#    tau    - Optional argument containing the decay factor.
def SGD(func, grad, start, budget=STEPS, alpha=SGD_ALPHA, tau=0.5):
    points = [start]
    # initialize
    x = start
    # main loop
    for t in range(0,budget):
        # update step
        x = x - alpha * grad(x)
        # decay 10 times over the course of 
        if ((t > 0) and (not t%(budget//10))):
            alpha = max(alpha * tau, np.finfo(np.float64).eps)
        points.append(x)
    return points


# ADAGRAD optimization algorithm as described by Duchi, Hazan, and Singer.
#   grad() - Function that returns the gradient at any point.
#   start  - Contains the starting iterate.
#   budget - Optional argument containing the gradient evaluation budget.
#   alpha  - Optional argument containing the step size in the space induced by G.
#   eps    - Optional argument containing the "Fudge factor" for adjusting G.
def ADAGRAD(func, grad, start, budget=STEPS, alpha=ADA_ALPHA, eps=10.0**(-6.0)):
    points = [start]
    # initialize matrix norm
    x = start
    G = 0
    # main loop
    for t in range(0, budget):
        # get the gradient
        g = grad(x)
        # update the norm
        G = G + (g ** 2.0)
        # take the adjusted trust region step
        g = g / (eps + np.sqrt(G))
        x = x - (alpha * g)
        points.append(x)
    return points


# ADaM optimization algorithm as proposed by D. P. Kingma and J. L. Ba.
#   grad() - Function that returns the gradient at any point.
#   start  - Contains the starting iterate.
#   budget - Optional argument that sets the number of gradient evaluations.
#   alpha  - Optional argument containing the step size.
#   beta1  - Optional argument containing the decay rate for first moment.
#   beta2  - Optional argument containing the decay rate for second moment.
#   eps    - Optional argument containing the fudge factor.
def ADAM(func, grad, start, budget=STEPS, alpha=ADA_ALPHA, beta1=0.9,
         beta2=0.99, eps=10.0**(-8.0)):
    points = [start]
    # initialize
    x = start
    m = 0.0
    v = 0.0
    # main loop
    for t in range(0, budget):
        # get gradient
        g = grad(x)
        # compute biased first and second moments
        m = beta1 * m + (1.0 - beta1) * g
        v = beta2 * v + (1.0 - beta2) * g * g
        # correct for bias
        m_hat = m / (1.0 - beta1 ** float(t+1))
        v_hat = v / (1.0 - beta2 ** float(t+1))
        # update step
        x = x - alpha * m_hat / (np.sqrt(v_hat) + eps)
        points.append(x)
    return points


