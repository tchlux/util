import numpy as np

STEPS = 100000 // 2
SGD_ALPHA = .1
ADA_ALPHA = .01

# BFGS minimization.
def L_BFGS(func, grad, start):
    from scipy.optimize import fmin_l_bfgs_b
    points = [start]
    def wrapped_grad(x):
        points.append(x)
        return grad(x)
    fmin_l_bfgs_b(func, start, fprime=wrapped_grad, maxiter=STEPS)
    return points[:10000]


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


