from util.approximate.base import Approximator

# Support Vector Regressor
SVR_KERNEL = "poly"


# ==================================
#      Support Vector Regressor     
# ==================================
class SVR(Approximator):
    def __init__(self, *args, kernel=SVR_KERNEL, **kwargs):
        # Disable the future warnings that SVR gives.
        import warnings
        warnings.filterwarnings("ignore", category=FutureWarning)
        # Import the model from sklearn.
        from sklearn.svm import SVR
        self.SVR = SVR
        kwargs.update(dict(kernel=kernel))
        self.svr = self.SVR(*args, **kwargs)
    def _fit(self, x, y, *args, **kwargs):
        if (y.shape[1] > 1): raise(Exception("SVR only supports 1D output."))
        y = y[:,0]
        return self.svr.fit(x, y, *args, **kwargs)
    def _predict(self, *args, **kwargs):
        output = self.svr.predict(*args, **kwargs)
        if (len(output.shape) == 1): output = output[:,None]
        return output

class SVR_RBF(Approximator):
    def __init__(self, *args, kernel="rbf", **kwargs):
        # Disable the future warnings that SVR gives.
        import warnings
        warnings.filterwarnings("ignore", category=FutureWarning)
        # Import the model from sklearn.
        from sklearn.svm import SVR
        self.SVR = SVR
        kwargs.update(dict(kernel=kernel))
        self.svr = self.SVR(*args, **kwargs)
    def _fit(self, x, y, *args, **kwargs):
        if (y.shape[1] > 1): raise(Exception("SVR only supports 1D output."))
        y = y[:,0]
        return self.svr.fit(x, y, *args, **kwargs)
    def _predict(self, *args, **kwargs):
        output = self.svr.predict(*args, **kwargs)
        if (len(output.shape) == 1): output = output[:,None]
        return output


if __name__ == "__main__":
    from util.approximate.testing import test_plot
    m = SVR()
    p, x, y = test_plot(m, random=True, N=200)
    p.show()
