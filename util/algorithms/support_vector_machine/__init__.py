# Support Vector Regressor
SVR_KERNEL = "poly"


    
# ==================================
#      Support Vector Regressor     
# ==================================
class SVR(Approximator):
    def __init__(self, *args, kernel=SVR_KERNEL, **kwargs):
        from sklearn.svm import SVR
        self.SVR = SVR
        kwargs.update(dict(kernel=kernel))
        self.svr = self.SVR(*args, **kwargs)
    def fit(self, *args, **kwargs):
        return self.svr.fit(*args, **kwargs)
    def predict(self, *args, **kwargs):
        return self.svr.predict(*args, **kwargs)
