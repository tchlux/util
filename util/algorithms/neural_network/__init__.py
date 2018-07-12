# Makes available all custom algorithms that I have worked with
from util.algorithms import Approximator

# MLP Regressor
MLP_R_ACTIVATION_FUNC = "relu"
MLP_R_SOLVER = "sgd"
MLP_ZERO_MEAN_UNIT_VAR = True

# Force many iterations of SGD with a chosen random state and negative tolerance.
# model = MLP(hidden_layer_sizes=(hidden_nodes,), solver='sgd',
#             max_iter=30000, random_state=10, tol=-float('inf'),
#             verbose=True, shuffle=False)


# ==========================================
#      Multi-Layer Perceptron Regressor     
# ==========================================
class MLPRegressor(Approximator):
    def __init__(self, *args, activation=MLP_R_ACTIVATION_FUNC,
                 solver=MLP_R_SOLVER, zmuv=MLP_ZERO_MEAN_UNIT_VAR, **kwargs):
        from sklearn.neural_network import MLPRegressor
        self.MLPRegressor = MLPRegressor
        kwargs.update(dict(activation=activation, solver=solver))
        self.mlp = self.MLPRegressor(*args, **kwargs)
        self.zmuv = zmuv # Zero mean, unit variance
        self.mean = 0
        self.var = 1
    def _fit(self, *args, **kwargs):
        # By default, zero-mean and unit-variance the data.
        if self.zmuv and (len(args) > 0):
            x = args[0]
            self.mean = x.mean(axis=0)
            self.var = (x.var(axis=0))**(1/2)
            x = (x - self.mean) / self.var
            args = (x,) + args[1:]
        return self.mlp.fit(*args, **kwargs)
    def _predict(self, *args, **kwargs):
        # By default, zero-mean and unit-variance the data.
        if len(args) > 0:
            x = args[0]
            x = x - self.mean
            x = x / self.var
            args = (x,) + args[1:]
        output = self.mlp.predict(*args, **kwargs)
        if len(output.shape) == 1: output = output[:,None]
        return output

