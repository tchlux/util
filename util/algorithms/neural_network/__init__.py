# MLP Regressor
MLP_R_ACTIVATION_FUNC = "relu"
MLP_R_SOLVER = "lbfgs"

# ==========================================
#      Multi-Layer Perceptron Regressor     
# ==========================================
class MLPRegressor(Approximator):
    def __init__(self, *args, activation=MLP_R_ACTIVATION_FUNC,
                 solver=MLP_R_SOLVER, **kwargs):
        from sklearn.neural_network import MLPRegressor
        self.MLPRegressor = MLPRegressor
        kwargs.update(dict(activation=activation, solver=solver))
        self.mlp = self.MLPRegressor(*args, **kwargs)
    def fit(self, *args, **kwargs):
        return self.mlp.fit(*args, **kwargs)
    def predict(self, *args, **kwargs):
        return self.mlp.predict(*args, **kwargs)

