# Makes available all custom algorithms that I have worked with
from util.approximate import Approximator

# Use keras and tensorflow to make a multi-layer perceptron regressor.
class NeuralNetwork(Approximator):
    def __init__(self, *args, **kwargs): pass
    def _fit(self, x, y, *args, layers=(32,)*10, activation='relu',
             optimizer='sgd', loss='mse', epochs=5000,
             batch_size=None, verbose=0, **kwargs):
        import numpy as np
        # Silence the "backend" message from keras on import.
        import os, sys
        stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
        from keras.models import Sequential
        from keras.layers import Dense, Activation
        sys.stderr = stderr
        # Normalize the data points.
        self.x_shift = np.mean(x, axis=0)
        x -= self.x_shift
        self.x_scale = np.sqrt(np.var(x, axis=0))
        self.x_scale[self.x_scale == 0] = 1.0
        x /= self.x_scale
        # Normalize the data values.
        self.y_shift = np.mean(y, axis=0)
        y -= self.y_shift
        self.y_scale = np.sqrt(np.var(y, axis=0))
        self.y_scale[self.y_scale == 0] = 1.0
        y /= self.y_scale
        # Construct the model layers, starting with the input layer.
        model_layers = [Dense(layers[0], input_shape=(x.shape[1],)), Activation(activation)]
        # Add all internal layers.
        for size in layers[1:]:
            model_layers += [Dense(size), Activation(activation)]
        # Add output layer.
        model_layers += [Dense(y.shape[1])]
        # Construct and fit the model.
        self.mlp = Sequential(model_layers)
        # Use mean squared error, solve with sgd.
        self.mlp.compile(optimizer=optimizer, loss=loss)
        if type(batch_size) == type(None): batch_size=min(200, len(x))
        # Turn off the tensorflow warning messages (about sub-optimality).
        #   0 = all messages are logged (default behavior)
        #   1 = INFO messages are not printed
        #   2 = INFO and WARNING messages are not printed
        #   3 = INFO, WARNING, and ERROR messages are not printed
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
        # Fit the model.
        self.mlp.fit(x, y, epochs=epochs, batch_size=batch_size, verbose=verbose)

    def _predict(self, x, *args, **kwargs):
        # Normalize the input data and denormalize the output prediction.
        x = (x - self.x_shift) / self.x_scale
        y = self.mlp.predict(x, *args, **kwargs)
        return (y * self.y_scale) + self.y_shift


if __name__ == "__main__":
    from util.approximate.testing import test_plot
    m = NeuralNetwork()
    p, x, y = test_plot(m, random=True, N=200)
    p.show()
