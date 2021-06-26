# Makes available all custom algorithms that I have worked with
from util.approximate import Approximator

# Use keras and tensorflow to make a multi-layer perceptron regressor.
# 
#  Fit / Init Parameters:
# 
#    layers      (32,)*10    -- Dense layer structure of this neural network.
#    activation  'relu'      -- Activation function for each layer.
#    optimizer   'sgd'       -- Method of minimizing 'loss' function.
#    loss        'mse'       -- The function measuring error.
#    epochs      500         -- Number of epochs of optimizer iterations to take.
#    batch_size  None        -- The size of each batch for loss minimization.
#    verbose     0           -- The level of Keras verbosity.
#    log_dir     '~/.tf_log' -- The log directory for training.
# 
class NeuralNetwork(Approximator):
    def __init__(self, *args, **kwargs):
        import os
        from tensorflow.keras.optimizers import SGD
        self.optimizers = {
            'sgd' : SGD(learning_rate=0.1, momentum=0.9)
        }
        # Set some default parameters.
        home_dir = os.path.expanduser('~')
        log_dir = os.path.join(home_dir, '.tf_log')
        defaults = dict(layers=(64,)*8, activation='relu',
                        optimizer='adam', loss='mse', epochs=500,
                        batch_size=None, log_dir=log_dir, verbose=0,
                        show=False, workers=None)
        # Update the defaults with the user specs.
        defaults.update(kwargs)
        # Store the default settings for this neural network for fit-time.
        self.default_settings = defaults

    def _fit(self, x, y, log='2', *args, **user_kwargs):
        # Import numpy module.
        import numpy as np
        # Set the backend for Keras to be Tensorflow.
        import os
        # Turn off the tensorflow warning messages (about sub-optimality).
        #   0 = all messages are logged (default behavior)
        #   1 = INFO messages are not printed
        #   2 = INFO and WARNING messages are not printed
        #   3 = INFO, WARNING, and ERROR messages are not printed
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = log
        # Set the defaults, then overwrite them with user keyword arguments.
        kwargs = self.default_settings.copy()
        kwargs.update(user_kwargs)
        # Convert the needed keyword arguments into local variables.
        layers = kwargs.pop('layers')
        activation = kwargs.pop('activation')
        loss = kwargs.pop('loss')
        epochs = kwargs.pop('epochs')
        verbose = kwargs.pop('verbose')
        batch_size = kwargs.pop('batch_size')
        log_dir = kwargs.pop('log_dir')
        optimizer = kwargs.pop('optimizer')
        show = kwargs.pop('show')
        workers = kwargs.pop('workers')
        # Retrieve a default optimizer if it exists.
        if (optimizer in self.optimizers):
            optimizers = self.optimizers[optimizer]
        # Go for a larger batch size (for speed), but don't use all
        # data (assume some stochasticity in the problem).
        if type(batch_size) == type(None):
            batch_size = max(1000, (len(x) // 4) + 1)
        # Check to see if this is partially trained..
        if (hasattr(self, 'mlp') and
            (hasattr(self, 'x_shift') and (len(self.x_shift) == x.shape[1])) and
            (hasattr(self, 'y_shift') and (len(self.y_shift) == y.shape[1]))):
            # Rescale the incoming data.
            x -= self.x_shift
            x /= self.x_scale
            y -= self.y_shift
            y /= self.y_scale
        # Otherwise create this model from scratch and start training.
        else:
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
            # --------------------------------------------------------
            # Build the MLP.
            #   Import keras layers for training.
            from tensorflow.keras.models import Sequential
            from tensorflow.keras.layers import Dense, Activation
            from tensorflow.keras.constraints import max_norm
            # Construct the model layers, starting with the input layer.
            model_layers = [Dense(layers[0], input_shape=(x.shape[1],),
                                  activation=activation)]
            # Add all internal layers.
            for size in layers[1:]:
                model_layers += [
                    Dense(size, activation=activation,
                          kernel_constraint=max_norm(size)),
                ]
            # Add output layer.
            model_layers += [Dense(y.shape[1])]
            # Construct and fit the model.
            self.mlp = Sequential(model_layers)
            # --------------------------------------------------------
        # Enable logging of progress for visualization..
        if (show):
            from tensorflow.keras.callbacks import TensorBoard
            num_updates = max(1, epochs // 100)
            import time
            log_dir = os.path.join(log_dir, time.strftime('%Y-%m-%d_%H-%M-%S'))
            tensorboard_callback = dict(callbacks=TensorBoard(
                log_dir=log_dir, histogram_freq=num_updates))
        else: tensorboard_callback = {}
        # Use mean squared error, compile loss function.
        self.mlp.compile(optimizer=optimizer, loss=loss)
        # Fit the model on the new data.
        import time
        fit_start = time.time()
        self.mlp.fit(x, y, epochs=epochs, batch_size=batch_size,
                     verbose=verbose, workers=workers, **tensorboard_callback)
        fit_end = time.time()
        self.fit_time = fit_end - fit_start
        # Kill the tensorboard visual.
        # Open a tensorboard visualizer to show how training went..
        if (show):
            from subprocess import Popen, PIPE
            from threading import Timer
            import sys, webbrowser
            # Open a process that runs the tensorboard server.
            p = Popen([sys.executable, '-m', 'tensorboard.main',
                       '--logdir', os.path.abspath(log_dir)],
                      stderr=PIPE)
            # Create a timeout in case this fails for some reason..
            t = Timer(15, p.kill)
            # Retreive the local web server address from the output.
            try:
                t.start()
                for line in p.stderr:
                    line = str(line, 'utf-8')
                    if ('http' in line):
                        address = line[line.index('http'):]
                        address = address[:address.index('(')-1]
                        break
                else: raise(Exception("Failed to find address."))
                # Open the tensorboard server page.
                webbrowser.open(address)
                time.sleep(5)
            # Cancel the timer.
            finally: t.cancel()
            # Kill the server (hopefully user isn't upset!)
            p.kill()
        

    def _predict(self, x, *args, **kwargs):
        # Normalize the input data and denormalize the output prediction.
        x = (x - self.x_shift) / self.x_scale
        y = self.mlp(x, *args, **kwargs)
        return (y * self.y_scale) + self.y_shift


    # Given a path to save, create a save directory with multiple files.
    def save(self, directory):
        import os
        if (not os.path.exists(directory)):
            os.mkdir(directory)
        elif (not os.path.isdir(directory)):
            raise(IOError("Gave path to existing file, need a directory.\n  '"+directory+"'"))
        # Save the model.
        model_path = os.path.join(directory, 'model')
        self.mlp.save(model_path)
        # Save the scaling factors.
        scale_path = os.path.join(directory, 'scaling_factors.npz')
        import numpy as np
        np.savez_compressed(scale_path, self.x_shift, self.x_scale, self.y_shift, self.y_scale)


    # Given a path to load, load all directory contents.
    def load(self, directory):
        import os
        if (not os.path.exists(directory)):
            raise(IOError("No directry exists.\n  '"+directory+"'"))            
        elif (not os.path.isdir(directory)):
            raise(IOError("Gave path to existing file, need a directry.\n  '"+directory+"'"))
        # Disable info and warning logs.
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
        # Load the model.
        model_path = os.path.join(directory, 'model')
        from tensorflow.keras.models import load_model
        self.mlp = load_model(model_path)
        # Load the scaling factors.
        scale_path = os.path.join(directory, 'scaling_factors.npz')
        import numpy as np
        data = np.load(scale_path)
        self.x_shift, self.x_scale, self.y_shift, self.y_scale = \
            data['arr_0'], data['arr_1'], data['arr_2'], data['arr_3']


if __name__ == '__main__':
    from util.approximate.testing import test_plot
    m = NeuralNetwork()
    print('Building approximation to test data..')
    p, x, y = test_plot(m, random=True, N=200)
    print('  showing result..')
    p.show()


# 2021-01-24 21:08:32
# 
##################################################
# # dropout=0.2, l2_penalty=1e-5,                #
#                                                #
# # l2_penalty = kwargs.pop('l2_penalty')        #
# # dropout = kwargs.pop('dropout')              #
#                                                #
# # from tensorflow.keras.regularizers import l2 #
#       # kernel_regularizer=l2(l2_penalty)),    #
# # Dropout(dropout),                            #
##################################################
