# A small wrapper model for building and training multi-layered
# perceptrons with Keras and Tensorflow. Allows training by giving
# data to match, or by giving observations and a loss function to
# minimize.
# 
# Automatically handles tensorboard when the 'show=True' option is given.

import os

# Turn off the tensorflow warning messages (about sub-optimality).
#   0 = all messages are logged (default behavior)
#   1 = INFO messages are not printed
#   2 = INFO and WARNING messages are not printed
#   3 = INFO, WARNING, and ERROR messages are not printed
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'


# Transform train_on_batch return value
# to dict expected by on_batch_end callback
def named_losses(model, losses):
    result = {}
    # Make sure the assumptions of this function hold.
    assert hasattr(model, "metrics_names"), "Provided model must have 'metrics_names' attribute."
    # If multiple metrics and losses are given, store them all.
    if (hasattr(losses, "__iter__")):
        assert hasattr(model.metrics_names, "__iter__"), "'model.metrics_names' must be iterable."
        for metric, log in zip(model.metrics_names, losses):
            result[metric] = log
    # If only a single loss is provided, then store it.
    if (issubclass(type(losses),float)):
        result[model.metrics_names[0]] = losses
    # Return the dictionary of metric names and losses.
    return result

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
class NeuralNetwork:
    def __init__(self, *args, **kwargs):
        import os
        # Turn off the tensorflow warning messages (about sub-optimality).
        #   0 = all messages are logged (default behavior)
        #   1 = INFO messages are not printed
        #   2 = INFO and WARNING messages are not printed
        #   3 = INFO, WARNING, and ERROR messages are not printed
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
        # Define default settings for some optimizers.
        from tensorflow.keras.optimizers import SGD
        self.optimizers = {
            'sgd' : SGD(learning_rate=0.1, momentum=0.9)
        }
        # Set some default parameters.
        home_dir = os.path.expanduser('~')
        log_dir = os.path.join(home_dir, '.tf_log')
        defaults = dict(layers=(32,)*10, activation='relu',
                        optimizer='adam', loss='mse', epochs=500,
                        batch_size=None, log_dir=log_dir, log='2',
                        verbose=0, show=False)

        # Update the defaults with the user specs.
        defaults.update(kwargs)
        # Store the default settings for this neural network for fit-time.
        self.default_settings = defaults


    # Build an MLP and calculate shift/scale vectors for inputs to
    #  automatically make them zero mean and unit variance.
    # 
    #  x -> numpy ndarray or equivalent
    #  y -> (numpy ndarray or equivalent) OR (integer, number of outputs)
    # 
    def build_model(self, x, y, *args, **user_kwargs):
        # Set the defaults, then overwrite them with user keyword arguments.
        kwargs = self.default_settings.copy()
        kwargs.update(user_kwargs)
        layers = kwargs.pop('layers')
        activation = kwargs.pop('activation')
        # --------------------------------------------------------
        # Identify normalization constants (zero mean, unit variance).
        import numpy as np
        from tensorflow import convert_to_tensor
        self.x_shift = convert_to_tensor(np.mean(x, axis=0))
        x -= self.x_shift
        self.x_scale = np.sqrt(np.var(x, axis=0))
        self.x_scale[self.x_scale == 0] = 1.0
        self.x_scale = convert_to_tensor(self.x_scale)
        x /= self.x_scale
        if (type(y) == int):
            y_shape = y
            from tensorflow import zeros, ones
            self.y_shift = zeros(y_shape, dtype=x.dtype)
            self.y_scale = ones(y_shape, dtype=x.dtype)
        else:
            y_shape = y.shape[-1]
            # Normalize the data values.
            self.y_shift = convert_to_tensor(np.mean(y, axis=0))
            y -= self.y_shift
            self.y_scale = np.sqrt(np.var(y, axis=0))
            self.y_scale[self.y_scale == 0] = 1.0
            self.y_scale = convert_to_tensor(self.y_scale)
            y /= self.y_scale
        # --------------------------------------------------------
        # Build the MLP.
        #   Import keras layers for training.
        from tensorflow.keras.models import Sequential
        from tensorflow.keras.layers import Dense, Activation
        from tensorflow.keras.constraints import max_norm
        # Construct the model layers, starting with the input layer.
        model_layers = [Dense(layers[0], input_shape=(x.shape[-1],),
                              activation=activation)]
        # Add all internal layers.
        for size in layers[1:]:
            model_layers += [
                Dense(size, activation=activation,
                      kernel_constraint=max_norm(size)),
            ]
        # Add output layer.
        model_layers += [Dense(y_shape)]
        # Construct and fit the model.
        self.mlp = Sequential(model_layers)
        # --------------------------------------------------------


    # Fit this MLP model to the given data (or loss function).
    # 
    #  x -> numpy ndarray or equivalent
    #  y -> (numpy ndarray or equivalent) OR (integer, number of outputs)
    # 
    def fit(self, x, y, **user_kwargs):
        # Build the internal model if it doesn't exist.
        if (not hasattr(self, "mlp")):
            self.build_model(x,y,**user_kwargs)
        # Set the defaults, then overwrite them with user keyword arguments.
        kwargs = self.default_settings.copy()
        kwargs.update(user_kwargs)
        # Convert the needed keyword arguments into local variables.
        loss = kwargs.pop('loss')
        optimizer = kwargs.pop('optimizer')
        batch_size = kwargs.pop('batch_size')
        epochs = kwargs.pop('epochs')
        show = kwargs.pop('show')
        verbose = kwargs.pop('verbose')
        log_dir = kwargs.pop('log_dir')
        log = kwargs.pop('log')
        # Set the tensorflow warning messages (about sub-optimality).
        #   0 = all messages are logged (default behavior)
        #   1 = INFO messages are not printed
        #   2 = INFO and WARNING messages are not printed
        #   3 = INFO, WARNING, and ERROR messages are not printed
        import os
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = log

        # Retrieve a default optimizer if it exists.
        if (optimizer in self.optimizers):
            optimizers = self.optimizers[optimizer]
        # Go for a larger batch size (for speed), but don't use all
        # data (assume some stochasticity in the problem).
        if type(batch_size) == type(None):
            batch_size = (len(x) // 4) + 1
        # Normalize the incoming data.
        x -= self.x_shift
        x /= self.x_scale
        # Rescale the y, if appropraite.
        if (hasattr(y, "shape")):
            y -= self.y_shift
            y /= self.y_scale

        # Enable logging of progress for visualization..
        tensorboard_callback = {}
        if (show):
            from tensorflow.keras.callbacks import TensorBoard
            import time
            log_dir = os.path.join(log_dir, time.strftime('%Y-%m-%d_%H-%M-%S'))
            tensorboard = TensorBoard(log_dir=log_dir, histogram_freq=epochs)
            tensorboard_callback["callbacks"] = tensorboard
            tensorboard.set_model(self.mlp)

        # Use mean squared error, compile loss function.
        self.mlp.compile(optimizer=optimizer, loss=loss)

        # If y is only a shape
        if ((type(y) == int) or (type(y) == tuple)):
            # Train the network on single batches.
            for i in range(epochs):
                # - need to reduce to batch-size here
                losses = self.mlp.train_on_batch(x, [x])
                print("% 5d  loss %.3e"%(i+1,losses),flush=True)
                # if (not i % 5): print("% 5d  loss %.3e"%(i,losses),flush=True)
                if (show): tensorboard.on_epoch_end(i, named_losses(self.mlp, losses))
            if (show): tensorboard.on_train_end(None)
        else:
            # Fit the model on the new data.
            #   (batches and epochs handled by tensorflow)
            self.mlp.fit(x, y, epochs=epochs, batch_size=batch_size,
                         verbose=verbose, **tensorboard_callback)

        # Open a tensorboard visualizer to show how training went..
        if (show): self.tensorboard(log_dir)


    # Call this model, doing appropriate normalization (and denormalization)
    # of inputs and outputs. Matches type of input, passes that to the
    # MLP, then matches natural output type of MLP.
    def __call__(self, x, *args, **kwargs):
        from tensorflow import cast
        # Normalize the input data and denormalize the output prediction.
        x_shift = cast(self.x_shift, dtype=x.dtype)
        x_scale = cast(self.x_scale, dtype=x.dtype)
        x = (x - x_shift) / x_scale
        y = self.mlp(x, *args, **kwargs)
        y_shift = cast(self.y_shift, dtype=y.dtype)
        y_scale = cast(self.y_scale, dtype=y.dtype)
        return (y * y_scale) + y_shift


    # Run tensorboard on the given log directory.
    def tensorboard(self, log_dir, server_run_seconds=5, server_start_seconds=15):
        from subprocess import Popen, PIPE
        from threading import Timer
        import sys, time, webbrowser
        # Open a process that runs the tensorboard server.
        p = Popen([sys.executable, '-m', 'tensorboard.main',
                   '--logdir', os.path.abspath(log_dir)],
                  stderr=PIPE)
        # Create a timeout in case this fails for some reason..
        t = Timer(server_start_seconds, p.kill)
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
            time.sleep(server_run_seconds)
        # Cancel the timer.
        finally: t.cancel()
        # Kill the server (hopefully user isn't upset!)
        p.kill()

        
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
        scale_path = os.path.join(directory, 'scaling_factors.pkl')
        import pickle
        with open(scale_path, "wb") as f:
            pickle.dump((self.x_shift, self.x_scale, self.y_shift, self.y_scale), f)


    # Given a path to load, load all directory contents.
    def load(self, directory):
        import os
        if (not os.path.exists(directory)):
            raise(IOError("No directry exists.\n  '"+directory+"'"))            
        elif (not os.path.isdir(directory)):
            raise(IOError("Gave path to existing file, need a directry.\n  '"+directory+"'"))
        # Load the model.
        model_path = os.path.join(directory, 'model')
        from keras.models import load_model
        try:    self.mlp = load_model(model_path)
        except: self.mlp = load_model(model_path, compile=False)
        # Load the scaling factors.
        scale_path = os.path.join(directory, 'scaling_factors.pkl')
        import pickle
        with open(scale_path, "rb") as f:
            self.x_shift, self.x_scale, self.y_shift, self.y_scale = pickle.load(f)


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
