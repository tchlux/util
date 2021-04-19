from numpy import zeros, ones, dot, sum, abs, max, argmax, clip, \
    random, prod, asarray, set_printoptions, unravel_index

# Generate a random uniform number (array) in range [0,1].
def zero(*shape):    return zeros(shape)
def randnorm(*shape): return random.normal(size=shape)
def randuni(*shape): return random.random(size=shape)
def randint(*shape, min=-3, max=9):
    data = asarray(random.randint(min+1,max+1,size=shape), dtype=float)
    data[data <= 0] -= 1
    return data

# Build a model given four integers:
#   di - dimension of input
#   ds - dimension for each internal state
#   ns - number of internal states
#   do - dimension of output               
def build_model(di, ds, ns, do):
    # Use random normal vectors.
    input_params    = randnorm(di+1, ds)
    internal_params = randnorm(ds+1, ds, ns-1)
    output_params   = randnorm(ds+1, do)
    # Normalize the length of all random normal vectors for input.
    input_params[:,:] /= ((input_params[:,:]**2).sum(axis=0))**(1/2)
    internal_params[:,:,:] /= ((internal_params[:,:,:]**2).sum(axis=0))**(1/2)
    # Set the bias values.
    input_params[-1,:] = 1
    internal_params[-1,:,:] = 0
    # Set the scratch space for storing internal values to zero.
    internal_values = zero(ds,  ns)
    return input_params, internal_params, output_params, internal_values

# Get the shape of a model (when provided the arrays).
def get_shape(*model):
    di, ds = model[0].shape
    di -= 1
    ns = model[1].shape[-1] + 1
    do = model[2].shape[-1]
    return di, ds, ns, do

# Function for pushing values forward through a dense MLP.
def forward(inputs, input_params, internal_params, output_params, internal_values, display=False):
    di, ds, ns, do = get_shape(input_params, internal_params, output_params)
    # Compute the input layer.
    internal_values[:,0] = clip(dot(inputs, input_params[:di,:]) +
                                input_params[di,:], 0.0, float('inf'))
    if display:
        print("^"*70)
        print("input: ",inputs)
        print()
        for n in range(ds):
            print(f"0.{n}    ", input_params[:di,n], '+', input_params[di,n], '=', internal_values[n,0])
        print(" 0 out ", internal_values[:,0])
    # Compute the next set of internal values with a rectified activation.
    for i in range(ns-1):
        internal_values[:,i+1] = internal_params[ds,:,i] + \
                                 dot(internal_values[:,i], 
                                     internal_params[:ds,:,i])
        if display:
            print()
            for n in range(ds):
                print(f"{i+1}.{n}    ", internal_params[:ds,n,i], '+', internal_params[ds:ds+1,n,i], '=', internal_values[n,i+1])
        internal_values[:,i+1] = clip(internal_values[:,i+1], 0.0, float('inf'))
        if display: print(f" {i+1} out ", internal_values[:,i+1])
    # compute the output.
    output = dot(internal_values[:,ns-1], output_params[:ds]) + output_params[ds]
    if display:
        print()
        for n in range(do):
            print(f"{ns}.{n}    ", output_params[:ds,n],'+', output_params[ds,n], '=', output[n])
        print(f" {ns} out ", output[:])
        print()
        print("output:", output)
        print("_"*70)
    return output

# Compute the gradient with respect to all parameters using finite differences.
def gradient(grad, inputs, *model, display=False):
    # Get the model shape.
    di, ds, ns, do = get_shape(*model)
    # Initialize storage for the gradients.
    input_grad = zeros(model[0].shape)
    internal_grad = zeros(model[1].shape)
    output_grad = ones(model[2].shape)
    # Retrieve the model parameters.
    internal_params = model[1]
    output_params = model[2]
    # Retreive the internal values of the model (after executing forwards).
    internal_values = model[-1]
    # Compute the gradient of the last parameters.
    nonzero = internal_values[:,-1].nonzero()
    output_grad[ds,:] = grad[:]
    for i in range(do):
        output_grad[:ds,i] = internal_values[:,-1] * grad[i]
    internal_values[nonzero,-1] = dot(output_params[:ds,:][nonzero], grad)
    if display:
        print("^"*70)
        print("Output grad:")
        print("",output_grad.T)
        print("",nonzero, internal_values[:,-1])
    # Compute the gradient of all internal parameters.
    for i in range(ns-2,-1,-1):
        # Compute the gradient for all weights.
        #  set the bias gradient.
        internal_grad[ds,:,i] = internal_values[:,i+1]
        #  set the gradient for each column of connections
        #    (to a single output in next layer).
        nonzero = internal_values[:,i].nonzero()
        for j in range(ds):
            if (internal_values[j,i+1] == 0): continue
            internal_grad[:ds,j,i][nonzero] = internal_values[nonzero,i] * internal_values[j,i+1]
            if display:
                print(f"layer {i} -> {i+1}, output node {j}")
                print("  ",internal_grad[:,j,i])
        # Compute the next preceding layer of internal values.
        internal_values[nonzero,i] = dot(internal_params[:ds,:,i][nonzero], internal_values[:,i+1])
        if display:
            print("Grads for next layer:")
            print("",nonzero, internal_values[:,i])
    # Compute the gradient for the input parameters.
    input_grad[di,:] = internal_values[:,0]
    for i in range(ds):
        input_grad[:di,i] = inputs[:] * internal_values[i,0]
    if display:
        print("Input grad:")
        print(input_grad.T)
        print("_"*70)
    # Return the gradients.
    return input_grad, internal_grad, output_grad


# Compute the gradient with respect to all parameters using finite differences.
def finite_difference(inputs, *model, diff=0.0001, display=False):
    # Shift matrices (used for computing finite differences).
    input_shift    = zeros(model[0].shape)
    internal_shift = zeros(model[1].shape)
    output_shift   = zeros(model[2].shape)
    # Function for producting the shifted model.
    shifted_model = lambda: (model[0]+input_shift, model[1]+internal_shift, model[2]+output_shift, model[3])
    # Gradient matrices.
    input_grad    = zeros(model[0].shape)
    internal_grad = zeros(model[1].shape)
    output_grad   = zeros(model[2].shape)
    # Total number of outputs.
    output_shape = forward(inputs, *model).shape
    num_outputs  = prod(output_shape)
    # Compute the expected set of nonzero internal activations.
    forward(inputs, *model)
    expected_nonzero = tuple(model[-1].nonzero()[0])
    # Function for testing the effect that a shift 
    def measure_layer(layer, grad, shift, name):
        for j in range(layer.size):
            curr_idx = unravel_index(j, layer.shape)
            shift[curr_idx] = diff/2
            out_high = forward(inputs, *shifted_model())[out_index]
            nonzero_high = tuple(model[3].nonzero()[0])
            shift[curr_idx] = -diff/2
            out_low  = forward(inputs, *shifted_model())[out_index]
            nonzero_low = tuple(model[3].nonzero()[0])
            shift[curr_idx] = 0
            # If a zero became nonzero (or vice versa), then the
            # finite different approximation is unstable.
            if ((len(nonzero_high) <= len(expected_nonzero)) and
                (len(nonzero_low) <= len(expected_nonzero))):
                # Compute the gradient
                grad[curr_idx] += sum(out_high - out_low) / diff
            if display:
                print(f"{name:14s}{str(curr_idx):10s} {grad[curr_idx]: .3f}")
                print(f"              {float(out_high)}")
                print(f"              {float(out_low)}")
                print(f"              {float(diff)}")
    # Display information.
    if display:
        print("^"*70)
        print("shifted_model: ",[v.shape for v in shifted_model()])
        print("output shape, size: ", output_shape, num_outputs)
    # Cycle over each output.
    for i in range(num_outputs):
        out_index = unravel_index(i, output_shape)
        if display: print("out_index: ",out_index)
        # Cycle over all model parameters, testing effect on output.
        #  input layer
        measure_layer(model[0], input_grad, input_shift, "input idx:")
        #  internal layers
        measure_layer(model[1], internal_grad, internal_shift, "internal idx:")
        #  output layer
        measure_layer(model[2], output_grad, output_shift, "output idx:")
    if display: print("_"*70)
    # Done computing finite difference gradient!
    return input_grad, internal_grad, output_grad


def test():
    print("Testing..")
    di_vals = (1,2,3)
    ds_vals = (1,2,3)
    ns_vals = (1,2,3)
    do_vals = (1,2,3)
    seeds = list(range(5))
    # Cycle all combination of tests.
    from itertools import product
    for (di, ds, ns, do, seed) in product(di_vals, ds_vals, ns_vals, do_vals, seeds):
        # --------------------------------------------------------------------
        # di - dimension of input                
        # ds - dimension for each internal state 
        # ns - number of internal states         
        # do - dimension of output               
        # --------------------------------------------------------------------
        random.seed(seed)

        # Create the model.
        # 
        model = build_model(di, ds, ns, do)
        # Call the "forward" function.
        # inputs = randuni(di)
        inputs = randuni(di)
        # inputs = randint(di)
        # Run the model forward to compute internal values.
        output = forward(inputs, *model, display=False)
        # Compute the gradients with a finite difference.
        approx_model_grad = finite_difference(inputs, *model, display=False)
        # Run the model again (fresh) to get the "internal values".
        output = forward(inputs, *model, display=False)
        # Print the model gradients that were directly computed.
        model_grad = gradient(ones(do), inputs, *model)

        # Check the correctness of the gradient function.
        for i,(app, true) in enumerate(zip(approx_model_grad, model_grad)):
            diff = (abs(app - true) / (abs(true) + 1)).T
            # Skip "internal params" if that is empty.
            if (len(diff) == 0): continue
            # Check for the difference.
            if (max(diff) > .01):
                set_printoptions(precision=3, sign=" ")
                print()
                print("ERROR ON TEST")
                print(" seed =",seed)
                print()
                print("di, ds, ns, do: ",di, ds, ns, do)
                print("input_params:    ",model[0].shape)
                print("internal_params: ",model[1].shape)
                print("output_params:   ",model[2].shape)
                print("internal_values: ",model[3].shape)
                print()
                # forward(inputs, *model, display=True)
                finite_difference(inputs, *model, display=True)
                print()
                print("model[0]:")
                print(model[0].T)
                print()
                print("model[1]:")
                print(model[1].T)
                print()
                print("model[2]:")
                print(model[2].T)
                print()
                print("internals:")
                print(model[-1].T)
                print()
                print()
                print("approx_model_grad[0]:")
                print(approx_model_grad[0].T)
                print()
                print("approx_model_grad[1]:")
                print(approx_model_grad[1].T)
                print()
                print("approx_model_grad[2]:")
                print(approx_model_grad[2].T)
                print()
                print()
                print("model_grad[0]:")
                print(model_grad[0].T)
                print()
                print("model_grad[1]:")
                print(model_grad[1].T)
                print()
                print("model_grad[2]:")
                print(model_grad[2].T)
                print()
                print()
                print("Phase",i,"(0 = input, 1 = internal, 2 = output)")
                print("",max(diff))
                print("",unravel_index(argmax(diff), diff.shape))
                print()
                print("Finite differene gradient:")
                print(app.T)
                print()
                print("Directly computed gradient:")
                print(true.T)
                print()
                print("Difference")
                print(diff)
                print()
                print("ERROR ON TEST")
                exit()

    print(" all passed!")

if __name__ == "__main__": 
    test()






class NN:
    def __init__(self, di, do, ds=16, ns=4):
        self.di = di
        self.ds = ds
        self.ns = ns
        self.do = do
        self.model = list(build_model(di, ds, ns, do))

    def fit(self, x, y, steps=1000, step_factor=0.01, display=False,
            show=False, **kwargs):
        # Make sure that the given data is the right shape.
        assert (self.di == x.shape[-1])
        assert (self.do == y.shape[-1])
        if (show and (self.do == 1) and (self.di == 1)):
            show_interval = max([1, steps // 100])
            from util.plot import Plot
            p = Plot()
            p.add("Data", *(x.T), y.flatten(), group='d', frame=-1)
            p.add_func("Model", self, [x.min(), x.max()], group='m', frame=-1)
            loss_values = []
        # For the number of training steps..
        for s in range(steps):
            if (not s%10): print(s, end="\r")
            if (show): loss_values.append( ((y - self(x))**2).sum()**(1/2) )
            grads = [zeros(l.shape) for l in self.model]
            # Average gradient from all data points.
            for i, (d_in, d_out) in enumerate(zip(x,y)):
                m_out = forward(d_in, *self.model, display=False)
                loss_grad = m_out - d_out
                grad_step = gradient(loss_grad, d_in, *self.model, display=False)
                # Dynamically update the average (of the gradients).
                for j in range(len(grad_step)):
                    grads[j] += (grad_step[j] - grads[j]) / (i+1)
            if display:
                yhat = self(x).reshape(y.shape)
                loss = ((y - yhat)**2).sum(axis=-1).mean()
            # Take a step in the gradient direction.
            for j in range(len(grads)):
                self.model[j] -= grads[j] * step_factor
            # Display progress.
            if display:
                print()
                print("Step:", s)
                print("loss:", loss)
                print("model:")
                for l in self.model[:-1]:
                    print("",l.T)
                print("grads: ")
                for l in grads[:-1]:
                    print("",-l.T)
                print()
            # Update the model plot, if appropriate.
            if (show and (s%show_interval == 0)):
                p.add("Data", *(x.T), y.flatten(), group='d', frame=s)
                p.add_func("Model", self, [x.min(), x.max()], group='m', frame=s)
        # Add the last frame, if it wasn't already added.
        if (show):
            print(" showing plot..")
            # Show the plot of the model.
            p.show(show=False)
            p = Plot("","Step","Loss value")
            p.add("Loss", list(range(len(loss_values))), loss_values, 
                  mode="markers+lines", color=1)
            p.show(append=True, show_legend=False)

    # Return predictions for new data.
    def predict(self, x):
        if (len(x.shape) == 2):
            outputs = []
            for x_in in x:
                outputs.append( forward(x_in, *self.model)[0] )
            return asarray(outputs)
        else: return forward(x, *self.model)
        
    # Wrapper for the "__call__".
    def __call__(self, *args):
        return self.predict(*args)
