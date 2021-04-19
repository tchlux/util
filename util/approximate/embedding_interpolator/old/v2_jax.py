# Function for pushing values forward through a dense MLP.
import jax
import jax.numpy as jnp
from jax.ops import index_update, index
from jax.numpy import zeros, asarray, vstack

def jax_fwd(inputs, input_params, internal_params, output_params, internal_values):
    from v1_numpy import get_shape
    # Get the shape of the arrays.
    di, ds, ns, do = get_shape(input_params, internal_params, output_params)
    # Compute the input layer.
    out = jnp.clip(jnp.dot(inputs, input_params[:di,:]) +
                   input_params[di,:], 0.0, float('inf'))
    # Compute the next set of internal values with a rectified activation.
    for i in range(ns-1):
        out = jnp.clip(internal_params[ds,:,i] + 
                       jnp.dot(out, internal_params[:ds,:,i]),
                       0.0, float('inf'))
    # compute the output.
    output = output_params[ds] + jnp.dot(out, output_params[:ds])
    return output[0]

# Compute the loss as the sum of squared errors.
def jax_loss(x, y, input_params, internal_params, output_params, internal_values):
    total_loss = 0
    for i in range(len(x)):
        total_loss = total_loss + (y[i,0] - jax_fwd(
            x[i], input_params, internal_params, 
            output_params, internal_values))**2
    # Return mean squared error (divided by 2 to make gradient the diff).
    loss = total_loss / (2 * len(x))
    return loss

jax_loss_grad = jax.grad(jax_loss, argnums=[2,3,4])


class NN:
    def __init__(self, di, do, ds=16, ns=4):
        from v1_numpy import build_model
        self.di = di
        self.ds = ds
        self.ns = ns
        self.do = do
        self.model = list(build_model(di, ds, ns, do))

    def fit(self, x, y, steps=1000, step_factor=0.01, display=False, show=False, **kwargs):
        # Make sure that the given data is the right shape.
        assert (self.di == x.shape[-1])
        assert (self.do == y.shape[-1])
        min_max = vstack((x.min(axis=0), x.max(axis=0))).T
        if (show and (self.do == 1) and (self.di == 1)):
            show_interval = max([1, steps // 100])
            from util.plot import Plot
            p = Plot()
            p.add("Data", *(x.T), y.flatten(), group='d', frame=-1)
            p.add_func("Model", self, *min_max, group='m', frame=-1)
            loss_values = []
        # For the number of training steps..
        for s in range(steps):
            if (not s%10): print(s, end="\r")
            if (show): loss_values.append( ((y - self(x))**2).sum()**(1/2) )
            grads = jax_loss_grad(x, y, *self.model)
            # Compute the jax gradients of loss.
            if display:
                print()
                print("Step:", s)
                loss = ((y - self(x))**2).sum()**(1/2)
                print("loss:", loss)
                print("model:")
                for l in self.model:
                    print("",l.T)
                print("grads: ")
                for l in grads[:-1]:
                    print("",-l.T)
                print()
            # Take a step in the gradient direction.
            for j in range(len(grads)):
                self.model[j] -= grads[j] * step_factor
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
                outputs.append( jax_fwd(x_in, *self.model) )
            return asarray(outputs)
        else:
            if (len(x.shape) == 0): x = asarray([x])
            return jax_fwd(x, *self.model)
        
    # Wrapper for the "__call__".
    def __call__(self, *args):
        return self.predict(*args)


