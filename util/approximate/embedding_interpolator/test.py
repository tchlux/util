
# import fmodpy

# f_compiler_args = ['-fPIC', '-shared', '-O2','-lblas', '-fopenmp']
# f_compiler_args = ['-fPIC', '-shared', '-O2',
#                    '-L', '/usr/local/Cellar/openblas/0.3.13/lib',
#                    '-lopenblas', '-fopenmp']
# f_compiler_args = ['-fPIC', '-shared', '-O2', '-fopenmp',
#                    '/Users/thomaslux/Downloads/ATLAS-10.3/Macman/lib/libatlas.a',
#                    '/Users/thomaslux/Downloads/ATLAS-10.3/Macman/lib/libf77blas.a']

# plrm = fmodpy.fimport("plrm.f90", verbose=True,
#                       f_compiler_args=f_compiler_args).plrm

import os
if os.path.exists("plrm/plrm.so"):
    os.remove("plrm/plrm.so")
from plrm import plrm

# import fmodpy
# plrm = fmodpy.fimport("plrm.f90", blas=True, omp=True, verbose=False).plrm


if __name__ == "__main__":
    from numpy import sin, linspace, pi, random, where, append, \
        asarray, zeros, linalg, vstack
    from util.plot import Plot

    random.seed(0)
    N = 2**14
    # N = 2**15
    # x = linspace(.001, 2*pi-.001, N)[:,None]
    D = 2**10
    # D = 2**9
    # D = 2**12
    x = 2 * pi * random.random(size=(N, D))
    y = sin(linalg.norm(x,axis=-1)).reshape(-1,1)

    # # Random normal data.
    # x = random.normal(size=(1000,3))
    # y = x[:,0:1]

    # 2097152

    x = asarray(x, dtype="float32")
    y = asarray(y, dtype="float32")

    # # Normalize the data.
    x -= x.mean(axis=0)[0]
    x /= where(x.std(axis=0)[0] == 0, 1.0, x.std(axis=0)[0])
    y -= y.mean(axis=0)[0]
    y /= where(y.std(axis=0)[0] == 0, 1.0, y.std(axis=0)[0])

    # Create a new model.
    di = x.shape[-1]
    do = y.shape[-1]
    ds = 8
    ns = 4
    steps = 500
    steps_per_plot = 1000
    # ds = 3
    # ns = 2

    print("creating model..")
    # plrm.new_model(di, ds, ns, do)
    # plrm.init_model(inputs=x.T, outputs=y.T, seed=0)
    print("(x.T).shape: ",(x.T).shape)
    print("(y.T).shape: ",(y.T).shape)


    COMPARE_WITH_TF = True
    TEST_PLOT_MODEL = True and (x.shape[1] == 1)
    TEST_MODEL_CONSISTENCY = False
    frames = min(steps, 500)


    if COMPARE_WITH_TF:
        print("Importing tensorflow..", flush=True)
        # Get a tensorflow neural network.
        from util.approximate import NeuralNetwork as nn
        # Get a pure gradient descent optimizer.
        # from tensorflow.compat.v1.train import GradientDescentOptimizer as GD
        # Optimizer = GD(learning_rate=0.01)
        Optimizer = "adam"
        # Get a timer.
        from util.system import Timer
        t = Timer()
        steps_1 = steps
        # Measure training time..
        print()
        print("Training models..", flush=True)
        print()
        # --------------------
        t.start()
        plrm.new_model(di, ds, ns, do)
        plrm.init_model(inputs=x.T, outputs=y.T, seed=0)
        plrm.minimize_mse(x.T, y.T, steps=steps_1)
        t.stop()
        f_time_1 = t.total
        print("Fortran:   ", f_time_1)
        # --------------------
        # t.start()
        model = nn()
        model.fit(x, y, optimizer=Optimizer, batch_size=len(x),
                  layers=(ds,)*ns,
                  show=False, epochs=steps_1)
        t_time_1 = model.fit_time
        # t.stop()
        # t_time_1 = t.total
        print("Tensorflow:", t_time_1)
        # --------------------
        # Measure execution time..
        out = zeros(y.shape, dtype="float32")
        # plrm.evaluate_one(x[0], out)
        t.start()
        plrm.evaluate(x.T, out.T)
        t.stop()
        print("F eval:", t.total)        
        t.start()
        model.mlp(x)
        t.stop()
        print("T eval:", t.total)
        print()
        # Plot the resulting models if possible.
        if (x.shape[-1] in {1,2}):
            from util.plot import Plot
            p = Plot()
            p.add("Data", *(x.T), y.flatten())

            def eval_fort(x):
                x = asarray(x.reshape((di,)), dtype="float32")
                output = zeros(y.shape[-1], dtype="float32")
                plrm.evaluate_one(x, output)
                return output[0]
            def eval_tens(x):
                x = asarray(x.reshape((di,)), dtype="float32")
                return model(x)
            min_max = vstack((x.min(axis=0), x.max(axis=0))).T
            p.add_func("Fort", eval_fort, *min_max)
            p.add_func("Tens", eval_tens, *min_max)
            p.show()


    if TEST_PLOT_MODEL:
        # Initialize the model if it has not been initialized yet.
        if (plrm.input_params is None):
            plrm.new_model(di, ds, ns, do)
            plrm.init_model(inputs=x.T, outputs=y.T, seed=0)

        # Create a function for evaluating the model.
        def f(xi):
            xi = asarray(xi, dtype="float32").reshape((di,))
            yi = zeros((do,), dtype="float32")
            plrm.evaluate_one(xi,yi)
            return yi[0]

        # Create a plot.
        from util.plot import Plot
        p = Plot()
        print("x: ",x.T[0])
        print("y: ",y.T[0])
        print("fitting model..")

        show_interval = max([1, steps // frames])
        mse = []
        for i in range(0, steps+1, steps_per_plot):
            if (not i%show_interval) or (i == steps-1):
                print(f"i: {i} {(mse+[0])[-1]:.4f}",end="\r")
                p.add("data", x.T[0], y.T[0], group="data", frame=i)
                p.add_func("f", f, [x.min(), x.max()], plot_points=100, group='f', frame=i)
            result = plrm.minimize_mse(x.T, y.T, steps=steps_per_plot, keep_best=True)
            mse.append( result )
        print()
        p.show(show=False)
        # Add loss function plot.
        p = Plot()
        p.add("MSE", list(range(len(mse))), mse, mode="lines", color=1)
        p.show(append=True)



    if TEST_MODEL_CONSISTENCY:
        def get_plrm_model():
            # Get the input parameters (in numpy / jax model format)
            params = plrm.get_input_params()
            bias = plrm.get_input_bias()
            input_params = append(params, bias.reshape(1,*bias.shape), axis=0)
            # Get the internal parameters (in numpy / jax model format)
            params = plrm.get_internal_params()
            bias = plrm.get_internal_bias()
            internal_params = append(params, bias.reshape(1,*bias.shape), axis=0)
            # Get the output parameters (in numpy / jax model format)
            params = plrm.get_output_params()
            bias = plrm.get_output_bias()
            output_params = append(params, bias.reshape(1,*bias.shape), axis=0)
            # Return the paremeters.
            return input_params, internal_params, output_params, zeros((ds,ns), dtype="float32")

        def get_plrm_grads(x, y):
            # Initialize
            sse = 0.0
            ip = zeros((di,ds),      dtype="float32", order='F')
            ib = zeros((ds,),        dtype="float32", order='F')
            mp = zeros((ds,ds,ns-1), dtype="float32", order='F')
            mb = zeros((ds,ns-1),    dtype="float32", order='F')
            op = zeros((ds,do),      dtype="float32", order='F')
            ob = zeros((do,),        dtype="float32", order='F')
            # Compute the gradient of all parameters with respect to SSE.
            sse = plrm.sse_gradient(x, y, sse, ip, ib, mp, mb, op, ob)[0]
            # Concatenate parameters to same shape as numpy.
            input_params = append(ip, ib.reshape(1,*ib.shape), axis=0)
            internal_params = append(mp, mb.reshape(1,*mb.shape), axis=0)
            output_params = append(op, ob.reshape(1,*ob.shape), axis=0)
            # Return the three matrices of weights.
            return input_params, internal_params, output_params

        # Create known models with the same parameters as the fortran code.
        from v1_numpy import NN as numpy_nn
        m_numpy = numpy_nn(di,do,ds,ns)
        m_numpy.model = [p.copy() for p in get_plrm_model()]

        steps = 1

        m_numpy.fit(x, y, steps=steps, display=True, show=False)

        for s in range(steps):
            plrm_grads = get_plrm_grads(x.T, y.T)
            loss = plrm.minimize_mse(x.T, y.T, steps=1)
            plrm_model = get_plrm_model()
            print()
            print("Step:", s)
            print("loss:", loss)
            print("model:", [l.shape for l in plrm_model[:-1]])
            for l in plrm_model[:-1]:
                print("",l.T)
            print("grads: ", [l.shape for l in plrm_grads])
            for l in plrm_grads:
                print("",-l.T)
            print()


 







#2021-02-15 11:55:08
#
    ############################################
    # # plrm.fit_data(x.T, y.T, 10000)         #
    # # print("evaluating model..")            #
    # # p.add("data", x.T[0], y.T[0])          #
    # # p.add_func("f", f, [x.min(), x.max()]) #
    ############################################


#2021-02-15 11:56:09
#

    #########################################################################
    # # Initialize the model (with the shapes).                             #
    # m = NN(x.shape[-1], y.shape[-1], ds=8, ns=4)                          #
    #                                                                       #
    # print()                                                               #
    # print("model: ")                                                      #
    # for l in m.model:                                                     #
    #     if (l.size == 0): continue                                        #
    #     print()                                                           #
    #     print(l.T)                                                        #
    # print()                                                               #
    #                                                                       #
    # # Show the output of the initial forward computation.                 #
    # out = forward(x[len(x)//2], *m.model, display=False)                  #
    # g = gradient(ones(*out.shape), x[len(x)//2], *m.model, display=False) #
    #                                                                       #
    # # out2 = jax_fwd(x[len(x)//2], *m.model)                              #
    # # g2 = jax_grad(x[len(x)//2], *m.model)                               #
    # # print("out2: ",out2)                                                #
    # # print("g2: ")                                                       #
    # # print(g2[0])                                                        #
    # # print(g2[1])                                                        #
    # # print(g2[2])                                                        #
    # # exit()                                                              #
    #                                                                       #
    # print()                                                               #
    # # Test the model before training.                                     #
    # my = m(x).reshape(y.shape)                                            #
    # error = y - my                                                        #
    # print()                                                               #
    # print("my:    ",my.T[0])                                              #
    # print("error: ",error.T[0])                                           #
    # print(sum(abs(error)) / len(error))                                   #
    # print()                                                               #
    #                                                                       #
    # # Train the model, making plots.                                      #
    # m.fit(x,y, steps=1000, display=False, show=True)                      #
    #                                                                       #
    # # Test the model after training.                                      #
    # my = m(x).reshape(y.shape)                                            #
    # error = y - my                                                        #
    # print()                                                               #
    # print("my:    ",my.T[0])                                              #
    # print("error: ",error.T[0])                                           #
    # print(sum(abs(error)) / len(error))                                   #
    # print()                                                               #
    #                                                                       #
    # # Show the output of the initial forward computation.                 #
    # out = forward(x[len(x)//2], *m.model, display=False)                  #
    # gradient(ones(*out.shape), x[len(x)//2], *m.model, display=False)     #
    # print()                                                               #
    #                                                                       #
    # print()                                                               #
    # print("model: ")                                                      #
    # for l in m.model:                                                     #
    #     if (l.size == 0): continue                                        #
    #     print()                                                           #
    #     print(l.T)                                                        #
    # print()                                                               #
    #########################################################################


#2021-02-15 19:28:04
#
    ##################################################
    # from util.approximate.testing import test_plot #
    # p, x, y = test_plot(model)                     #
    # p.show()                                       #
    ##################################################


# 2021-02-20 15:11:53
# 
    ################################################################################
    # TEST_TRANSFORMATIONS = False                                                 #
    # if TEST_TRANSFORMATIONS:                                                     #
    #     if (1 <= ds <= 3):                                                       #
    #         p = Plot()                                                           #
    #         p.add("Data", *(x.T), marker_line_width=2, marker_size=4)            #
    #     # Holder for output values for "evaluate".                               #
    #     _ = y[0].copy()                                                          #
    #     # Holder for transformation.                                             #
    #     ls = []                                                                  #
    #     for i in range(ns):                                                      #
    #         ls.append( zeros((len(x),ds)) )                                      #
    #     # Get all the transformed "x" values.                                    #
    #     for i in range(len(x)):                                                  #
    #         plrm.evaluate(x[i], _)                                               #
    #         internals = plrm.get_internal_values()                               #
    #         for j in range(ns):                                                  #
    #             ls[j][i,:] = internals[:,j]                                      #
    #     # Print out stats on each layer.                                         #
    #     for j in range(ns):                                                      #
    #         print()                                                              #
    #         print(f"l{j+1}:")                                                    #
    #         print("", "mean", ls[j].mean(axis=0))                                #
    #         print("", "var ", ls[j].var(axis=0))                                 #
    #         if (1 <= ds <= 3):                                                   #
    #             p.add(f"l{j+1}", *(ls[j]).T, marker_line_width=2, marker_size=4) #
    #     if (1 <= ds <= 3):                                                       #
    #         p.show()                                                             #
    ################################################################################


# 2021-02-24 09:08:28
# 
        ##################
        # # i = 8        #
        # # x = x[i:i+1] #
        # # y = y[i:i+1] #
        ##################



# 2021-03-06 15:17:38
# 
################################################################################
# additional_steps = 100
# # --------------------                                                       #
# # Measure training time for more points..                                    #
# # --------------------                                                       #
# steps_2 = steps_1 + additional_steps                                         #
# t.start()                                                                    #
# plrm.new_model(di, ds, ns, do)                                               #
# plrm.init_model()                                                            #
# plrm.minimize_mse(x.T, y.T, steps_2)                                         #
# t.stop()                                                                     #
# f_time_2 = t.total                                                           #
# print("Fortran per iteration:", (f_time_2 - f_time_1) / (steps_2 - steps_1)) #
# # --------------------                                                       #
# model = nn()                                                                 #
# model.fit(x, y, optimizer="sgd", batch_size=len(x),                          #
#           layers=(ds,)*ns,                                                   #
#           show=False, epochs=steps_2)                                        #
# t_time_2 = model.fit_time                                                    #
# print("TF per iteration:     ", (t_time_2 - t_time_1) / (steps_2 - steps_1)) #
################################################################################


# 2021-04-19 08:27:14
# 
    ###############################################################
    # # plrm.new_model(di, ds, ns, do, rectifier=-5.0, bias=-5.0) #
    # # plrm.get_input_params()[:,:] = 1.0                        #
    # # plrm.get_internal_params()[:,:,:] = 1.0                   #
    # # plrm.get_output_params()[:,:] = 1.0                       #
    # # plrm.get_input_bias()[:] = 1.0                            #
    # # plrm.get_internal_bias()[:,:] = 1.0                       #
    # # plrm.get_output_bias()[:] = 1.0                           #
    ###############################################################
