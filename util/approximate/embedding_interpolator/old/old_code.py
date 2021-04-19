
# From Fortran source code, when checking how model was initialized.
# 
# ! Before the conditional SGEMM (second to last) in model initialization:
# ! CALL PRINT_DISTRIBUTION(SIZE(INTERNAL_VALUES(:,:)), INTERNAL_VALUES(:,:))
# 
# ! SUBROUTINE PRINT_DISTRIBUTION(N, VALUES)
# !   USE ISO_FORTRAN_ENV, ONLY: REAL32, REAL64, INT64
# !   USE FAST_SORT, ONLY: ARGSORT_R64
# !   IMPLICIT NONE
# !   INTEGER, INTENT(IN) :: N
# !   REAL(KIND=REAL32), INTENT(IN), DIMENSION(1:N) :: VALUES
# !   REAL(KIND=REAL64), DIMENSION(1:N) :: SORTED
# !   INTEGER(KIND=INT64), DIMENSION(1:N) :: INDICES
# !   REAL(KIND=REAL64) :: TEMP
# !   INTEGER :: I, J
# !   ! Copy the values into a temporary array.
# !   SORTED(:) = REAL(VALUES(:), REAL64)
# !   CALL ARGSORT_R64(SORTED, INDICES)
# !   ! Compute the mean.
# !   TEMP = SUM(SORTED(:)) / REAL(N, REAL32)
# !   WRITE (*,"('(1:',I6,')')") N
# !   WRITE (*,"('  min:   ', E10.3)") SORTED(1)
# !   WRITE (*,"('  mean:  ', E10.3)") TEMP
# !   WRITE (*,"('  median:', E10.3)") SORTED(N/2)
# !   WRITE (*,"('  max:   ', E10.3)") SORTED(N)
# !   WRITE (*,*) ''
# ! END SUBROUTINE PRINT_DISTRIBUTION



# 2021-02-06 23:48:16
# 
###############################################
# def rand(*shape):                           #
#     global rng_key                          #
#     key, rng_key = random.split(rng_key)    #
#     return random.uniform(key, shape=shape) #
# 
# backward = grad(forward)
# output_gradient = backward(input_params, internal_params, output_params)
# 
# rng_key = random.PRNGKey(0)
# 
###############################################


# 2021-02-07 10:47:37
# 
############################################
# # Disable JAX logging messages.          #
# import os                                #
# os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2" #
############################################


# 2021-02-07 23:42:02
# 
###########################
# print()                 #
# print("model_grad[0]:") #
# print(model_grad[0].T)  #
# print()                 #
# print("model_grad[1]:") #
# print(model_grad[1].T)  #
# print()                 #
# print("model_grad[2]:") #
# print(model_grad[2].T)  #
# print()                 #
###########################


# 2021-02-07 23:42:24
# 
###########################################################
# # Print the model gradients with respect to parameters. #
# print()                                                 #
# print("model_grad[0]:")                                 #
# print(model_grad[0].T)                                  #
# print()                                                 #
# print("model_grad[1]:")                                 #
# print(model_grad[1].T)                                  #
# print()                                                 #
# print("model_grad[2]:")                                 #
# print(model_grad[2].T)                                  #
###########################################################


# 2021-02-07 23:42:44
# 
#################################
# # Print the model parameters. #
# print()                       #
# print("model[0]:")            #
# print(model[0].T)             #
# print()                       #
# print("model[1]:")            #
# print(model[1].T)             #
# print()                       #
# print("model[2]:")            #
# print(model[2].T)             #
#                               #
# # Print the internal values.  #
# print()                       #
# print("internals:")           #
# print(model[-1].T)            #
#################################


# 2021-02-10 16:45:40
# 
                ###############################
                # print("true in:   ",d_in)   #
                # print("true out:  ",d_out)  #
                # print("model out: ",m_out)  #
                # print("grad:      ",grad)   #
                # print("shape:", grad.shape) #
                ###############################


# 2021-02-10 16:45:44
# 
                #####################
                # print()           #
                # print("step: ")   #
                # for l in step:    #
                #     print()       #
                #     print("",l.T) #
                #####################


# 2021-02-10 16:45:48
# 
                ####################################
                # if display: print("*"*70,"\n",i) #
                # print('-'*70)                    #
                # print("i: ",i)                   #
                ####################################


# 2021-02-10 16:50:18
# 
            ########################################
            # if display: print("\n"*3,"Step:", s) #
            ########################################


# 2021-02-13 11:12:50
# 
#######################################
# from util.plot import Plot          #
# import numpy as np                  #
# from scipy.special import erfinv    #
#                                     #
# p = Plot()                          #
#                                     #
# eps = 0.00001                       #
# n = 1000                            #
# x = np.linspace(-1+eps, 1-eps, n)   #
# y1 = np.arctanh(x) / 1.2            #
# y2 = erfinv(x)                      #
#                                     #
# p.add("atanh", x, y1, mode="lines") #
# p.add("ierf", x, y2, mode="lines")  #
# p.show()                            #
#                                     #
# exit()                              #
#######################################


# 2021-02-14 10:24:12
# 
##############################################################################
# from util.math import Spline, Polynomial, fit_polynomial, Fraction         #
#                                                                            #
# # Make all values in a (nested) list exact.                                #
# def exact(l):                                                              #
#     try:    return [exact(v) for v in l]                                   #
#     except: return Fraction(l)                                             #
#                                                                            #
# # p1 = Polynomial(fit_polynomial([0,1], [0,1], dx0=0))                     #
#                                                                            #
# p1 = Polynomial(exact([ 1, 0, 0]))                                         #
# p2 = Polynomial(exact([-2, 6,-3]))                                         #
# p3 = Polynomial(exact([ 1,-6, 9]))                                         #
#                                                                            #
# xs = exact([(0,.5,1), (1, 1.5, 2), (2, 2.5, 3)])                           #
# ys = exact([map(p1, xs[0]), map(p2, xs[1]), map(p3, xs[2])])               #
# # Translate x and y values to desired scale.                               #
# xs = [[v / 3 for v in xi] for xi in xs]                                    #
# ys = [[(v*3)/2 for v in yi] for yi in ys]                                  #
# # Create new polyhnomials.                                                 #
# p1, p2, p3 = [Polynomial(fit_polynomial(xi,yi)) for (xi,yi) in zip(xs,ys)] #
#                                                                            #
# print("xs: ",xs)                                                           #
# print("ys: ",ys)                                                           #
# # exit()                                                                   #
#                                                                            #
# # f = Spline(exact([0,1,2,3]), functions=[p1,p2,p3])                       #
# f = Spline([Fraction(0),Fraction(1,3),Fraction(2,3),Fraction(1)],          #
#            functions=[p1, p2, p3])                                         #
# g = f.integral()                                                           #
# print("f: ",f)                                                             #
# print("g: ",g)                                                             #
#                                                                            #
# p.add_func("f", f, [0,1], color=(0,0,0), line_width=3, dash="dash")        #
# p.add_func("p1", p1, [0,1/3])                                              #
# p.add_func("p2", p2, [1/3,2/3])                                            #
# p.add_func("p3", p3, [2/3,1])                                              #
# p.add_func("g", g, [0,1])                                                  #
# p.show()                                                                   #
# exit()                                                                     #
##############################################################################


#2021-02-14 21:26:40
#
##########################################################################################
# import numpy as np                                                                     #
#                                                                                        #
# d = 1000                                                                               #
# n = 50000                                                                              #
# vecs = np.zeros((n,d))                                                                 #
#                                                                                        #
# from util.system import Timer                                                          #
# # out = nn.plrm.random_unit_vectors(vecs.T)                                            #
# t = Timer()                                                                            #
# t.start()                                                                              #
# out = nn.plrm.random_unit_vectors(vecs.T)                                              #
# print(t.stop())                                                                        #
# exit()                                                                                 #
# # np.random.seed(0)                                                                    #
# # vecs = 2 * np.random.random(size=(n,d)) - 1                                          #
# # vecs = np.arctanh(vecs)                                                              #
# # vecs = (vecs.T / np.linalg.norm(vecs, axis=1)).T                                     #
#                                                                                        #
#                                                                                        #
# # print("vecs: ",vecs)                                                                 #
# # print("out: ",out.T)                                                                 #
# from util.plot import Plot                                                             #
# p = Plot()                                                                             #
#                                                                                        #
# p.add("", *(vecs.T), color=4, marker_line_width=2, marker_line_color=0, marker_size=3) #
# p.show(show_legend=False)                                                              #
# exit()                                                                                 #
##########################################################################################


#2021-02-14 21:26:44
#
##################################################################
# # # Order 2 approximation                                      #
# # shifted_box = 1.5_REAL64 + 1.5_REAL64 * shifted_box          #
# # WHERE ((1.0_REAL64 .LE. shifted_box) .AND. &                 #
# #      (shifted_box .LT. 2.0_REAL64)) shifted_box = (&         #
# #      -(2*shifted_box**2 - 6*shifted_box + 3) )               #
# # WHERE ((2.0_REAL64 .LE. shifted_box) .AND. &                 #
# #      (shifted_box .LE. 3.0_REAL64)) shifted_box = (&         #
# #      (shifted_box - 3)**2)                                   #
# # WHERE (shifted_box .GT. 3.0_REAL64) shifted_box = 0.0_REAL64 #
# # box_vals(pt_num, box_num) = PRODUCT(shifted_box)             #
##################################################################


# 2021-02-15 08:42:50
# 
#######################################################################
# # Add an expected (good) approximation.                             #
# from util.approximate import NeuralNetwork                          #
# print("Getting good model..")                                       #
# good_model = NeuralNetwork()                                        #
# # try:                                                              #
# #     good_model.load("good_test_model")                            #
# #     print(" loaded existing!")                                    #
# # except:                                                           #
# print(" training new!")                                             #
# good_model.fit(x, y, layers=(self.ds,)*self.ns, optimizer='sgd',    #
#                batch_size=len(x), epochs=steps)                     #
# # print(" saving..")                                                #
# # good_model.save("good_test_model")                                #
# s += 1                                                              #
# p.add("Data", *(x.T), y.flatten(), group='d', frame=s)              #
# # p.add_func("Model", self, [x.min(), x.max()], group='m', frame=s) #
# print(" adding to plot..")                                          #
# p.add_func("Good model", good_model, [x.min(), x.max()],            #
#            vectorized=True, group='g', frame=s)                     #
#######################################################################

#2021-02-15 11:53:42
#
########################################################################################################
# from v2_jax import NN as jax_nn                                                                      #
# m_jax = jax_nn(di,do,ds,ns)                                                                          #
# m_jax.model   = [input_params.copy(), internal_params.copy(), output_params.copy(), m_jax.model[-1]] #
# print()                                                                                              #
# m_jax.fit(x, y, steps=2, display=True, show=True)                                                    #
########################################################################################################

