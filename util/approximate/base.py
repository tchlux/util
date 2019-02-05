import numpy as np
from util.math import is_none, is_numeric
from util.algorithms import class_name

# Exceptions for raising useful messages to users.
class ImproperUsage(Exception): pass
class UnexpectedType(Exception): pass
class UnexpectedShape(Exception): pass
class MissingOperator(Exception): pass
class UnderdefinedModel(Exception): pass

# Generic class definition for creating an algorithm that can perform
# regression in the multi-variate setting. In these classes 'x' refers
# to a potentially multi-dimensional set of euclydian coordinates with
# associated single-dimensional response values referred to as 'y'.
# 
# The "fit" function is given both 'x' and associated 'y' and creates
# any necessary (or available) model for predicting 'y' at new 'x'.
# 
# The "predict" function predicts 'y' at potentially novel 'x' values.
# 
# Subclasses are expected to implement _fit and _predict methods which
# take 2D numpy arrays of points as input.
class Approximator:
    _response_dim = 2
    classifier = False

    def __str__(self):
        return class_name(self)

    # Fit a 2D x and a 2D y with this model.
    def fit(self, x, y, classifier=None, *args, **kwargs):
        if (not is_none(classifier)): self.classifier = classifier
        elif (not is_numeric(y[0])):  self.classifier = True
        if ((type(x) != np.ndarray) or (len(x.shape) != 2)):
            raise(UnexpectedType("Provided 'x' should be a 2D numpy array."))
        # If this is a classification problem, convert y values into
        # vertices of a regular simplex.
        if (self.classifier):
            from util.data import regular_simplex
            from util.system import sorted_unique
            self.class_map = sorted_unique(y)
            values = regular_simplex(len(self.class_map))
            y = np.array([values[self.class_map.index(v)] for v in y])
        if (type(y) == list): y = np.array(y)
        if (type(y) != np.ndarray):
            raise(UnexpectedType("Provided 'y' should be a 1D or 2D numpy array."))
        elif (len(y.shape) == 1):
            y = np.reshape(y, (y.shape[0], 1))
            self._response_dim = 1
        elif (len(y.shape) == 2): 
            pass
        else:
            raise(UnexpectedShape("Provided 'y' should be a 1D or 2D numpy array."))
        return self._fit(x.copy(), y.copy(), *args, **kwargs)
        
    # Predict y at new x locations from fit model.
    def predict(self, x, *args, **kwargs):
        if ((type(x) != np.ndarray) or (0 >= len(x.shape) < 2)):
            raise(UnexpectedType("Provided 'x' should be a 1D or 2D numpy array."))
        single_response = len(x.shape) == 1
        if single_response:
            x = np.reshape(x, (1,len(x)))
        response = np.asarray(self._predict(x.copy(), *args, **kwargs), dtype=float)
        if (self.classifier):
            from util.data import category_ratio
            class_response = []
            for guess in response:
                class_response.append(self.class_map[
                    np.argmax(category_ratio(guess)) ])
            response = np.array(class_response)
        # Reduce to one response value if that's what was trained.
        if self._response_dim == 1: response = response[:,0]
        # Reduce to one approximation point if that's what was provided.
        if single_response:         response = response[0]
        # Return the response
        return response

    # Wrapper for 'predict' that returns a single value for a single
    # prediction, or an array of values for an array of predictions
    def __call__(self, *args, **kwargs): return self.predict(*args, **kwargs)



# Generic class definition for an approximator that uses an algorithm
# which predicts new values as a linear combination of provided values.
# These types of approximators are convenient for their ability to
# predict any response value that has defined addition and
# multiplication operators.
# 
# Subclasses should implement "_fit" operator which takes 2D numpy
# arrays of row-points, and the "_predict" operator which takes 2D
# numpy arrays of row-points and returns the indices of source points
# and associated weights required to make a prediction.
class WeightedApproximator(Approximator):
    WeightedApproximator = True
    y = None

    # Fit a 2D x and a 2D y with this model.
    def fit(self, x, y=None, classifier=None, **kwargs):
        if (not is_none(classifier)):                     self.classifier = classifier
        elif (not is_none(y)) and (not is_numeric(y[0])): self.classifier = True
        if ((type(x) != np.ndarray) or (len(x.shape) != 2)):
            raise(UnexpectedType("Provided 'x' should be a 2D numpy array."))
        if (not is_none(y)):
            if (not hasattr(y, "__len__")):
                raise(MissingOperator("Provided 'y' to fit must have defined '__len__' operator."))
            elif (not hasattr(y, "__getitem__")):
                raise(MissingOperator("Provided 'y' to fit must have defined '__getitem__' operator."))
            elif (not hasattr(y[0], "__add__")):
                raise(MissingOperator("Elements of provided 'y' must have defined '__add__' operator."))
            elif (not hasattr(y[0], "__mul__")):
                raise(MissingOperator("Elements of provided 'y' must have defined '__mul__' operator."))
            self.y = y
        # Fit the provided x values.
        return self._fit(x.copy(), **kwargs)
        
    # Predict y at new x locations from fit model.
    def predict(self, x, *args, **kwargs):
        if ((type(x) != np.ndarray) or (0 >= len(x.shape) < 2)):
            raise(UnexpectedType("Provided 'x' should be a 1D or 2D numpy array."))
        single_response = len(x.shape) == 1
        if single_response:
            x = np.reshape(x, (1,len(x)))
        indices, weights = self._predict(x.copy(), *args, **kwargs)
        # Return the indices and weights if no y values were provided.
        if (is_none(self.y)): 
            response = [(ids,wts) for (ids,wts) in zip(indices, weights)]
        else:
            # Collect response values via weighted sums of self.y values
            response = []        
            for ids, wts in zip(indices, weights):
                if self.classifier:
                    val_weights = {}
                    # Sum the weights associated with each category.
                    for i, w in zip(ids, wts):
                        val_weights[self.y[i]] = val_weights.get(self.y[i],0.)+w
                    # Return the category with the largest sum of weight.
                    response.append( max(val_weights.items(), key=lambda i: i[-1])[0] )
                else:
                    # Return the weighted sum of predictions.
                    response.append( sum(self.y[i]*w for (i,w) in zip(ids, wts)) )
        # Reduce to one approximation point if that's what was provided.
        if single_response: response = response[0]
        # Return the response
        return response

    # Wrapper for 'predict' that returns a single value for a single
    # prediction, or an array of values for an array of predictions
    def __call__(self, *args, **kwargs): return self.predict(*args, **kwargs)    

