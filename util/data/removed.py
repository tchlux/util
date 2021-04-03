# TODO:  Remove builtin prediction code for anything other than
#        nearest neighbor prediction. Maybe remove prediction in
#        general, since it isn't something that can be picked beforehand.
# TODO:  Running 'predict' with nearest neighbor seg-faults after 4th
#        round, might involve conditioning, probably involves ball tree.
# TODO:  'predict' does not work with a view, picks wrong column, it
#        also appears to pick the wrong column even when not a view.
# TODO:  Rewrite "predict" to work given a data object. This will
#        break up the data into collections of unique missing values.
#        Then for each collection, build a model over the appropriate
#        matrix-format data and predict that collection. This will
#        work more nicely with the k-fold method.
# TODO:  Make "predict" not require that target columns be convertible
#        to numeric form, only that they operate under a weighted sum.
# TODO:  Remove the "fill" method in favor of "predict".
# 
# 


# DEPRECATION WARNING -- THIS IS ONLY FOR BACKWARDS COMPATIBILITY
def deflate(self, name):
    import warnings
    warnings.warn("\n  The 'deflate' function is being deprecated, use 'pack' instead.")
    return self.pack(name)

# DEPRECATION WARNING -- THIS IS ONLY FOR BACKWARDS COMPATIBILITY
def inflate(self, column):
    import warnings
    warnings.warn("\n  The 'inflate' function is being deprecated, use 'unpack' instead.")
    return self.unnpack(column)


# Use a prediction model to fill missing values in a provided row.
def fill(self, row, model=None, weights=None, bias={}):
    import numpy as np
    numeric = self.to_matrix()
    # Pick the filling model based on data size.
    if (model is None):
        # Use Voronoi if there are fewer than 10000 points.
        if (len(self) <= 5000):
            from util.approximate import Voronoi, unique
            model = unique(Voronoi)()
        # Otherwise use nearest neighbor (can't get points + weights from NN)
        else:
            from util.approximate import unique, condition
            samples = min(numeric.shape[0], 10000)
            dim = min(numeric.shape[1], 1000)
            # Use nearest neighbor (can't get points + weights from NN)
            from util.approximate import NearestNeighbor
            model = condition(unique(NearestNeighbor), dim=dim, samples=samples)()
    # Set the weights accordingly
    if (weights is None):
        weights = np.ones(numeric.shape[1])
    elif (len(weights) != self.shape[1]):
        raise(Data.BadValue("Weights vector is length {len(weights)}, expected length {self.shape[1]}."))
    else:
        # Convert the provided weight to the correct shape.
        col_weights, weights = weights, []
        col_inds = list(range(len(numeric.names)))
        while (len(col_inds) > 0):
            if (col_inds[0] not in numeric.cats):
                col_inds.pop(0)
                weights.append( col_weights.pop(0) )
            else:
                for i in range(1,len(col_inds)):
                    not_categorical = (col_inds[i] not in numeric.cats)
                    beginning_of_next = ("1" == 
                        numeric.names[col_inds[i]].split('-')[-1])
                    if (not_categorical or beginning_of_next): break
                else: i += 1
                col_inds = col_inds[i:]
                weights += [col_weights.pop(0)] * i
        weights = np.array(weights)
    # Get the indices of the known and unknown values in this row.
    row_known = [i for i in range(len(row)) if (row[i] is not None)]
    row_unknown = [i for i in range(len(row)) if (row[i] is None)]
    if len(row_unknown) == 0: raise(Data.BadValue("Provided row had no missing values."))
    # Get the numeric version of the provided row.
    numeric_row = numeric.to_real(row, fill=True)
    # Identify the indices of known values and unknown values.
    known_values = np.array([i for i in range(len(numeric_row))
                             if (numeric_row[i] is not None)])
    unknown_values = np.array([i for i in range(len(numeric_row))
                               if (numeric_row[i] is None)])
    # Get the part of the numeric row
    no_none_numeric_row = np.array([numeric_row[i] for i in known_values])
    # Normalize the numeric row (to debias numeric scale differences)
    no_none_numeric_row = ((no_none_numeric_row - numeric.shift[known_values])
                           / numeric.scale[known_values]) * weights[known_values]
    # Get the normalized training point locations.
    normalized_numeric_data = ((numeric - numeric.shift)
                               / numeric.scale) * weights
    # Fit the model and make a prediction for this row.
    if hasattr(model, "WeightedApproximator"):
        # Fit the model (just to the points)
        that_type = type(normalized_numeric_data[:,known_values])
        model.fit(normalized_numeric_data[:,known_values])
        del(normalized_numeric_data)
        # Get the source points and source weights for the prediction
        pts, wts = model.predict(np.array(no_none_numeric_row))
        # Compute the fill row by performing the weighted sum.
        fill_row = np.sum((numeric[pts,:][:,unknown_values].T * wts), axis=1)
        # Convert points and weights into python list (in case they're not)
        pts = list(pts)
        wts = list(wts)
    else:
        # The model makes predictions using all points.
        # Fit the model to the input points and output points.
        model.fit(normalized_numeric_data[:,known_values], 
                  normalized_numeric_data[:,unknown_values])
        del(normalized_numeric_data)
        pts = list(range(numeric.data.shape[0]))
        wts = [1 / len(pts)] * len(pts)
        # Use the model to make a prediction.
        fill_row = model.predict(no_none_numeric_row)
        # De-normalize the predicted output.
        fill_row = ((fill_row / weights[unknown_values]) * 
                    numeric.scale[unknown_values] - 
                    numeric.shift[unknown_values])
    # Pad fill row with 0's so that it is the correct length for conversion.
    fill_row = [0.] * len(known_values) + list(fill_row)
    # Transfer the filled numeric row into the original row
    fill_row = numeric.from_real(fill_row, bias=bias)
    for i,val in enumerate(fill_row):
        if (row[i] is None):
            # Set the value in the row to be the predicted fill value.
            row[i] = val
        else:
            # Backwards correct predicted output to have known values.
            fill_row[i] = row[i]
    prediction_model = str(model)
    # Construct a "Prediction" object with info about the prediction.
    class Prediction:
        parent = self
        model = prediction_model
        given = row_known
        filled = row_unknown
        predictors = pts
        weights = wts
        data = fill_row
        def __str__(self): 
            string =  f"Prediction made with {self.model}\n"
            string += f"  Given {[self.parent.names[i] for i in self.given]}\n"
            string += f"   {[self.data[i] for i in self.given]}\n"
            string += f"  Filled {[self.parent.names[i] for i in self.filled]}\n"
            string += f"   {[self.data[i] for i in self.filled]}\n"
            string += f"  Indices of source row (collections)\n"
            string += f"   {self.predictors}\n"
            string += f"  Weights of source row (collections)\n"
            string += f"   {self.weights}\n"
            return string
    return Prediction()

# Generate a compact numpy structured array from the contents of self.
def to_struct(self):
    import time
    start = time.time()
    import numpy as np
    # Generate the dtypes
    dtypes = []
    for i,(n,t) in enumerate(zip(self.names,self.types)):
        # Update user on progress if too much time has elapsed..
        if (time.time() - start) > self.max_wait:
            print(f" {100.*i/len(self):.2f}% to struct", end="\r", flush=True)
            start = time.time()
        if (t != str):
            pair = (n,) + NP_TYPES[t]
        else:
            max_len = max(len(s) for s in self[n] if type(s) == str)
            pair = (n,) + NP_TYPES[t][:1]+(max_len,)
        dtypes.append(pair)
    dtypes = np.dtype(dtypes)
    # Return the numpy structured array
    return np.array([tuple(row) for row in self if (None not in row)], 
                    dtype=dtypes)

# Given a set of columns of this data to try and predict, run a
# k-fold cross validation over predictions using 'fill'.
def predict(self, target_columns, model=None, weights=None, bias={}, k=10):
    import numpy as np
    # Get the numeric representation of self.
    numeric = self.to_matrix()
    # Pick the filling model based on data size.
    if (model is None):
        from util.approximate import unique, condition
        samples = min(numeric.data.shape[0], 10000)
        dim = min(numeric.data.shape[1], 1000)
        # Use nearest neighbor (can't get points + weights from NN)
        from util.approximate import NearestNeighbor
        model = condition(unique(NearestNeighbor), dim=dim, samples=samples)()
    # Set the weights accordingly
    if (weights is None):
        weights = np.ones(numeric.data.shape[1])
    elif (len(weights) != self.shape[1]):
        raise(Data.BadValue("Weights vector is length {len(weights)}, expected length {self.shape[1]}."))
    else:
        # Convert the provided weight to the correct shape.
        col_weights, weights = weights, []
        col_inds = list(range(len(numeric.names)))
        while (len(col_inds) > 0):
            if (col_inds[0] not in numeric.cats):
                col_inds.pop(0)
                weights.append( col_weights.pop(0) )
            else:
                for i in range(1,len(col_inds)):
                    not_categorical = (col_inds[i] not in numeric.cats)
                    beginning_of_next = ("1" == 
                        numeric.names[col_inds[i]].split('-')[-1])
                    if (not_categorical or beginning_of_next): break
                else: i += 1
                col_inds = col_inds[i:]
                weights += [col_weights.pop(0)] * i
        weights = np.array(weights)
    # Make sure the targets are in a list
    if type(target_columns) == str: 
        target_columns = [target_columns]
    elif type(target_columns) != list: 
        raise(Data.BadTarget("Target must either be a string (column name) or list of strings (column names)."))
    # Get the data all normalized (and weighted)
    normalized_data = weights * (numeric.data - numeric.shift) / numeric.scale
    print("normalized_data.shape: ",normalized_data.shape)
    # Identify those columns that are being predicted
    target_column_indices = [self.names.index(c) for c in target_columns]
    print("target_column_indices: ",target_column_indices)
    results = Data(names=[c + " Predicted" for c in target_columns]+["Prediction Index"],
                   types=[self.types[i] for i in target_column_indices]+[int])
    # Identify the numeric ranges of those columns
    indices = set(target_column_indices)
    print("target names: ",[self.names[i] for i in target_column_indices])
    sample = numeric.to_real([
        (v if i in indices else None)
        for i,v in enumerate(self[0])])
    del(indices)
    train_cols = np.array([i for (i,v) in enumerate(sample) if (v is None)])
    test_cols =  np.array([i for (i,v) in enumerate(sample) if (v is not None)])
    print("train_cols: ",train_cols)
    print("test_cols: ",test_cols)
    del(sample)
    # Cycle through and make predictions.
    for i,(train_rows, test_rows) in enumerate(
            Data.k_fold(range(numeric.data.shape[0]), k=k)):
        print(f"  Fold {i+1} of {k}.. ", end="", flush=True)
        train_data = normalized_data[train_rows]
        print(f"fitting {i+1} of {k}.. ", end="", flush=True)
        model.fit( normalized_data[train_rows,:][:,train_cols], 
                   normalized_data[train_rows,:][:,test_cols]   )
        print(f"predicting {i+1} of {k}.. ", end="", flush=True)
        predictions = model.predict( normalized_data[test_rows,:][:,train_cols] )
        # De-normalize predictions
        predictions = ( (predictions / weights[test_cols]) * 
                        numeric.scale[test_cols] + 
                        numeric.shift[test_cols] )
        print(f"storing.. {i+1} of {k}", end="\r", flush=True)
        # Store the prediction results.
        for j,i in enumerate(test_rows):
            row = list(numeric.data[i,train_cols]) + list(predictions[j,:])
            row = numeric.from_real(row, bias=bias)
            results.append( [row[i] for i in target_column_indices] + [i] )
        print(" "*70,end="\r",flush=True)
    # Sort results by the original index
    results.sort(key=lambda r: r[-1])
    results.pop("Prediction Index")
    # Return the results in a new Data object.
    return results


# Verify conversion into numpy structured array (without "None")
b = a.to_struct()
assert(type(b) == np.ndarray)
assert(type(b.dtype.names) == tuple)
assert(tuple(b['1']) == tuple(a['1'][:-1]))


# Test the 'fill' method.
b = a[:]
b[-2,1] = 'd'
b.append( b[0] )
assert(str(b.fill(b[-2]).data) == str([-1,'a',1.7999999999999998]))
assert(str(a.fill(a[-1][:]).data) == str([-1,'2',2.1]))

# Test cases for using "fill" with weights.
b.pop(-2)
b["3"] = b["2"]
b["2"] = b["1"]
b[-1,-2] = 'b'
b[-1,-1] = None
condense = lambda v: str(v)[:4]
assert(tuple(map(condense, b.fill(b[-1], weights=[1., 2., 3., 4.]).data))
       == ('1', 'a', 'b', '2.50'))
b[-1,-1] = None
assert(tuple(map(condense, b.fill(b[-1], weights=[100., 10., 1., 1000.]).data))
       == ('1', 'a', 'b', '1.21'))

