import numpy as np

# This subclass of a NumPy matrix has relevant information to
# the process used when transforming a Data object into a
# purely numerical matrix of values for prediction purposes.
# 
# Numeric.to_real   -- the function for processing regular
#                      vectors into associated real vectors
# Numeric.from_real -- the function for processing real vectors
#                      back into regular vector
# Numeric.in_names  -- names of columns before being mapped
# Numeric.in_types  -- types of columns before being mapped
# Numeric.mappings  -- the mappings of input column values to outputs
# Numeric.groups    -- the collections of real columns tied to input columns
# Numeric.dropped   -- columns that are dropped from input
# Numeric.names     -- meaningful column names for the real vectors
# Numeric.nums      -- indices of numerical value outputs
# Numeric.cats      -- indices of categorical mapping outputs
# Numeric.shift     -- How to shift this data to have min 0.
# Numeric.scale     -- How to scale this data (after shift) to 
#                      have a domain of [0., 1.]                    
class Numeric(np.ndarray):

    # ----------------------------------------------------------------
    # Initialize this Numeric object with special Numpy ndarray
    # subclassing methods (these are needed for correct behavior).
    def __new__(cls, data, in_names, in_types, mappings, groups,
                dropped, names, nums, cats, print_column_width):
        # Either an actual matrix or a shape can be provided.
        if (type(data) is tuple):
            self = np.asarray(np.zeros(data)).view(cls)
        else:
            self = np.asarray(data).view(cls)
        # Make sure the shape is correct.
        assert (len(names) == self.shape[1]), f"Data number of columns ({self.shape[1]}) must equal length of names ({len(names)})."
        # Assign extra local information.
        self.in_names = in_names
        self.in_types = in_types
        self.groups = groups
        self.mappings = mappings
        self.dropped = dropped
        self.names = names
        self.nums = nums
        self.cats = cats
        self.print_column_width = print_column_width
        # Store an internal lookup table that maps back to original values.
        self.np_mappings = {}
        for col in mappings:
            values = []
            encodings = []
            for (val,encoding) in mappings[col].items():
                values.append(val)
                encodings.append(encoding)
            self.np_mappings[col] = (np.asarray(encodings), values)
        # Compute the shift and scale for normalization (of numeric values)
        self.shift = np.zeros(self.shape[1])
        self.scale = np.ones(self.shape[1])
        if (type(data) is not tuple):
            self.compute_shift_and_scale()
        return self
    def __array_finalize__(self, obj):
        if obj is None: return
        # Assign extra local information.
        self.in_names = getattr(obj, "in_names", None)
        self.in_types = getattr(obj, "in_types", None)
        self.groups = getattr(obj, "groups", None)
        self.mappings = getattr(obj, "mappings", None)
        self.dropped = getattr(obj, "dropped", None)
        self.names = getattr(obj, "names", None)
        self.nums = getattr(obj, "nums", None)
        self.cats = getattr(obj, "cats", None)
        self.print_column_width = getattr(obj, "print_column_width", None)
        self.np_mappings = getattr(obj, "np_mappings", None)
        # Compute the shift and scale for normalization (of numeric values)
        self.shift = getattr(obj, "shift", None)
        self.scale = getattr(obj, "scale", None)
    # Define a reduce function that includes internal attributes.
    def __reduce__(self):
        cls, internal, state = super().__reduce__()
        state = state + (
            self.in_names, self.in_types, self.groups, self.mappings,
            self.dropped, self.names, self.nums, self.cats,
            self.print_column_width, self.np_mappings, self.shift, self.scale)
        return (cls, internal, state)
    # Define a setstate function that captures internal attributes.
    def __setstate__(self, state):
        np_state, custom_state = state[:-12], state[-12:]
        self.in_names, self.in_types, self.groups, self.mappings, \
            self.dropped, self.names, self.nums, self.cats, \
            self.print_column_width, self.np_mappings, self.shift, \
            self.scale = custom_state
        super().__setstate__(np_state)
    # ----------------------------------------------------------------

    # Update the calculation of the shift and scale for this Numeric object.
    def compute_shift_and_scale(self):
        self.shift[self.nums] = np.min(self[:,self.nums], axis=0)
        self.scale[self.nums] = np.max(self[:,self.nums], axis=0) - self.shift[self.nums]

    # Takes a row in the original space and generates a row in the
    # associated real-valued space.
    def to_real(self, row, fill=False):
        # Check for row length
        if (len(row) != len(self.in_names)):
            raise(Data.BadElement(f"Provided row has {len(row)} elements, expected {len(names)}."))
        # Bulid the real-valued vector
        real_row = []
        for i,(n,t,v) in enumerate(zip(self.in_names, self.in_types, row)):
            # Skip dropped columns.
            if (i in self.dropped): continue
            # If there is no map, then we must exit here.
            if (n not in self.mappings):
                from .exceptions import BadValue
                raise(Data.BadValue(f"No known mapping for column '{n}' to real values."))
            n_map = self.mappings[n]
            # Check if this incoming value is missing.
            if (v is None):
                if (None not in n_map):
                    from .exceptions import BadValue
                    raise(BadValue(f"Value '{v}' for column '{n}' was unexpected and cannot be mapped to real."))
                real_row += n_map[None]
            # Check if the type if float or int (literally translate).
            elif (t is float):
                real_row += [v]
            elif (t is int):
                real_row += [float(v)]
            # Otherwise, this value needs to be mapped.
            else:
                # If the value is not known..
                if (v not in n_map):
                    from .exceptions import BadValue
                    raise(BadValue(f"Value '{v}' for column '{n}' has no known real encoding."))
                real_row += n_map[v]
        # Check for row length
        if (len(real_row) != len(self.names)):
            from .exceptions import BadElement
            raise(BadElement(f"Provided row has {len(row)} elements, translation to {len(self.names)} unsuccesful, only created {len(real_row)}."))
        # Return the real-valued vector
        return real_row


    # Takes a row in the real-valued space and generates a row in the
    # original space.
    def from_real(self, real_row, bias={}):
        import numpy as np
        # Verify length of real-vector
        real_row = list(real_row)
        if (len(real_row) != len(self.names)):
            from .exceptions import BadElement
            raise(BadElement(f"Provided real vector has {len(real_row)} elements, expected {len(self.names)}."))
        # Build the original row
        row = []
        for (n,t) in zip(self.in_names, self.in_types):
            cols = self.groups[n]
            if (t is float):
                row += [real_row[cols[0]]]
            elif (t is int):
                row += [int(round(real_row[cols[0]]))]
            else:
                encoding = np.asarray([real_row[c] for c in cols])
                encodings, values = self.np_mappings[n]
                nearest = np.argmin(np.linalg.norm(encodings - encoding, axis=1))
                row += [values[nearest]]
        # Check for row length
        if (len(row) != len(self.in_names)):
            from .exceptions import BadElement
            raise(BadElement(f"Provided row has {len(real_row)} elements, translation to {len(self.in_names)} unsuccesful, only created {len(row)}."))
        return row


    # Produce a string that includes semantic information about this object.
    def __str__(self):
        opener = "  ["
        str_real_names = [opener]
        for n in self.names:
            if not ((len(str_real_names[0]) == len(opener)) or
                    (len(str_real_names[-1]) + len(n) < self.print_column_width)):
                str_real_names.append("   ")
            str_real_names[-1] += f"'{n}', "
        str_real_names[-1] = str_real_names[-1][:-2] + "]"
        str_real_names = "\n".join(str_real_names)
        return f"Numeric object for data with {len(self.names)} columns named:\n{str_real_names}\n\n{super().__str__()}"


# import inspect
# # Assign the comments to the official documentation spot.
# Numeric.__doc__ = inspect.getcomments(Numeric).replace("# ","")
