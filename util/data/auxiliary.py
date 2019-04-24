from util.data.exceptions import BadValue, ImproperUsage

# This class is a protected type of list where the contents are
# supposed to be of a specific type. It enforces that assignment and
# append operations do not violate the type specifiation of the object.
class Descriptor:
    # Initialize this descriptor, make sure all elements of 'values'
    # have the same type if provided.
    def __init__(self, data=None, t=type(None), indices=None):
        # Store the provided data (or make a new set)
        if (type(data) == type(None)): data = []
        self.data = list(data)
        # Store the provided indices (or make a new list of them)
        if (type(indices) == type(None)): indices = list(range(len(data)))
        self.indices = indices
        # Initialize the type.
        self.type = t
        for value in self.data:
            if (self.type == type(None)): self.type = type(value)
            if (self.type != type(value)):
                raise(BadValue(f"Only type '{self.type}' objects can be added to this descriptor."))

    # Check for the type of the added value.
    def append(self, value):
        if (self.type == type(None)): self.type = type(value)
        if (self.type != type(value)):
            raise(BadValue(f"Only type '{self.type}' objects can be added to this descriptor."))
        return self.data.append(value)

    # Check for type in the "setitem" operation.
    def __setitem__(self, index, value):
        if (self.type == type(None)): self.type = type(value)
        if (self.type != type(value)):
            raise(BadValue(f"Only type '{self.type}' objects can be added to this descriptor."))
        self.data[self.indices[index]] = value

    # When "getitem" is called, make a "view" over this Descriptor.
    def __getitem__(self, index):
        # If this is a slice, then return a 
        if (type(index) == slice):
            return type(self)(data=self.data, t=self.type, indices=self.indices[index])
        # Use the standard getitem operator.
        return self.data[self.indices[index]]

    # Define a "len" operator.
    def __len__(self): return len(self.data)

    # Define a "str" operator.
    def __str__(self): return str(self.data)

    # Define an "insert" function.
    def insert(self, idx, value):
        if (self.type == type(None)): self.type = type(value)
        if (self.type != type(value)):
            raise(BadValue(f"Only type '{self.type}' objects can be added to this descriptor."))
        if (tuple(set(self.indices)) != tuple(range(len(self.data)))):
            raise(ImproperUsage("Not allowed to insert into a Descriptor view."))
        # Add the additional index.
        self.indices.append(len(self.data))
        # Add the new "value" at a desired "index".
        return self.data.insert(self.indices[idx], value)

    # Identify the index (column number) of a value.
    def index(self, value): return self.indices.index(self.data.index(value))

    # Pop a value out of this descriptor.
    def pop(self, index):
        self.indices.pop(index)
        return self.data.pop(index)
        

# Local mutable Column class (for column item assignment,
# retrieval, and iteration over a column in a data object).
class Column:
    index = 0
    def __init__(self, data, column, indices=None):
        # Generate indices if they were not provided
        if (type(indices) == type(None)):
            indices = list(range(len(data)))
        # Store the parent data and the column number 
        self.data = data
        self.column = column
        self.indices = indices

    #     Indexing
    # -----------------
    def __setitem__(self, index, value):
        index = self.indices[index]
        # Default setitem for the parent data
        self.data[index, self.column] = value
    def __getitem__(self, index):
        # Default getitem for the parent data
        if (type(index) == int):
            index = self.indices[index]
            return self.data[index, self.column]
        elif (type(index) == slice):
            indices = range(len(self.indices))[index]
        elif (hasattr(index, "__iter__")):
            # Convert provided index into a list of integers.
            indices = []
            for i in index:
                if (type(i) != int):
                    raise(self.data.BadIndex(f"Data column index iterable contained type {type(i)}, expected {int}."))                    
                indices.append( self.indices[i] )
        else:
            raise(self.data.BadIndex(f"Data column does not support indexing with type {type(index)}."))
        # Return a new column with modified indices.
        return self.data.Column(self.data, self.column, indices)

    #     Iterator
    # ----------------
    # Define this as an interator
    def __iter__(self):
        self.index = 0
        return self
    def __next__(self):
        if (self.index >= len(self.indices)):
            raise(StopIteration)
        self.index += 1
        return self.data[self.indices[self.index-1], self.column]

    #    Descriptors
    # -----------------
    # Get the length of this Column
    def __len__(self): return len(self.indices)
    # Use the iteration to generate a string of this column
    def __str__(self):
        return str(list(self))

    #     Operators
    # -----------------
    # Generic boolean operator function.
    def _gen_operator(self, other, op_name):
        # First check for length.
        if hasattr(other, "__len__"):
            # If it is not a singleton, raise an error.
            if (type(other) != self.data.types[self.column]):
                # If the length doesn't match, check if it is a singleton.
                if (len(other) != len(self)):
                    raise(self.data.BadValue(f"Length '{len(other)}' does not have the same number of elements as this Column, {len(self)}."))
                # Return pairwise comparison for equal-length objects.
                for (v1, v2) in zip(self, other):
                    if (type(v1) == type(None)) or (type(v2) == type(None)): yield None
                    else:                                                    yield getattr(v1, op_name)(v2)
            else:
                # Otherwise return singleton comparison.
                for v in self:
                    if (type(v) == type(None)): yield None
                    else:                       yield getattr(v, op_name)(other)
        # Without a length, it must be a singleton comparison.
        else:
            for v in self:
                if (type(v) == type(None)): yield None
                else:                       yield getattr(v, op_name)(other)

    # Custom equality operator.
    def _custom_equality(self, other, equality=True):
        # Try to convert given group to an iterable type for comparison.
        if ((type(other) == self.data.types[self.column]) or (type(other) == type(None))):
            # Return indices where single value equals column values.
            for i,v in enumerate(self):
                if (equality and (v == other)):       yield i
                elif (not equality) and (v != other): yield i
        elif (type(other) == type(self)):
            # Return indices where the two columns are equal (or not).
            for i,(v1, v2) in enumerate(zip(self, other)):
                if (equality and (v1 == v2)):       yield i
                elif (not equality) and (v1 != v2): yield i
        elif (type(other) == set):
            # Iterate over rows, generating indices that match the equality.
            for i,v in enumerate(self):
                if (equality and (v in other)):             yield i
                elif ((not equality) and (v not in other)): yield i
        else:
            raise(self.data.BadData(f"There is no defined equality operator for the provided type {type(other)}, when it does not match expected type {self.data.types[self.column]}."))

    # Define the custom inequality operator using custom equality.
    def _custom_inequality(self, other): return self._custom_equality(other, equality=False)

    # Define the "in" method to return True for contained values.
    def __in__(self, value): return any(v == value for v in self)
    # Equality operator.
    def __eq__(self, other): return self._custom_equality(other)
    # Inequality operator.
    def __ne__(self, other): return self._custom_inequality(other)
    # Less than operator.
    def __lt__(self, other): return self._gen_operator(other, "__lt__")
    # Greater than operator.
    def __gt__(self, other): return self._gen_operator(other, "__gt__")
    # Less than equal operator.
    def __le__(self, other): return self._gen_operator(other, "__le__")
    # Greater than equal operator.
    def __ge__(self, other): return self._gen_operator(other, "__ge__")
    # Addition operator.
    def __add__(self, other): return self._gen_operator(other, "__add__")
    # Subtraction operator.
    def __sub__(self, other): return self._gen_operator(other, "__sub__")
    # Right-side addition operator.
    def __radd__(self, other): return self._gen_operator(other, "__radd__")
    # Right-side subtraction operator.
    def __rsub__(self, other): return self._gen_operator(other, "__rsub__")
    # Multiplication operator.
    def __mul__(self, other): return self._gen_operator(other, "__mul__")
    # True division operator ( a / b ).
    def __truediv__(self, other): return self._gen_operator(other, "__truediv__")
    # Floor division operator ( a // b ).
    def __floordiv__(self, other): return self._gen_operator(other, "__truediv__")


# Local class for iterating over 'rows' of this data object,
# prevents each row from being modified directly by users.
class Row:
    def __init__(self, data, values):
        self.data = data
        self.values = values

    # Convert an acceptable index into a list of numbers.
    def _index_to_nums(self, index):
        # Declare the indices of this row based on its parent data object.
        if (type(self.data.col_indices) != type(None)): indices = self.data.col_indices
        else:                                           indices = list(range(len(self.values)))
        # Return integer indices.
        if   type(index) == int:   nums = [indices[index]]
        # Return slice indices.
        elif type(index) == slice: nums = [indices[i] for i in range(len(self.values))[index]]
        # Return string indices (based on parent data).
        elif type(index) == str:   nums = [self.data.names.index(index)]
        # Convert an iterable into a list of numeric indices.
        elif hasattr(index, "__iter__"):
            nums = []
            for i,v in enumerate(index):
                if   type(v) == int:  nums.append(indices[v])
                elif type(v) == str:  nums.append(self.data.names.index(v))
                elif (type(v) == bool) and v: nums.append(indices[i])
                else: raise(self.data.BadIndex(f"Index '{index}' not understood by {type(self)} object."))
            # Make sure that boolean indices are full length.
            if (type(v) == bool) and (i != len(self)-1):
                raise(self.data.BadValue(f"Length '{i+1}' boolean index does not have the same number of elements as this Row, {len(self)}."))
        else: raise(self.data.BadIndex(f"Index '{index}' not understood by {type(self)} object."))
        return nums
    # Get an item from this row.
    def __getitem__(self, index):
        idx = self._index_to_nums(index)
        if (len(idx) == 1): return self.values[idx[0]]
        else:               return [self.values[i] for i in idx]
    # Set an item in this row.
    def __setitem__(self, index, value):
        idx = self._index_to_nums(index)
        if (len(idx) > 1): raise(self.data.BadAssignment(f"{type(self)} object can only assign one position at a time."))
        i = idx[0]
        # Assign the type of the column of data if it is unassigned.
        if self.data.types[i] == type(None):
            self.data.types[i] = type(value)
        elif type(value) == type(None): pass
        elif (self.data.types[i] != type(value)):
            raise(self.data.BadValue(f"Index {i} can only accept {self.data.types[i]}, received {type(value)}."))
        # Assign the value if it is allowable.
        self.values[i] = value
    # Define a "length" for this row.
    def __len__(self): return min(len(self.data.names), len(self.values))
    def __str__(self):
        # Declare the indices of this row based on its parent data object.
        if (type(self.data.col_indices) != type(None)): indices = self.data.col_indices
        else:                                           indices = list(range(len(self.values)))
        return str([self.values[i] for i in indices])
    # Define a "insert" function.
    def insert(self, i, value):
        # Declare the indices of this row based on its parent data object.
        if (type(self.data.col_indices) != type(None)): indices = self.data.col_indices
        else:                                           indices = list(range(len(self.values)))
        if (tuple(set(indices)) != tuple(range(len(self.values)))):
            raise(self.data.ImproperUsage(f"This row is only a view and does not support insertion."))
        if (len(self) != len(self.data.names)-1):
            print()
            print("self.data.names: ",len(self.data.names))
            print("self.data.names: ",self.data.names)
            print("len(self):       ",len(self))
            print("self.values:     ",self.values)
            raise(self.data.ImproperUsage(f"Invalid insertion operation on {type(self)}."))
        # Return the insertion of the new value.
        return self.values.insert(i, value)
    # Define a "pop" function.
    def pop(self, i):
        # Declare the indices of this row based on its parent data object.
        if (type(self.data.col_indices) != type(None)): indices = self.data.col_indices
        else:                                           indices = list(range(len(self.values)))
        if (tuple(set(indices)) != tuple(range(len(self.values)))):
            raise(self.data.ImproperUsage(f"This row is only a view and does not support 'pop'."))
        if (len(self) != len(self.data.names)):
            raise(self.data.ImproperUsage(f"Invalid pop operation on {type(self)}."))
        return self.values.pop(i)

    #    Operators
    # ---------------

    # Generic boolean operator function.
    def _gen_operator(self, other, op_name):
        # First check for length.
        if hasattr(other, "__len__"):
            # If the length doesn't match, check if it is a singleton.
            if (len(other) != len(self)):
                # If it is not a singleton, raise an error.
                if (type(other) != self.data.types[self.column]):
                    raise(self.data.BadValue(f"Length '{len(other)}' does not have the same number of elements as this Row, {len(self)}."))
                # Otherwise return singleton comparison.
                for v in self:
                    if (type(v) == type(None)): yield None
                    else:                       yield getattr(v, op_name)(other)
            # Return pairwise comparison for equal-length objects.
            else:
                for (v1, v2) in zip(self, other):
                    if (type(v1) == type(None)) or (type(v2) == type(None)): yield None
                    else:                                                    yield getattr(v1, op_name)(v2)
        # Without a length, it must be a singleton comparison.
        else:
            for v in self:
                if (type(v) == type(None)): yield None
                else:                       yield getattr(v, op_name)(other)

    # Inequality operator.
    def __ne__(self, other): return self._gen_operator(other, "__ne__")
    # Equality operator.
    def __eq__(self, other): return self._gen_operator(other, "__eq__")
    # Less than operator.
    def __lt__(self, other): return self._gen_operator(other, "__lt__")
    # Greater than operator.
    def __gt__(self, other): return self._gen_operator(other, "__gt__")
    # Less than equal operator.
    def __le__(self, other): return self._gen_operator(other, "__le__")
    # Greater than equal operator.
    def __ge__(self, other): return self._gen_operator(other, "__ge__")
    # Addition operator.
    def __add__(self, other): return self._gen_operator(other, "__add__")
    # Subtraction operator.
    def __sub__(self, other): return self._gen_operator(other, "__sub__")
    # Right-side addition operator.
    def __radd__(self, other): return self._gen_operator(other, "__radd__")
    # Right-side subtraction operator.
    def __rsub__(self, other): return self._gen_operator(other, "__rsub__")
    # Multiplication operator.
    def __mul__(self, other): return self._gen_operator(other, "__mul__")
    # True division operator ( a / b ).
    def __truediv__(self, other): return self._gen_operator(other, "__truediv__")
    # Floor division operator ( a // b ).
    def __floordiv__(self, other): return self._gen_operator(other, "__truediv__")
