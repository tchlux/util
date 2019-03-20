from util.data.exceptions import BadValue

# This class is a protected type of list where the contents are
# supposed to be of a specific type. It enforces that assignment and
# append operations do not violate the type specifiation of the object.
class Descriptor(list):
    # Initialize this descriptor, make sure all elements of 'values'
    # have the same type if provided.
    def __init__(self, values, t=type(None)):
        self.type = t
        super().__init__(values)
        for value in values:
            if (self.type == type(None)): self.type = type(value)
            if (self.type != type(value)):
                raise(BadValue(f"Only type '{self.type}' objects can be added to this descriptor."))

    # Check for the type of the added value.
    def append(self, value):
        if (self.type == type(None)): self.type = type(value)
        if (self.type != type(value)):
            raise(BadValue(f"Only type '{self.type}' objects can be added to this descriptor."))
        return super().append(value)

    # Check for type in the "setitem" operation.
    def __setitem__(self, index, value):
        if (self.type == type(None)): self.type = type(value)
        if (self.type != type(value)):
            raise(BadValue(f"Only type '{self.type}' objects can be added to this descriptor."))
        return super().__setitem__(index, value)


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
    # Define the "in" method to return True for contained values.
    def __in__(self, value): return any(v == value for v in self)
    # Inequality operator
    def __ne__(self, group):
        return self.__eq__(group, equality=False)
    # Iterate over the indices that exist in a set of values
    def __eq__(self, group, equality=True):
        # Try to convert given group to an iterable type for comparison.
        if ((type(group) == self.data.types[self.column]) or (type(group) == type(None))):
            # Return indices where single value equals column values.
            for i,v in enumerate(self):
                if (equality and (v == group)):       yield i
                elif (not equality) and (v != group): yield i
        elif (type(group) == type(self)):
            # Return indices where the two columns are equal (or not).
            for i,(v1, v2) in enumerate(zip(self, group)):
                if (equality and (v1 == v2)):       yield i
                elif (not equality) and (v1 != v2): yield i
        elif (type(group) in {set, list}):
            # Iterate over rows, generating indices that match the equality.
            for i,v in enumerate(self):
                if (equality and (v in group)):             yield i
                elif ((not equality) and (v not in group)): yield i
        else:
            raise(self.data.BadData(f"There is no defined equality operator for the provided type {type(group)}, when it does not match expected type {self.data.types[self.column]}."))
    # Iterate over the indices based on a comparison operator
    def __lt__(self, val):
        # Iterate over values in this column
        for i,v in enumerate(self):
            if (v < val): yield i
    # Iterate over the indices based on a comparison operator
    def __gt__(self, val):
        # Iterate over values in this column
        for i,v in enumerate(self):
            if (v > val): yield i
    # Iterate over the indices based on a comparison operator
    def __le__(self, val):
        # Iterate over values in this column
        for i,v in enumerate(self):
            if (v <= val): yield i
    # Iterate over the indices based on a comparison operator
    def __ge__(self, val):
        # Iterate over values in this column
        for i,v in enumerate(self):
            if (v >= val): yield i


# Local class for iterating over 'rows' of this data object,
# prevents each row from being modified directly by users.
class Row:
    def __init__(self, data, values):
        self.data = data
        self.values = values
    # Convert an acceptable index into a list of numbers.
    def _index_to_nums(self, index):
        # Return integer indices.
        if   type(index) == int:   nums = [index]
        # Return slice indices.
        elif type(index) == slice: nums = range(len(self.values))[index]
        # Return string indices (based on parent data).
        elif type(index) == str:   nums = [self.data.names.index(index)]
        # Convert an iterable into a list of numeric indices.
        elif hasattr(index, "__iter__"):
            nums = []
            for v in index:
                if   type(v) == int: nums.append(v)
                elif type(v) == str: nums.append(self.data.names.index(v))
                else: raise(self.data.BadIndex(f"Index '{index}' not understood by {type(self)} object."))
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
        elif type(value) == type(None):
            self.data.missing.add(id(self))
        elif (self.data.types[i] != type(value)):
            raise(self.data.BadValue(f"Index {i} can only accept {self.data.types[i]}, received {type(value)}."))
        # Assign the value if it is allowable.
        self.values[i] = value
        # Remove this row from the "missing" list if all values are assigned.
        if (id(self) in self.data.missing) and (None not in self):
            self.data.missing.remove(id(self))
    # Define a "length" for this row.
    def __len__(self): return len(self.values)
    def __str__(self): return str(self.values)
    # Define a "insert" function.
    def insert(self, i, value):
        if (len(self) != len(self.data.names)-1):
            raise(self.data.ImproperUsage(f"Invalid insertion operation on {type(self)}."))
        return self.values.insert(i, value)
    # Define a "pop" function.
    def pop(self, i):
        if (len(self) != len(self.data.names)+1):
            raise(self.data.ImproperUsage(f"Invalid pop operation on {type(self)}."))
        return self.values.pop(i)

    # #    Operators
    # # ---------------
    # # Define the "in" method to return True for contained values.
    # def __in__(self, value): return any(v == value for v in self.values)
    # # Inequality operator
    # def __ne__(self, group):
    #     return self.__eq__(group, equality=False)
    # # Iterate over the indices that exist in a set of values
    # def __eq__(self, group, equality=True):
    #     # Try to convert given group to an iterable type for comparison.
    #     if (not hasattr(group, "__len__")):
    #         raise(self.data.ImproperUsage(f"Type '{type(group)}' cannot be compared with Row."))
    #     elif (len(group) != len(self.values)):
    #         raise(self.data.BadValue(f"Length '{len(group)}' does not have the same number of elements as this Row."))
    #     elif (type(group) == type(self)):
    #         # Return indices where the two columns are equal (or not).
    #         for i,(v1, v2) in enumerate(zip(self, group)):
    #             if (equality and (v1 == v2)):       yield i
    #             elif (not equality) and (v1 != v2): yield i
    #     elif (type(group) in {set, list}):
    #         # Iterate over rows, generating indices that match the equality.
    #         for i,v in enumerate(self):
    #             if (equality and (v in group)):             yield i
    #             elif ((not equality) and (v not in group)): yield i
    #     else:
    #         raise(self.data.BadData(f"There is no defined equality operator for the provided type {type(group)}, when it does not match expected type {self.data.types[self.column]}."))
    # # Iterate over the indices based on a comparison operator
    # def __lt__(self, val):
    #     # Iterate over values in this column
    #     for i,v in enumerate(self):
    #         if (v < val): yield i
    # # Iterate over the indices based on a comparison operator
    # def __gt__(self, val):
    #     # Iterate over values in this column
    #     for i,v in enumerate(self):
    #         if (v > val): yield i
    # # Iterate over the indices based on a comparison operator
    # def __le__(self, val):
    #     # Iterate over values in this column
    #     for i,v in enumerate(self):
    #         if (v <= val): yield i
    # # Iterate over the indices based on a comparison operator
    # def __ge__(self, val):
    #     # Iterate over values in this column
    #     for i,v in enumerate(self):
    #         if (v >= val): yield i
