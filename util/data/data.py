import os
from util.data.read import read_data
from util.data.categories import regular_simplex
from util.data import DEFAULT_DISPLAY_WAIT, GENERATOR_TYPE, NP_TYPES


# TODO:  Make saving to CSV files automatically put quotes around
#        strings that have the separator character in them.
# TODO:  Make printout shorter for data with lots of columns.
# TODO:  Add tests for empty data object.
# TODO:  Check the __iadd__ for data object columns, breaking "missing".
# TODO:  Update __iadd__ to behave differently when adding data:
#        when adding data with subset of known columns, add rows and
#        columns (for unknown) and None for missing.
# TODO:  Make "to_matrix" automatically flatten elements of columns
#        that contain iterables into their base types. If the types
#        are numeric, use float, otherwise use string. Use "is numeric"
#        to determine if the target needs to be classified or interpolated.
# TODO:  Make "to_matrix" automatically consider 'None' a unique value
#        for categorical columns and add a dimension 'is present' otherwise.
# TODO:  Remove the default storage of "to_matrix" output.
# TODO:  Wrap the necessary "Numeric" attributes into the matrix,
#        using a custom ndarray object from numpy.
# TODO:  Make data only deep copy when "copy()" is called or "[:,:]"
#        is given. Otherwise return "DataView" object (limited class 
#        that allows for deep copying if desired).
# TODO:  Write "collect" that does the combined operations of 
#        "unique().collect()" in a more memory efficient manner (by
#        cannabalizing the existing Data object).
# TODO:  Write "expand" that does the opposite operation of "collect",
#        expanding a data set out according to list columns.
# TODO:  Check for error in "collect" when column names are
#        duplicated. Raise errror for duplicate column names?
# TODO:  Speed up the read of a file. Put the read operation into the
#        background? Store everything as a string until the type is
#        determined, then update the types at the end.
# TODO:  Rewrite "predict" to work given a data object. This will
#        break up the data into collections of unique missing values.
#        Then for each collection, build a model over the appropriate
#        matrix-format data and predict that collection. This will
#        work more nicely with the k-fold method.
# TODO:  Remove the "fill" method in favor of "predict".

# TODO:  Make "predict" not require that target columns be convertible
#        to numeric form, only that they operate under a weighted sum.
# TODO:  Implment tracking for the "size" of a Row & Data object, and
#        the recently accessed rows (keep an ordered set of them).
# TODO:  Implement caching for rows when a Data object becomes very large.
# TODO:  Make reading data asynchronous, so waiting is only ever done
#        when an object is accessed by index.

# TODO:  Implmenent 'IndexedAudio' reader that mimics the
#        'IndexedVideo' reader and allows for access to peaks.
# TODO:  Implemenet 'IndexedText' reader that mimics the
#        'IndexedVideo' reader and allows for access to characters.

# TODO:  Subclass Data object so that the methods are broken up by the
#        type of function they provide. Stats should be an additional
#        set of methods outside of summary / operators.



# ======================================================
#      Python 'Data' with named and typed columns     
# ======================================================

# Define a python data matrix structure (named and typed columns)
class Data:
    # Holder for data.
    data = None
    missing = None
    # Boolean for declaring if this is a "view" or actual data.
    row_indices = None
    col_indices = None
    # Default values for 'self.names' and 'self.types'
    _names = None
    _types = None
    # Default maximum printed rows of this data.
    max_display = 10
    # Default maximum string length of printed column values.
    max_str_len = 50
    # Define the maximum display width for a table (with many columns)
    # before the data is automatically displayed in different format.
    max_print_width = 100

    # Import all of the exceptions.
    from util.data.exceptions import NameTypeMismatch, \
        BadSpecifiedType, BadSpecifiedName, NoNamesSpecified, \
        ImproperUsage, UnknownName, Unsupported, BadAssignment, \
        BadElement, BadTarget, BadValue, BadIndex, BadData, Empty
    # Import the classes used to protect this data object.
    from util.data.auxiliary import Descriptor, Column, Row

    # Protect the "names" and "types" attributes by using properties
    # to control the access and setting of these values.
    @property
    def names(self): return self._names
    @names.setter
    def names(self, value): self._names = Data.Descriptor(value)
    @property
    def types(self): return self._types
    @types.setter
    def types(self, value): self._types = Data.Descriptor(value)
    # Return True if empty.
    @property
    def empty(self):
        return ((type(self.names) == type(None)) or 
                (type(self.types) == type(None)))
    # Declare the "start" property / function to re-initialized when called.
    @property
    def shape(self):
        if (type(self.col_indices) != type(None)): return (len(self), len(self.col_indices))
        else:                                      return (len(self), len(self.names))

    #      Overwritten list Methods     
    # ==================================

    def __init__(self, data=None, names=None, types=None, rows=None, cols=None):
        # Verify the length of the types and names
        if ((type(types) != type(None)) and
            (type(names) != type(None)) and
            (len(types) != len(names))):
            raise(self.NameTypeMismatch(f"Length of provided names {len(names)} does not match length of provided types {len(types)}."))
        # Get the names (must be strings)
        if (type(names) != type(None)):
            if any(type(n) != str for n in names):
                raise(self.BadSpecifiedName(f"An entry of provided names was not of {str} class."))
            self.names = list(names)
        # Store the types if provided
        if (type(types) != type(None)):
            if any(type(t) != type for t in types):
                raise(self.BadSpecifiedType(f"An entry of provided types was not of {type(type)} class."))
            self.types = list(types)
            # Construct default names for each column if names were not provided.
            if (type(self.names) == type(None)):
                self.names = [str(i) for i in range(len(types))]
        # Initialize the contents of this Data to be an empty list.
        self.data = []
        self.missing = set()
        # If this is a "view" object, store a pointer.
        if (type(rows) != type(None)) or (type(cols) != type(None)):
            self.data = data
            self.row_indices = rows
            self.col_indices = cols
        # Otherwise, build this data as a deep copy of provided data.
        elif type(data) != type(None):
            for row in data: self.append(row)
    # Overwrite the "len" operator.
    def __len__(self):
        if (type(self.row_indices) != type(None)): return len(self.row_indices)
        else:                                      return len(self.data)

    # Overwrite the standard "[]" get operator to accept strings as well.
    def __getitem__(self, index):
        # Identify whether or not this is a data view.
        view = ((type(self.row_indices) != type(None)) or (type(self.col_indices) != type(None)))
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            raise(self.Empty("Cannot get item from empty data."))
        if ((type(index) == tuple) and (len(index) == 2) and 
            (all(type(i) == int for i in index))):
            if view: index = (self.row_indices[index[0]], self.col_indices[index[1]])
            # This index is accessing a (row,col) entry
            return self.data[index[0]][index[1]]
        elif ((type(index) == tuple) and (len(index) == 2) and 
              (type(index[0]) == int and type(index[1]) == str)):
            if (index[1] not in self.names):
                raise(self.BadIndex(f"Column index {index[1]} is not a recognized column name."))
            # Update the index based on the subset of "row" and "column"
            # indices that are stored in this data (if it is a view).
            if view: index = (row_indices[index[0]], index[1])
            # This index is accessing a (row,col) entry
            return self.data[index[0]][self.names.index(index[1])]
        elif (type(index) == int):
            # Update the index based on the subset of "row"
            # indices that are stored in this data (if it is a view).
            if (type(self.row_indices) != type(None)): index = self.row_indices[index]
            # Return a view of a row.
            return self.data[index]
        elif (type(index) == str):
            # This index is accessing a column
            if (index not in self.names):
                raise(self.UnknownName(f"This data does not have a column named '{index}'."))
            col_number = self.names.index(index)
            # Return a mutable "Column" object
            return Data.Column(self, col_number)
        else:
            # This index is retrieving a sliced subset of self.
            rows, cols = self._index_to_rows_cols(index)
            # If they are extracting the full data then assume a copy is desired.
            if ((len(set(rows)), len(set(cols))) == self.shape): return self.copy()
            # Otherwise make an alias data set (a data "view") and return it.
            names = Data.Descriptor(self.names, self.names.type, cols)
            types = Data.Descriptor(self.types, self.types.type, cols)
            new_data = Data(data=self.data, names=names, types=types,
                            rows=rows, cols=cols)
            return new_data

    # Overwrite the standard "[]" get operator to accept strings as well.
    def __setitem__(self, index, value):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            raise(self.Empty("Cannot set item for empty data."))
        # Special case assignment of new column (*RETURN* to exit function).
        if ((type(index) == str) and (index not in self.names)):
            # Check to see if this is a data view object.
            if ((type(self.row_indices) != type(None)) or
                (type(self.col_indices) != type(None))):
                raise(self.ImproperUsage("This is a data alias and does not support assignment."))
            # Otherwise perform the normal column addition.
            return self.add_column(value, name=index)
        # Identify whether or not this is a data view.
        view = ((type(self.row_indices) != type(None)) or (type(self.col_indices) != type(None)))
        # Get the list of rows and list of columns being assigned.
        rows, cols = self._index_to_rows_cols(index)
        # Assume if it has an iterator that is not a singleton.
        # Assume that if the value type matches all assigned types it is a singleton.
        singleton = (not hasattr(value, '__iter__')) or \
                    (all(self.types[c] in {type(value), type(None)} for c in cols))
        # If it has a "length" then use that to verify singleton status.
        if (singleton and hasattr(value, "__len__")):
            singleton = len(value) != (len(rows)*len(cols))
        # Assume if the "value" is a string, then this is a singleton.
        if (not singleton) and (type(value) == str): singleton = True
        # Assume if the value is a generator and the column(s) it is
        # assigning to are not *all* generators, it is not a singleton.
        if ((type(value) == GENERATOR_TYPE) and
            any(self.types[c] != GENERATOR_TYPE for c in cols)): singleton = False
        # Assume singleton status is correct.
        if (not singleton): value_iter = value.__iter__()
        # Reset type entirely if the new data assigns to the full column,
        # only if this data views the entire column.
        if (len(rows) == len(self.data)):
            for col in cols: self.types[col] = type(None)
        # Iterate over rows and perform assignment.
        for step, row in enumerate(rows):
            if type(row) != int:
                raise(self.BadIndex(f"The provided row index of type {type(row)}, {row}, is not understood."))
            for col in cols:
                # Retreive next value from iterator if necessary.
                if not singleton:
                    try:    value = next(value_iter)
                    except: raise(self.BadValue(f"Provided iterable only contained {step+1} elements, expected more."))
                # The existing column doesn't have a type, assign it.
                if (self.types[col] == type(None)):
                    self.types[col] = type(value)
                # Mark this row as contianing a missing value if assigned 'None'.
                elif (type(value) == type(None)):
                    if (id(self[row]) not in self.missing):
                        self.missing.add(id(self[row]))
                # Check against the existing column type (for verification).
                elif (type(value) != self.types[col]):
                    raise(self.BadAssignment(f"Provided value {value} of type {type(value)} does not match expected type {self.types[col]}."))
                # Make the assignment
                self.data[row][col] = value
                # Remove entry from missing if appropriate.
                if (type(value) != type(None)):
                    if (id(self[row]) in self.missing):
                        if (None not in self[row]):
                            self.missing.remove(id(self[row]))
        # Verify that the iterable has been exhausted.
        if (not singleton):
            try:
                value = next(value_iter)
                raise(self.BadAssignment(f"Column assignment requires a column of length {len(self)}, provided values iterable was too long."))
            except StopIteration: pass

    # Printout a brief table-format summary of this data.
    def __str__(self, short=False):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            return "This Data has no contents.\n"
        # Custom short string for a 'data' type object.
        if short: return f"Data ({len(self.names)} x {len(self)}) -- {self.names}"
        # Identify whether or not this is a data view.
        view = ((type(self.row_indices) != type(None)) or (type(self.col_indices) != type(None)))
        # Make a pretty table formatted output
        num_rows, num_cols = self.shape
        rows = []
        rows += [ list(self.names) ]
        rows += [ list(map(lambda t: str(t)[8:-2],self.types)) ]
        for i,row in enumerate(self):
            if view: row = [row[col] for col in self.col_indices]
            rows += [ row ]
            if i >= self.max_display: break
        rows = [list(map(self._val_to_str,r)) for r in rows]
        lens = [max(len(r[c]) for r in rows)
                for c in range(num_cols)]
        rows = [[v + " "*(lens[i]-len(v)) for i,v in enumerate(row)]
                for row in rows]
        rows = [" " + (" | ".join(row)) + "\n" for row in rows]
        string = "\n" + "="*len(rows[0]) + "\n"
        string += f"Size: ({num_rows} x {num_cols})\n\n"
        string += rows[0]
        string += rows[1]
        string += "-"*len(rows[0]) + "\n"
        for row in rows[2:]:                string += row
        if (len(self) > self.max_display): string += "  ...\n"
        string += "\n"
        # Print some info about missing values if they exist.
        if len(self.missing) > 0:
            # Get the indices of missing elements
            indices = []
            for i in range(len(self)):
                if any(v == type(None) for v in map(type,self[i])): 
                    indices.append(i)
                # WARNING: This ^^ is inefficient and is being done
                #          because the "missing" id's are bugged.
                if len(indices) > self.max_display:
                    break
            # Print out the indices of the rows with missing values
            print_indices = ', '.join(map(str, indices[:self.max_display]))
            string += f" missing values ({len(self.missing)}) at rows: \n  [{print_indices}"
            # Add elipses if there are a lot of missing values.
            if (len(indices) > self.max_display): string += ", ..."
            string += "]\n"
        string += "="*len(rows[0]) + "\n"
        return string

    # Define a convenience funciton for concatenating another
    # similarly typed and named data to self.
    def __iadd__(self, data):
        # Check to see if this is a data view object.
        if ((type(self.row_indices) != type(None)) or
            (type(self.col_indices) != type(None))):
            raise(self.ImproperUsage("This is a data alias and does not support assignment."))
        # Check for improper usage
        if type(data) != Data:
            raise(self.Unsupported(f"In-place addition only supports type {Data}, but {type(data)} was given."))
        # Special case for being empty.
        if data.empty: return self
        # Special case for being empty.
        elif self.empty:
            self.names = data.names[:]
            self.types = data.types[:]
            # Load in the data
            for row in data:
                self.append(row)
        # Add rows to this from the provided Data
        elif all(n1 == n2 for (n1,n2) in zip(self.names, data.names)):
            # Load in the data
            for row in data:
                self.append(row)            
        # Add columns to this from the provided Data
        elif len(data) == len(self):
            for col in data.names:
                self[col] = data[col]
        # Return self
        return self

    # Define a convenience funciton for concatenating another
    # similarly typed and named data to self.
    def __add__(self, data):
        self = self.copy()
        self += data
        return self

    # Redefined append operator that 
    def append(self, element):
        # Check to see if this is a data view object.
        if ((type(self.row_indices) != type(None)) or
            (type(self.col_indices) != type(None))):
            raise(self.ImproperUsage("This is a data alias and does not support assignment."))
        # Convert the element into a python list
        try:    element = element
        except: raise(self.BadValue(f"Invalid appended element, failed conversion to list."))
        missing_values = False
        # Check length in case names already exist
        if (type(self.names) != type(None)):
            if (len(element) != len(self.names)):
                raise(self.BadElement(f"Only elements of length {len(self.names)} can be added to this data."))
            # Try type casing the new element
            for i, (val, typ) in enumerate(zip(element, self.types)):
                # Record 'None' types as missing entries
                if type(val) == type(None):
                    # This is a missing value
                    missing_values = True
                elif (typ == type(None)):
                    # Update unassigned types with the new values' type
                    self.types[i] = type(val)
                elif (type(val) != typ):
                    # If not missing, then check the type
                    try:
                        # Try casting the element as the expected type
                        element[i] = typ(val)
                    except ValueError:
                        # Otherwise raise an error to the user
                        raise(self.BadValue(f"Value '{val}' of type '{type(val)}' could not successfully be cast as '{typ}'."))
                    except TypeError:
                        print("Data error:", typ, typ==type(None))
        else:
            # Construct the expected "types" of this data
            self.types = [type(val) for val in element]
        if (type(self.names) == type(None)):
            # Construct default names for each column
            self.names = [str(i) for i in range(len(element))]
        # Call the standard append operator, adding the element to self
        self.data.append(Data.Row(self, element))
        # Add to the list of missing values if necessary
        if (missing_values):
            self.missing.add(id(self.data[-1]))


    # Make sure the copy of this data is a deep copy.
    def copy(self):
        # Identify whether or not this is a data view.
        view = ((type(self.row_indices) != type(None)) or (type(self.col_indices) != type(None)))
        # Construct a new data set that is a copy of the data in this object.
        if view: rows, cols = self.row_indices, self.col_indices
        else:    rows, cols = map(list,map(range,self.shape))
        data = [[self.data[row][col] for col in cols] for row in rows]
        # Return a new data object.
        return type(self)(data=data, names=list(self.names), types=list(self.types))

    # Overwrite the 'pop' method to track missing values.
    def pop(self, index=-1):
        # Identify whether or not this is a data view.
        view = ((type(self.row_indices) != type(None)) or (type(self.col_indices) != type(None)))
        if view: raise(self.ImproperUsage("Cannot 'pop' from a data view. Copy this object to remove items."))
        # Popping of a column
        if type(index) == int:
            if (id(self[index]) in self.missing):
                self.missing.remove( id(self[index]) )
            return self.data.pop(index)
        # Popping of a column
        elif (type(index) == str):
            if index not in self.names:
                raise(self.BadSpecifiedName(f"There is no column named {index} in this Data."))
            col = self.names.index(index)
            self.names.pop(col)
            self.types.pop(col)
            values = []
            for row in self:
                values.append( row.pop(col) )
                # Check for removal from missing values.
                if (id(row) in self.missing):
                    if all(type(v) != type(None) for v in row):
                        self.missing.remove(id(row))
            return values
        else:
            raise(self.BadIndex(f"Index {index} of type {type(index)} is not recognized."))

    # Overwrite the 'remove' method to track missing values.
    def remove(self, index):
        # Identify whether or not this is a data view.
        view = ((type(self.row_indices) != type(None)) or (type(self.col_indices) != type(None)))
        if view: raise(self.ImproperUsage("Cannot 'remove' from a data view. Copy this object to remove items."))
        # See if this row is in the "missing" set.
        if (id(self[index]) in self.missing):
            self.missing.remove( id(self[index]) )
        return self.data.remove( index )

    # ========================
    #      Custom Methods     
    # ========================

    # Custom function for mapping values to strings in the printout of self.
    def _val_to_str(self, v): 
        # Get the string version of the value.
        if not ((type(v) == type(self)) or (issubclass(type(v),type(self)))):
            string = str(v).replace("\n"," ")
        else:
            # Custom string for a 'data' type object.
            string = v.__str__(short=True)
        # Shorten the string if it needs to be done.
        if len(string) > self.max_str_len: 
            string = string[:self.max_str_len]+".."
        return string

    # Given an index (of any acceptable type), convert it into an
    # iterable of integer rows and a list of integer columns.
    def _index_to_rows_cols(self, index):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            raise(self.Empty("Cannot get rows and cols from empty data."))
        # Standard usage.
        if (type(index) == int):
            rows = [index]
            cols = list(range(len(self.names)))
        # Case for column index access.
        elif (type(index) == str):
            rows = list(range(len(self)))
            if index not in self.names:
                raise(self.UnknownName(f"This data does not have a column named '{index}'."))
            cols = [self.names.index(index)]
        # Case for slice-index access.
        elif (type(index) == slice):
            # Construct a new "data" with just the specified rows (deep copy)
            rows = list(range(len(self))[index])
            cols = list(range(len(self.names)))
        # Case for tuple-string index access where thera are two selected columns.
        elif ((type(index) == tuple) and (len(index) == 2) and
              (type(index[0]) == str) and (type(index[1]) == str)):
            rows = list(range(len(self)))
            cols = []
            for i in index:
                if i not in self.names:
                    raise(self.UnknownName(f"This data does not have a column named '{i}'."))
                cols.append( self.names.index(i) )
        # Case where (row indices, col indices) were provided.
        elif ((type(index) == tuple) and (len(index) == 2)):
            if all(type(i) == int for i in index):
                # This index is accessing a (row, col) entry.
                rows = [index[0]]
                cols = [index[1]]
            elif ((type(index[0]) == int) and (type(index[1]) == str)):
                # This index is accessing a (row, col) entry with a named column.
                rows = [index[0]]
                cols = [self.names.index(index[1])]
            else:
                rows, cols = index

                # Handling the "rows"
                if (type(rows) == int):     rows = [rows]
                elif (type(rows) == slice): rows = range(len(self))[rows]
                elif hasattr(rows, "__iter__"):
                    rows = list(rows)
                    if   all(type(i)==int for i in rows):  pass
                    elif all(type(i)==bool for i in rows):
                        rows = [i for (i,b) in enumerate(rows) if b]
                    else:
                        type_printout = sorted(map(str,set(map(type,rows))))
                        raise(self.BadIndex(f"The provided row index, {index}, is not understood.\n It has types {type_printout}."))
                else:
                    type_printout = sorted(map(str,set(map(type,rows))))
                    raise(self.BadIndex(f"The provided row index, {index}, is not understood.\n It is typed {type_printout}."))

                # Handling the "columns"
                if (type(cols) == int):     cols = [cols]
                elif (type(cols) == slice): cols = range(len(self))[cols]
                elif hasattr(cols, "__iter__"):
                    col_ids = list(cols)
                    cols = []
                    for i,v in enumerate(col_ids):
                        if (type(v) not in {str,int,bool}):
                            raise(self.BadIndex(f"The provided column index of type {type(v)}, {v}, is not understood."))
                        elif (type(v) == str):
                            if (v not in self.names):
                                raise(self.UnknownName(f"This data does not have a column named '{v}'."))
                            i = self.names.index(v)
                        # If it is an int, use that as the index.
                        elif (type(v) == int): i = v
                        # If it is a boolean and it is False, skip this index.
                        elif (type(v) == bool) and (not v): continue
                        # Append the column index
                        cols.append( i )
                else:
                    type_printout = sorted(map(str,set(map(type,cols))))
                    raise(self.BadIndex(f"The provided column index, {index}, is not understood.\n It has types {type_printout}."))
        # Iterable index access.
        elif hasattr(index, "__iter__"):
            index = list(index)
            # Case for iterable access of rows.
            if all(type(i)==int for i in index):
                rows = list(index)
                cols = list(range(len(self.names)))
            # Case for iterable access of columns.
            elif all(type(i)==str for i in index):
                rows = list(range(len(self)))
                cols = []
                for i in index:
                    if i not in self.names:
                        raise(self.UnknownName(f"This data does not have a column named '{i}'."))
                    cols.append( self.names.index(i) )
            # Case for boolean-index access to rows.
            elif all(type(i)==bool for i in index):
                rows = [i for (i,b) in enumerate(index) if b]
                cols = list(range(len(self.names)))
            # Undefined behavior for this index.
            else:
                type_printout = sorted(map(str,set(map(type,index))))
                raise(self.BadIndex(f"The provided index, {index}, is not understood.\n It is typed {type_printout}."))
        # Undefined behavior for this index.
        else:
            type_printout = sorted(map(str,set(map(type,index))))
            raise(self.BadIndex(f"The provided index, {index}, is not understood.\n It is typed {type_printout}."))
        # Return the final list of integer-valued rows and columns.
        return rows, cols

    # =======================================
    #      Saving and Loading from Files     
    # =======================================

    # Convenience method for loading from file.
    def load(path, *args, **read_data_kwargs):
        import pickle, dill, gzip
        # Handle two different usages of the 'load' function.
        if (type(path) != str): 
            # If the user called "load" on an initialized data, ignore
            self = path
            if (len(args) == 1):
                path = args[0]
            else:
                raise(Data.ImproperUsage("'load' method for Data must be given a path."))
        else:
            # Otherwise define 'self' as a new Data object.
            self = Data()
        # Check for extension, add default if none is provided.
        if "." not in path: path += ".pkl"
        # Check for compression
        compressed = path[-path[::-1].index("."):] == "gz"
        file_opener = open
        if compressed:
            file_opener = gzip.open
            base_path = path[:-path[::-1].index(".")-1]
            ext = base_path[-base_path[::-1].index("."):]
        else:
            ext = path[-path[::-1].index("."):]
        # Handle different base file extensions
        if (ext == "pkl"):
            with file_opener(path, "rb") as f:
                d = pickle.load(f)
                self += d
        elif (ext == "dill"):
            with file_opener(path, "rb") as f:
                d = pickle.load(f)
                # WARNING: Getting errors with old files..
                # print([list(r) for r in d][:10])
                # print(dir(d))
                self += d
        elif (ext in {"csv", "txt", "tsv"}):
            mode = "r" + ("t" if compressed else "")
            with file_opener(path, mode) as f:
                read_data_kwargs["opened_file"] = f
                self += read_data(**read_data_kwargs)
        else:
            raise(self.Unsupported(f"Cannot load file with extension '{ext}'."))
        return self

    # Convenience method for saving to a file
    def save(self, path, create=True):
        import pickle, dill, gzip
        # Create the output folder if it does not exist.
        output_folder = os.path.split(os.path.abspath(path))[0]
        if create and not os.path.exists(output_folder): os.makedirs(output_folder)
        # Check for compression
        if "." not in path: path += ".pkl"
        compressed = path[-path[::-1].index("."):] == "gz"
        file_opener = open
        if compressed:
            file_opener = gzip.open
            base_path = path[:-path[::-1].index(".")-1]
            ext = base_path[-base_path[::-1].index("."):]
        else:
            ext = path[-path[::-1].index("."):]
        # Handle different base file extensions
        if (ext == "pkl"):
            with file_opener(path, "wb") as f:
                pickle.dump(self, f)
        elif (ext == "dill"):
            with file_opener(path, "wb") as f:
                dill.dump(self, f)
        elif (ext in {"csv", "tsv"}):
            sep = "," if (ext == "csv") else "\t"
            mode = "w" + ("t" if compressed else "")
            with file_opener(path, mode) as f:
                print(sep.join(self.names), file=f)
                for row in self:
                    print_row = [str(v) if type(v) != type(None) else "" for v in row]
                    print(sep.join(print_row), file=f)
        else:
            raise(self.Unsupported(f"Cannot save {'compressed ' if compressed else ''}file with base extension '{ext}'."))
        return self

    # ====================================
    #      Reorganizing Existing Data     
    # ====================================

    # Given a (sub)set of column names in this data, reorder the data
    # making those names first (with unprovided names remaining in
    # same order as before).
    def reorder(self, names):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            raise(self.Empty("Cannot reorder empty data."))
        # Check for proper usage.
        my_names = set(self.names)
        for n in names:
            if (n not in my_names):
                raise(self.BadSpecifiedName(f"No column named '{n}' exists in this data."))
        # Construct the full reordering of every row of data.
        order = [self.names.index(n) for n in names]
        taken = set(order)
        order += [i for i in range(len(self.names)) if i not in taken]
        # Re-order names, types, and all rows of data.
        self.names = [self.names[i] for i in order]
        self.types = [self.types[i] for i in order]
        for row in range(len(self)):
            self[row] = [self[row,i] for i in order]


    # Given a new list of types, re-cast all elements of columns with
    # changed types into the newly specified type.
    def retype(self, types, columns=None):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            raise(self.Empty("Cannot retype empty data."))
        # Handle a single column or type being provided.
        if (type(columns) == str): columns = [columns]
        if (type(types) == type): types = [types] * len(columns)
        # Initialize "columns" if it was not provided to be the whole Data.
        if (type(columns) == type(None)):
            columns = self.names[:]
        if (len(types) != len(columns)):
            raise(self.BadSpecifiedType(f"{Data.retype} given {len(types)} types, epxected {len(columns)}."))
        for (new_t, col) in zip(types, columns):
            if (col not in self.names):
                raise(self.BadSpecifiedName(f"No column named '{col}' exists in this data."))
            c = self.names.index(col)
            old_t = self.types[c]
            if (new_t == old_t): continue
            # Update the stored type
            self.types[c] = new_t
            # Retype all non-missing elements in that column (in place)
            for i in range(len(self)):
                if (type(self[i,c]) == type(None)): continue
                try:
                    self[i,c] = new_t(self[i,c])
                except ValueError:
                    raise(self.BadSpecifiedType(f"Type casting {new_t} for column {c} is not compatible with existing value '{self[i,c]}' on row {i}."))

    # Given a new column, add it to this Data.
    def add_column(self, column, name=None, idx=None):
        # Set the default index to add a column to be the end.
        if (idx == None): idx = len(self.names)
        # Verify the name
        if (type(name) == type(None)):
            num = 0
            while (str(num) in self.names): num += 1
            name = str(num)
        elif (type(name) != str):
            raise(self.BadSpecifiedName(f"Only string names are allowed. Received name '{name}' with type {type(name)}."))
        # Set the name
        self.names.insert(idx, name)
        # Verify the column type dynamically. Add new values to all rows.
        new_type = type(None);
        for i,val in enumerate(column):
            # Verify valid index first..
            if (i >= len(self)):
                # Remove the added name.
                self.names.pop(idx)
                # Remove the added elements if the length was not right
                for j in range(len(self)): self[j].pop(idx)
                # Raise error for too long of a column
                raise(self.BadData(f"Provided column has at least {i+1} elements, more than the length of this data ({len(self)})."))
            # Append the value to this row..
            self[i].insert(idx, val)
            if (type(val) == type(None)): # Only add to missing values entry if it's not already there
                self.missing.add(id(self.data[i]))
            elif (new_type == type(None)):
                # Capture the new type (type must be right)
                new_type = type(val)
            elif (type(val) != new_type):
                # Remove the added name.
                self.names.pop(idx)
                # Remove the added elements because a problem was encountered.
                for j in range(i+1): self[j].pop(idx)
                # This is a new type, problem!
                raise(self.BadValue(f"Provided column has multiple types. Original type {new_type}, but '{val}' has type {type(val)}."))
        # Verify the column length
        if (i < len(self)-1):
            # Remove the added name.
            self.names.pop(idx)
            # Remove the added elements if the length was not right
            for j in range(i+1): self[j].pop(idx)
            # Raise error for too short of a column
            raise(self.BadData(f"Provided column has length {i+1}, less than the length of this data ({len(self)})."))
        # Finally, record the new type
        self.types.insert(idx, new_type)

    # Generate a pair of functions. The first function will map rows
    # of Data to a real-valued list. The second function will map
    # real-valued lists back to rows in the original space.
    def generate_mapping(self):
        if self.empty: raise(self.ImproperUsage("Cannot map empty Data."))
        # Make sure we can handle this type of data
        if not all(t in {int,float,str} for t in self.types):
            raise(self.Unsupported("A mapping is only designed for Data composed of <float>, <int>, and <str> typed values."))
        # Get those categoricals that need to be converted to reals
        names = self.names[:]
        types = self.types[:]
        cat_constants = {}
        cat_names = [n for (n,t) in zip(names,types) if (t == str)]
        cat_values = {n:sorted(set(self[n])-{None}) for n in cat_names}
        # Get the categorical simplex mappings.
        cat_mappings = {n:regular_simplex(len(cat_values[n]))
                        for n in cat_names}
        # Get the real values that have more than one unique value.
        real_values = {n:None for n in names if n not in cat_values}
        for n in real_values:
            unique = set()
            # Identify those columns (of reals) that have >1 unique value.
            for v in self[n]:
                unique.add(v)
                if len(unique) > 1: break
            else: real_values[n] = unique.pop()
        # Initialize stoarge for the indices
        num_inds = []
        cat_inds = []
        # Calculate the dimension of the real space (and the column names)
        real_names = []
        counter = 0
        for (n,t) in zip(names, types):
            if (t == str):
                real_names += [f"{n}-{i+1}" for i in range(len(cat_values[n])-1)]
                # Record the indices of categorical mappings
                cat_inds += list(counter+i for i in range(len(cat_values[n])-1))
                if (len(cat_inds) > 0): counter = cat_inds[-1] + 1
            elif (real_values[n] == None):
                real_names.append( n )
                # Record the indices of numerical values
                num_inds.append( counter )
                counter += 1

        # Generate a function that takes a row in the original space
        # and generates a row in the associated real-valued space.
        def to_real(row, fill=False):
            # Check for row length
            if (len(row) != len(names)):
                raise(self.BadElement(f"Provided row has {len(row)} elements, expected {len(names)}."))
            # Bulid the real-valued vector
            real_row = []
            for (n,t,v) in zip(names, types, row):
                # Check for the type of the elements of the row
                if (type(v) == type(None)):
                    if (t == str):
                        real_row += [None] * len(cat_mappings[n][0])
                    elif (real_values[n] == None):
                        real_row += [None]
                    continue
                elif (type(v) != t) and (t != type(None)):
                    try:    v = t(v)
                    except: raise(self.BadValue(f"Value '{v}' in provided row has wrong type. Got {type(v)}, expected {t} and could not succesfully cast."))
                # Handle the addition of a proper value
                if (t == str):
                    # WARNING: Not handling "unrecognized categories"
                    if v in cat_values[n]:
                        cat_index = cat_values[n].index(v)
                        mapped_vec = cat_mappings[n][cat_index]
                        real_row += list(mapped_vec)
                    elif (not fill):
                        # Only raise an error if it's wanted
                        raise(self.BadValue(f"Unrecognized categorical value {v} for column {n}."))
                    else:
                        # Put a filler value in for the unknown (middle category)
                        real_row += [0.] * (len(cat_values[n]) - 1)
                elif (real_values[n] == None):
                    real_row += [float(v)]
            # Check for row length
            if (len(real_row) != len(real_names)):
                raise(self.BadElement(f"Provided row has {len(row)} elements, translation to {len(real_names)} unsuccesful, only created {len(real_row)}."))
            # Return the real-valued vector
            return real_row

        # Generate a function that takes a row in the real-valued
        # space and generates a row in the original space.
        def real_to(real_row, bias={}):
            import numpy as np
            # Verify length of real-vector
            real_row = list(real_row)
            if (len(real_row) != len(real_names)):
                raise(self.BadElement(f"Provided real vector has {len(real_row)} elements, expected {len(real_names)}."))
            # Build the original row
            row = []
            for (n,t) in zip(names, types):
                if (t == str):
                    num_categories = len(cat_values[n])
                    vec = np.array(real_row[:num_categories-1])
                    # Get the weights for each of the categorical values.
                    cat_weights = np.array(
                        [(bias[n][c] if (n in bias and c in bias[n]) else 1.0)
                         for c in cat_values[n]])
                    # Handle reverse-categorical mapping
                    closest_cat = np.argmin(np.sum((cat_mappings[n] - vec)**2, axis=1) / cat_weights)
                    category = cat_values[n][closest_cat]
                    row += [category]
                    # Reduce the length of the real-valued row
                    real_row = real_row[num_categories-1:]
                elif (real_values[n] == None):
                    if (t == int):
                        # Use the rounded version of the integer
                        row += [int(real_row.pop(0) + .5)]
                    else:
                        # Use the float representation of the numpy floatf
                        row += [float(real_row.pop(0))]
                else:
                    # Fill with the only value that was provided.
                    row += [real_values[n]]
            if (len(real_row) > 0):
                raise("Left over values!")
            # Check for row length
            if (len(row) != len(names)):
                raise(self.BadElement(f"Provided row has {len(real_row)} elements, translation to {len(names)} unsuccesful, only created {len(row)}."))
            return row
        # Return the two mapping functions and some info
        return to_real, real_to, real_names, num_inds, cat_inds

    # Convert this Data automatically to a real-valued array.
    # Return real-valued array, function for going to real from
    # elements of this Data, function for going back to elements of
    # this data from a real vector.
    def to_matrix(self, print_column_width=70):
        import numpy as np
        to_real, real_to, real_names, num_inds, cat_inds = self.generate_mapping()
        data = []
        for row in self:
            # How to handle rows with missing values? Add column for "is missing"?
            if (id(row) in self.missing): continue
            data.append( to_real(row) )
        if (len(data) == 0):
            raise(self.BadData("No rows in this data are complete, resulting matrix empty."))
        array = np.array(data)
        # Container for important post-processing information. 
        # 
        # Numeric.to_real -- the function for processing regular
        #                    vectors into associated real vectors
        # Numeric.real_to -- the function for processing real vectors
        #                    back into regular vector
        # Numeric.names   -- meaningful column names for the real vectors
        # Numeric.nums    -- standard numerical indices
        # Numeric.cats    -- translated categorical indices
        # Numeric.data    -- The numeric version of the data.
        # Numeric.shift   -- How to shift this data to have min 0.
        # Numeric.scale   -- How to scale this data (after shift) to 
        #                    have a domain of [0., 1.]
        class Numeric:
            def __str__(self):
                opener = "  ["
                str_real_names = [opener]
                for n in real_names:
                    if not ((len(str_real_names[0]) == len(opener)) or
                            (len(str_real_names[-1]) + len(n) < print_column_width)):
                        str_real_names.append("   ")
                    str_real_names[-1] += f"'{n}', "
                str_real_names[-1] = str_real_names[-1][:-2] + "]"
                str_real_names = "\n".join(str_real_names)
                return f"Numeric object for data with {len(real_names)} columns named:\n{str_real_names}"
        # Generate a container object and fill it.
        numeric = Numeric()
        numeric.to_real = to_real
        numeric.real_to = real_to
        numeric.names = real_names
        numeric.nums = num_inds
        numeric.cats = cat_inds
        numeric.data = array
        # Compute the shift and scale for normalization (of numeric values)
        numeric.shift = np.min(array, axis=0)
        numeric.scale = np.max(array, axis=0) - numeric.shift
        numeric.shift[cat_inds] = 0.
        numeric.scale[cat_inds] = 1.
        # Return the container object
        return numeric

    # Generate a compact numpy structured array from the contents of self.
    def to_struct(self):
        import numpy as np
        # Generate the dtypes
        dtypes = []
        for i,(n,t) in enumerate(zip(self.names,self.types)):
            if (t != str):
                pair = (n,) + NP_TYPES[t]
            else:
                max_len = max(len(s) for s in self[n] if type(s) == str)
                pair = (n,) + NP_TYPES[t][:1]+(max_len,)
            dtypes.append(pair)
        dtypes = np.dtype(dtypes)
        # Return the numpy structured array
        return np.array([tuple(row) for row in self if (id(row) not in self.missing)], 
                        dtype=dtypes)

    # Collect the dictionaries of unique values (with counts) for each column.
    def counts(self, columns=None):
        if (type(columns) == type(None)): columns = self.names
        column_info = {n:{} for n in self.names}
        for row in self:
            for n,val in zip(columns, row):
                # If this column has been removed (because it is
                # unhashable), then skip it in the processing
                if n not in column_info: pass
                # Try to add the new value
                try:
                    column_info[n][val] = column_info[n].get(val,0) + 1
                except TypeError:
                    # Cannot hash an element of self[n], disable tracking
                    column_info.pop(n)
        return column_info

    # Generate a new data with only the unique rows (by content) from this data.
    def unique(self, display_wait_sec=DEFAULT_DISPLAY_WAIT):
        import time
        from util.system import hash
        unique = Data(names=self.names, types=self.types)
        found = set()
        start = time.time()
        for i,row in enumerate(self):
            # Update user on progress if too much time has elapsed..
            if (time.time() - start) > display_wait_sec:
                print(f" {100.*i/len(self):.2f}%", end="\r", flush=True)
                start = time.time()
            # Get the hashed value by consecutively hashing all of the
            # values within this row. Hashing the entire row sometimes
            # produces changing hashes (likely based on memory locations).
            hashed = ''
            for v in row: hashed = hash( hashed + hash(v) )
            if hashed not in found:
                found.add(hashed)
                unique.append(row)
        return unique

    # Given a Data that has at least one column with the same names as
    # a column in self, collect all those values in non-mutual names
    # in lists of values where the shared columns had the same content.
    # 
    # Example:
    #   Combined with "unique", this function can be used to convert a
    #   struct where some rows are repeats of each other having unique
    #   values in only one column to a Data object with only unique
    #   rows and one column contains lists of values associated with
    #   each unique row.
    # 
    #   densified_data = data[some_columns].unique().collect(data)
    # 
    def collect(self, data, display_wait_sec=DEFAULT_DISPLAY_WAIT):
        import time
        from util.system import hash
        # Get the names of shared columns and indices in each.
        match_names = set(n for n in self.names if n in data.names)
        self_columns = [i for (i,n) in enumerate(self.names) if n in match_names]
        # Sort the columns of self to be in the same order as that of other.
        self_columns = sorted(self_columns, key=lambda i: data.names.index(self.names[i]))
        other_columns = [i for (i,n) in enumerate(data.names) if n in match_names]
        # Collection columns
        names_to_collect = [n for n in data.names if n not in match_names]
        indices_to_collect  = [i for (i,n) in enumerate(data.names) if n not in match_names]
        collections = {n:[] for n in names_to_collect}
        # Add columns for holding the matches
        for n in names_to_collect:
            self.add_column(([] for i in range(len(self))), name=n)
        # Generate a set of lookups (hashed values even for mutable types)
        index_of_row = {}
        for i,row in enumerate(self):
            hash_value = hash([row[c] for c in self_columns])
            index_of_row[hash_value] = i
        # Paired list of source -> destination indices for copying.
        source_dest = list(zip(indices_to_collect, names_to_collect))
        start = time.time()
        # Perform the serach process
        for i,row in enumerate(data):
            # Update user on progress if too much time has elapsed..
            if (time.time() - start) > display_wait_sec:
                print(f" {100.*i/len(data):.2f}%", end="\r", flush=True)
                start = time.time()
            # Convert row into its hash value, lookup storage location.
            hash_value = hash([row[c] for c in other_columns])
            if hash_value in index_of_row:
                row_idx = index_of_row[hash_value]
                for source, dest in source_dest:
                    self[row_idx, dest].append(row[source])
            else:
                print("Could not find:", [row[c] for c in other_columns])
        # Return the self (that has been modified)
        return self

    # Return an iterator that provides "k" (rest, fold) sub-data
    # from this data. Randomly shuffle indices with seed "seed".
    def k_fold(self, k=10, seed=0, only_indices=False):
        # Generate random indices
        import random
        old_state = random.getstate()
        random.seed(seed)
        indices = list(range(len(self)))
        random.shuffle(indices)
        # Reset the random number generator as not to muck with
        # other processes that might be using it.
        random.setstate(old_state)
        # Store some variables for generating the folds
        total = len(self)
        if (type(self) != Data): only_indices = True
        if not only_indices:     source = self
        for batch in range(k):
            # Find the indices for training and testing
            first = int(.5 + batch * total / k)
            last = int(.5 + (batch + 1) * total / k)
            fold = indices[first : last]
            rest = indices[:first] + indices[last:]
            if only_indices: yield rest, fold
            else:            yield source[rest], source[fold]

    # Use a prediction model to fill missing values in a provided row.
    def fill(self, row, model=None, weights=None, bias={}):
        import numpy as np
        numeric = self.to_matrix()
        # Pick the filling model based on data size.
        if (type(model) == type(None)):
            # Use Voronoi if there are fewer than 10000 points.
            if (len(self) <= 5000):
                from util.approximate import Voronoi, unique
                model = unique(Voronoi)()
            # Otherwise use nearest neighbor (can't get points + weights from NN)
            else:
                from util.approximate import unique, condition
                samples = min(numeric.data.shape[0], 10000)
                dim = min(numeric.data.shape[1], 1000)
                # Use nearest neighbor (can't get points + weights from NN)
                from util.approximate import NearestNeighbor
                model = condition(unique(NearestNeighbor), dim=dim, samples=samples)()
        # Set the weights accordingly
        if (type(weights) == type(None)):
            weights = np.ones(numeric.data.shape[1])
        elif (len(weights) != len(self.names)):
            raise(self.BadValue("Weights vector is length {len(weights)}, expected length {len(self.names)}."))
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
        row_known = [i for i in range(len(row)) if (type(row[i]) != type(None))]
        row_unknown = [i for i in range(len(row)) if (type(row[i]) == type(None))]
        if len(row_unknown) == 0: raise(self.BadValue("Provided row had no missing values."))
        # Get the numeric version of the provided row.
        numeric_row = numeric.to_real(row, fill=True)
        # Identify the indices of known values and unknown values.
        known_values = np.array([i for i in range(len(numeric_row))
                                 if (type(numeric_row[i]) != type(None))])
        unknown_values = np.array([i for i in range(len(numeric_row))
                                   if (type(numeric_row[i]) == type(None))])
        # Get the part of the numeric row
        no_none_numeric_row = np.array([numeric_row[i] for i in known_values])
        # Normalize the numeric row (to debias numeric scale differences)
        no_none_numeric_row = ((no_none_numeric_row - numeric.shift[known_values])
                               / numeric.scale[known_values]) * weights[known_values]
        # Get the normalized training point locations.
        normalized_numeric_data = ((numeric.data - numeric.shift)
                                   / numeric.scale) * weights
        # Fit the model and make a prediction for this row.
        if hasattr(model, "WeightedApproximator"):
            # Fit the model (just to the points)
            model.fit(normalized_numeric_data[:,known_values])
            del(normalized_numeric_data)
            # Get the source points and source weights for the prediction
            pts, wts = model.predict(np.array(no_none_numeric_row))
            # Compute the fill row by performing the weighted sum.
            fill_row = np.sum((numeric.data[pts,:][:,unknown_values].T * wts), axis=1)
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
        # python3 -m pytest data.py
        # Transfer the filled numeric row into the original row
        fill_row = numeric.real_to(fill_row, bias=bias)
        for i,val in enumerate(fill_row):
            if type(row[i]) == type(None):
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

    # Given a set of columns of this data to try and predict, run a
    # k-fold cross validation over predictions using 'fill'.
    def predict(self, target_columns, model=None, weights=None, bias={}, k=10):
        import numpy as np
        # Get the numeric representation of self.
        numeric = self.to_matrix()
        # Pick the filling model based on data size.
        if (type(model) == type(None)):
            from util.approximate import unique, condition
            samples = min(numeric.data.shape[0], 10000)
            dim = min(numeric.data.shape[1], 1000)
            # Use nearest neighbor (can't get points + weights from NN)
            from util.approximate import NearestNeighbor
            model = condition(unique(NearestNeighbor), dim=dim, samples=samples)()
        # Set the weights accordingly
        if (type(weights) == type(None)):
            weights = np.ones(numeric.data.shape[1])
        elif (len(weights) != len(self.names)):
            raise(self.BadValue("Weights vector is length {len(weights)}, expected length {len(self.names)}."))
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
            raise(self.BadTarget("Target must either be a string (column name) or list of strings (column names)."))
        # Get the data all normalized (and weighted)
        normalized_data = weights * (numeric.data - numeric.shift) / numeric.scale
        print("normalized_data.shape: ",normalized_data.shape)
        # Identify those columns that are being predicted
        target_column_indices = [self.names.index(c) for c in target_columns]
        results = Data(names=[c + " Predicted" for c in target_columns]+["Prediction Index"],
                       types=[self.types[i] for i in target_column_indices]+[int])
        # Identify the numeric ranges of those columns
        indices = set(target_column_indices)
        sample = numeric.to_real([
            (v if i in indices else None)
            for i,v in enumerate(self[0])])
        del(indices)
        train_cols = np.array([i for (i,v) in enumerate(sample) if type(v) == type(None)])
        test_cols =  np.array([i for (i,v) in enumerate(sample) if type(v) != type(None)])
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
                row = numeric.real_to(row, bias=bias)
                results.append( [row[i] for i in target_column_indices] + [i] )
            print(" "*70,end="\r",flush=True)
        # Sort results by the original index
        results.sort(key=lambda r: r[-1])
        results.pop("Prediction Index")
        # Return the results in a new Data object.
        return results

    # Compute the pairwise effects between columns. Use the following:
    # 
    #    number vs. number     -- Correlation coefficient between the two sequences.
    #    category vs. number   -- "method" 1-norm difference between full distribution
    #                             and conditional distributions for categories.
    #    category vs. category -- "method" total difference between full distribution
    #                             of other sequence given value of one sequence.
    # 
    # Print out the sorted table of pairwise effects between columns,
    # showing the highest effects first, the smallest last.
    def effect(self, compare_with=None, method="mean", display=True, **kwargs):
        from itertools import combinations
        from util.stats import effect
        if (type(compare_with) == type(None)): compare_with = set(self.names)
        elif (type(compare_with) == list):     compare_with = set(compare_with)
        else:                                  compare_with = {compare_with}
        # Print out the effect in column form (not matrix form).
        effs = []
        for (col_1, col_2) in combinations(self.names, 2):
            if (col_1 in compare_with):
                eff = effect(list(self[col_1]), list(self[col_2]), method=method, **kwargs)
                effs.append( (col_1, col_2, eff) )
            elif (col_2 in compare_with):
                eff = effect(list(self[col_1]), list(self[col_2]), method=method, **kwargs)
                effs.append( (col_2, col_1, eff) )
        # Convert an effect to a sortable number (lowest most important).
        def to_num(eff): 
            if type(eff) == dict: return -sum(eff.values()) / len(eff)
            else:                 return -abs(eff)
        # Convert an effect to a printable string.
        def to_str(eff):
            if type(eff) == dict:
                eff_str = []
                for key in sorted(eff, key=lambda k: -abs(eff[k])):
                    eff_str += [f"'{key}':{eff[key]:5.2f}"]
                return f"{-to_num(eff):5.2f}  {{"+", ".join(eff_str)+"}"
            else:
                return f"{eff:5.2f}"
        # Sort the data by magnitude of effect. (most effected -> first)
        effs = sorted(effs, key=lambda row: to_num(row[-1]))
        # Do a quick return without printing if display is turned off.
        if (not display): return effs
        # Print out the effects between each column nicely.
        max_len_1 = max( max(map(lambda i:len(i[0]), effs)), len("Column 1") )
        max_len_2 = max( max(map(lambda i:len(i[1]), effs)), len("Column 2") )
        header = f"{'Column 1':{max_len_1}s}  |  {'Column 2':{max_len_2}s}  |  Effect"
        rows = []
        for (col_1, col_2, eff) in effs:
            eff = to_str(eff)
            rows += [f"{col_1:{max_len_1}s}  |  {col_2:{max_len_2}s}  |  {eff}"]
        row_len = max(len(header), max(map(len, rows)))
        rows = ["", '-'*row_len, header, '-'*row_len] + rows + ['-'*row_len, ""]
        for row in rows: print(row)
        # Return the sorted effect triple list [(col 1, col 2, effect)].
        return effs

    # Give an overview (more detailed than just "str") of the contents
    # of this data. Useful for quickly understanding how data looks.
    def summarize(self, max_display=None):
        # Special case for an empty data
        if len(self) == 0:
            print(self)
            return
        # Set the "max_display" to the default value for this class
        if type(max_display) == type(None): max_display = self.max_display
        print("SUMMARY:")
        print()
        print(f"  This data has {len(self)} row{'s' if len(self) != 1 else ''}, {len(self[0])} column{'s' if len(self[0]) != 1 else ''}.")
        num_ordinal = sum(1 for t in self.types if t in {float,int})
        print(f"    {num_ordinal} column{'s are' if num_ordinal != 1 else ' is'} recognized as ordinal, {len(self[0])-num_ordinal} {'are' if (len(self[0])-num_ordinal) != 1 else 'is'} categorical.")
        print(f"    {len(self.missing)} row{'s have' if len(self.missing) != 1 else ' has'} missing values.")
        print()
        print("COLUMNS:")
        print()
        name_len = max(map(len, self.names))
        type_len = max(map(lambda t: len(str(t)), self.types))
        # Describe each column of the data
        for c,(n,t) in enumerate(zip(self.names, self.types)):
            # Count the number of elements for each value
            counts = {}
            to_string = False
            for val in self[n]:
                if to_string: val = str(val)
                try: counts[val] = counts.get(val,0) + 1
                except TypeError:
                    to_string = True
                    val = str(val)
                    counts[val] = counts.get(val,0) + 1
            print(f"  {c:{len(str(len(self.names)))}d} -- \"{n}\"{'':{1+name_len-len(n)}s}{str(t):{type_len}s} ({len(counts)} unique value{'s' if (len(counts) != 1) else ''})")
            # Remove the "None" count from "counts" to prevent sorting problems
            none_count = counts.pop(None, 0)
            # For the special case of ordered values, reduce to ranges
            if (t in {int,float}) and (len(counts) > max_display):
                if (none_count > 0):
                    perc = 100. * (none_count / len(self))
                    print(f"    None                 {none_count:{len(str(len(self)))}d} ({perc:5.1f}%) {'#'*round(perc/2)}")
                # Order the values by intervals and print
                min_val = min(counts)
                max_val = max(counts)
                width = (max_val - min_val) / (max_display-1)
                for i in range(max_display-1):
                    lower = min_val + width*i
                    upper = min_val + width*(i+1)
                    if (i == (max_display - 2)):
                        num = sum(counts[v] for v in counts if lower <= v <= upper)
                        cap = "]"
                    else:
                        num = sum(counts[v] for v in counts if lower <= v < upper)
                        cap = ")"
                    perc = 100. * (num / len(self))
                    print(f"    [{lower:9.2e}, {upper:9.2e}{cap} {num:{len(str(len(self)))}d} ({perc:5.1f}%) {'#'*round(perc/2)}")
            else:
                if t in {int, float}:
                    # Order the values by their inate ordering
                    ordered_vals = sorted(counts)
                else:
                    # Order the values by their frequency and print
                    ordered_vals = sorted(counts, key=lambda v: -counts[v])
                val_len = max(map(lambda v: len(str(v)), counts))
                val_len = max(val_len, len("None"))
                val_len = min(val_len, self.max_str_len)
                if (t == str): val_len += 2
                if (none_count > 0):
                    perc = 100. * (none_count / len(self))
                    print(f"    {'None':{val_len}s}{'  ' if (t == str) else ''} {none_count:{len(str(len(self)))}d} ({perc:5.1f}%) {'#'*round(perc/2)}")
                for val in ordered_vals[:max_display]:
                    perc = 100. * (counts[val] / len(self))
                    if (t == str):
                        print(f"    \"{val}\"{'':{1+val_len-len(str(val))}s}{counts[val]:{len(str(len(self)))}d} ({perc:5.1f}%) {'#'*round(perc/2)}")
                    else:
                        print(f"    {self._val_to_str(val):{val_len}s} {counts[val]:{len(str(len(self)))}d} ({perc:5.1f}%) {'#'*round(perc/2)}")
                if (len(ordered_vals) > max_display):
                    print("    ... (increase 'max_display' to see more summary statistics).")
            print()




