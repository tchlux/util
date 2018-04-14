from util import COMMON_SEPERATORS, UPDATE_FREQUENCY, \
    MAX_ERROR_PRINTOUT, NP_TYPES, PY_TYPES, GENERATOR_TYPE
from util.system import AtomicOpen

# Given "d+1" categories that need to be converted into real
# space, generate a regular simplex in d-dimensional space. (all
# points are equally spaced from each other and the origin) This
# process is more efficient than "one-hot" encoding by one
# dimension. It also guarantees that all points are not placed
# explicitly on a sub-dimensional manifold.
def regular_simplex(num_categories):
    import numpy as np
    d = num_categories
    # Initialize all points to be zeros.
    points = np.zeros((d,d-1))
    # Set the initial first point as 1
    points[0,0] = 1
    # Calculate all the intermediate points
    for i in range(1,d-1):
        # Set all points to be flipped from previous calculation while
        # maintaining the angle "arcos(-1/d)" between vectors.
        points[i:,i-1] = -1/(d-i) * points[i-1,i-1]
        # Compute the new coordinate using pythagorean theorem
        points[i,i] = (1 - sum(points[i,:i]**2))**(1/2)
    # Set the last coordinate of the last point as the negation of the previous
    points[i+1,i] = -points[i,i]
    # Return the regular simplex
    return points

# Given a point, determine the affine combination of categories in the
# defining regular simplex that produces the point.
def category_ratio(point):
    import numpy as np
    categories = regular_simplex(len(point)+1)
    # Add the "constant" column to the matrix of categories
    categories = np.hstack((categories, np.ones((len(point)+1,1))))
    # Add the "sum to 1" column to the affinely created point
    point = np.hstack((point, 1))
    # Calculate and return the affine weights
    return np.linalg.solve(categories.T, point)


# =================================================
#      Automatically Reading Structured Arrays     
# =================================================

# Return the type of a string (try converting to int and float)
def get_type(string):
    if (len(string) == 0):
        return int
    try:
        int(string)
        return int
    except ValueError:
        try:
            float(string)
            return float
        except ValueError:
            return str

# Pre:  "filename" is the name of an existing file (including path)
# Post: "filename" is scanned for viable seperators, those are all
#       characters used the same number of times on each line. If at
#       the end of the file there are still multiple candidate
#       seperators, they are sorted by most-occuring first. If one of
#       COMMON_SEPERATORS exist in the set, it is used, otherwise the
#       most occurring seperator is used.
def detect_seperator(filename="<no_provided_file>", verbose=False, opened_file=None):
    lines = []
    not_viable = set()
    if not opened_file:
        f = AtomicOpen(filename)
    else:
        f = opened_file
    # Identify the potential seperators from the first line
    first_line = f.readline().strip()
    seperators = {char:first_line.count(char) for char in set(first_line)}
    # Cycle the file (narrowing down the list of potential seperators
    for line in f.readlines():
        line_chars = set(line.strip())
        # Make sure all tracked seperators appear the correct
        # number of times in the line
        for char in list(seperators.keys()):
            if (line.count(char) != seperators[char]):
                not_viable.add( char )
                seperators.pop( char )
        # WARNING: Do *not* break early, to verify that the
        #          seperators are used correctly throughout the
        #          entire file!
    if not opened_file:
        # Close the opened file object
        f.close()
    else:
        # Go back to the beginning of the file
        f.seek(0)
    # Get the seperators from the dictionary sorted by the decreasing
    # number of occurrences the character has on each line of the data.
    seperators = sorted(seperators.keys(), key=lambda k: -seperators[k])
    # Check to make sure there is a seperator detected
    if len(seperators) == 0:
        raise(Exception("ERROR: No consistently used seperator in file."))
    if verbose:
        print(" Detected the following possible seperators:\n  ", seperators)
    # Search through the common seperators first, using one of them if possible
    for sep in COMMON_SEPERATORS:
        if sep in seperators:
            if verbose: print(" Using '%s' as seperator."%sep)
            break
    else:
        # Otherwise, use the most occurring character that could be a seperator
        sep = seperators[0]
        print("WARNING: Using '%s' as a seperator, even though it is unexpected."%(sep))
    # Return the automatically detected seperator
    return sep

# Pre:  "filename" is the name of a txt file that has a header
#       "sep" is the seperator used in the data file
#       "display" is True if the user wants the detected types printed
# Post: A data.Struct class containing all data from the file.
def read_struct(filename="<no_provided_file>", sep=None, types=None,
                verbose=False, opened_file=None):
    import time
    if sep == None: sep = detect_seperator(filename, verbose, opened_file)
    # Get the opened file object
    if (not opened_file):
        # Open the file to get the data
        f = AtomicOpen(filename)
    else:
        f = opened_file
        file_name = f.name
    # Get the first line, assuming it's the header
    line = f.readline()
    header = [n.strip() for n in line.strip().split(sep)]
    # Get the first line of data and learn the types of each column
    line = f.readline().strip().split(sep)
    # Automatically detect types (if necessary)
    type_digression = False
    if type(types) == type(None):
        types = list(map(get_type, line))
        type_digression = True
    # Initialize our data holder
    data = Struct(names=header, types=types)
    # Add the first line to the data store (after inserting "None")
    while "" in line: line[line.index("")] = None
    data.append(line)
    # Get the rest of the raw data from the file and close it
    raw_data = f.readlines()
    if (not opened_file):
        # Close the file
        f.close()
    # WARNING: Provided file object left to be closed by caller
    # Print out a nicely formatted description of the detected types
    if verbose:
        print("Column names and types:")
        first_line  = ""
        second_line = ""
        for h,t in zip(header,types):
            length = max(len(h), len(t.__name__))
            first_line += f"  {h:{length}s}"
            second_line += f"  {t.__name__:{length}s}"
        print(first_line)
        print(second_line)
    # Now process the rest of the data, dynamically overwriting old data
    errors = []
    last_update = 0
    should_update = lambda: (verbose and ((time.time() - last_update) > UPDATE_FREQUENCY))
    for i,line in enumerate(raw_data):
        # Update the user on progress if appropriate
        if should_update():
            last_update = time.time()
            print("\r%0.1f%% complete"% (100.0*i/len(raw_data)), flush=True, end="")
        # Read the line of data into a list format (replace empty string with None)
        list_line = line.strip().split(sep)
        # Replace any empty strings with "None" for missing values
        while "" in list_line: list_line[list_line.index("")] = None
        # Try and add the line of data
        try:
            data.append(list_line)
        except Struct.BadValue:
            if type_digression:
                # Update the types based on digression order
                new_types = list(map(get_type, line.strip().split(sep)))
                type_order = {int:0, float:1, str:2}
                new_types = [max(old,new, key=lambda v: type_order[v])
                             for (old,new) in zip(types, new_types)]
                # Update the user if appropriate
                if verbose:
                    print(f"\nType digression because of line {i+3}."+
                          f"\n  old types: {types}"+
                          f"\n  new types: {new_types}")
                # Retype the existing data to match the new types
                data.retype(new_types)
                # Append line that matches new types
                data.append(list_line)
                types = new_types
            else:
                errors.append(i+3)
                if len(errors) <= MAX_ERROR_PRINTOUT:
                    print(f"\nERROR ON LINE {i+3}:\n  {line}")
                if len(errors) == MAX_ERROR_PRINTOUT:
                    print("Suppressing further error output. See all erroneous lines in later message.")
    if verbose and (last_update != 0): print("\r100.0% complete.")
    if (len(errors) > 0) and (not verbose):
        print( "WARNING: Encountered potentially erroneous lines while\n"+
              f"         reading '{filename}'.\n"+
               "         Consider setting 'verbose=True' and checking errors.")
    if verbose and (len(errors) > MAX_ERROR_PRINTOUT):
        print(f"ERRONEOUS LINES FROM DATA:\n  {errors}")
    # Return the struct
    return data

# ======================================================
#      Python 'Struct' with named and typed columns     
# ======================================================

# Define a python structure (named columns)
class Struct(list):
    # Default values for 'self.names' and 'self.types'
    names = None
    types = None
    # Default maximum printed rows of this struct
    _max_display = 10

    # Local descriptive exceptions
    class NameTypeMismatch(Exception): pass
    class BadSpecifiedType(Exception): pass
    class BadSpecifiedName(Exception): pass
    class NoNamesSpecified(Exception): pass
    class UnknownName(Exception): pass
    class Unsupported(Exception): pass
    class BadAssignment(Exception): pass
    class BadElement(Exception): pass
    class BadValue(Exception): pass
    class BadIndex(Exception): pass
    class BadData(Exception): pass
    # Local mutable Column class (for column item assignment,
    # retrieval, and iteration over a column in a struct).
    class Column:
        index = 0
        def __init__(self, struct, column):
            # Store the parent struct and the column number 
            self.struct = struct
            self.column = column
        def __setitem__(self, index, value):
            # Default setitem for the parent struct
            self.struct[index, self.column] = value
        def __getitem__(self, index):
            # Handle integer indices
            if (type(index) == int):
                return self.struct[index, self.column]
            # Handle slices as indices
            return [self.struct[i,self.column] for
                    i in range(len(self.struct))[index]]
        def __len__(self): return len(self.struct)
        # Define this as an interator
        def __iter__(self):
            self.index = 0
            return self
        def __next__(self):
            if (self.index >= len(self.struct)):
                raise(StopIteration)
            self.index += 1
            return self.struct[self.index-1, self.column]
        # Use the iteration to generate a string of this column
        def __str__(self):
            return str(list(self))
        # Iterate over the indices that are equivalent to a value
        def __eq__(self, val):
            for i,v in enumerate(self):
                if (v == val):
                    yield i

    #      Overwritten list Methods     
    # ==================================

    def __init__(self, data=None, names=None, types=None):
        super().__init__()
        # Verify the length of the types and names
        if ((type(types) != type(None)) and
            (type(names) != type(None)) and
            (len(types) != len(names))):
            raise(self.NameTypeMismatch(f"Length of provided names {len(names)} does not match length of provided types {len(types)}."))
        # Get the names (must be strings)
        if (type(names) != type(None)):
            if any(type(n) != str for n in names):
                raise(self.BadSpecifiedName(f"An entry of provided names was not of {str} class."))
            self.names = names[:]
        # Store the types if provided
        if (type(types) != type(None)):
            if any(type(t) != type(type) for t in types):
                raise(self.BadSpecifiedType(f"An entry of provided types was not of {type(type)} class."))
            self.types = types[:]
        # WARNING: self.missing never shrinks, only grows with life of struct.
        self.missing = set()
        # If data was provided, put it into this struct
        if type(data) != type(None):
            for row in data:
                self.append(row)

    # Redefined append operator that 
    def append(self, element):
        missing_values = False
        # Check length in case names already exist
        if (type(self.names) != type(None)):
            if (len(element) != len(self.names)):
                raise(self.BadElement(f"Only elements of length {len(self.names)} can be added to this struct."))
        # Check length in case types already exist
        if (type(self.types) != type(None)):
            if len(element) != len(self.types):
                raise(self.BadElement(f"Only elements of length {len(self.types)} can be added to this struct."))                
            # Try type casing the new element
            for i, (val, typ) in enumerate(zip(element, self.types)):
                # Record 'None' types as missing entries
                if type(val) == type(None):
                    missing_values = True
                elif (type(val) != typ):
                    # If not missing, then check the type
                    try:
                        # Try casting the element as the expected type
                        element[i] = typ(val)
                    except ValueError:
                        # Otherwise raise an error to the user
                        raise(self.BadValue(f"Value '{val}' of type '{type(val)}' could not successfully be cast as '{typ}'."))
        else:
            # Construct the expected "types" of this struct
            self.types = [type(val) for val in element]
        if (type(self.names) == type(None)):
            # Construct default names for each column
            self.names = [str(i) for i in range(len(element))]
        # Call the standard append operator, adding the element to self
        super().append(element)
        # Add to the list of missing values if necessary
        if (missing_values): self.missing.add(id(self[-1]))


    # Make sure the copy of this struct is a deep copy
    def copy(self): return self[:]

    # Overwrite the standard "[]" get operator to accept strings as well.
    def __getitem__(self, index):
        if (type(index) == str):
            # Get the column of information under the given name
            if type(self.names) == type(None):
                raise(self.NoNamesSpecified("This struct does not have assigned 'names'."))
            if (index not in self.names):
                raise(self.UnknownName(f"This struct does not have a column named '{index}'."))
            col_number = self.names.index(index)
            # Return a mutable "Column" object
            return type(self).Column(self, col_number)
        elif (type(index) == slice):
            # Construct a new "struct" with just the specified rows (deep copy)
            new_struct = type(self)(names=self.names, types=self.types)
            for row in super().__getitem__(index):
                new_struct.append(row[:])
            return new_struct                
        elif (type(index) in {tuple,list}):
            if all(type(i) == int for i in index):
                if (type(index) == tuple) and (len(index) == 2):
                    # Return the multi-dimensional index access
                    return self[index[0]][index[1]]
                else:
                    # Assume that the user is accessing many rows by index
                    new_struct = Struct(names=self.names, types=self.types)
                    for row in index:
                        new_struct.append( self[row][:] )
                    return new_struct
            elif all(type(i) in {int,slice} for i in index):
                if (len(index) > 2):
                    raise(self.BadIndex(f"The provided index, {index}, is not understood."))
                # Make sure that both indices are slices
                if type(index[0]) == int:
                    index = (slice(index[0],index[0]+1), index[1])
                if type(index[1]) == int:
                    index = (index[0], slice(index[1],index[1]+1))
                # Iterate over both slices
                new_struct = type(self)(names=self.names[index[1].start : index[1].stop : index[1].step],
                                        types=self.types[index[1].start : index[1].stop : index[1].step])
                for row in self[index[0].start : index[0].stop : index[0].step]:
                    new_struct.append(
                        row[index[1].start : index[1].stop : index[1].step]
                    )
                return new_struct
            elif all(type(t)==str for t in index):
                # The user is slicing by column name, get the requested
                # columns and put them into a new struct
                types = [self.types[self.names.index(n)] for n in index]
                new_struct = type(self)(names=list(index), types=types)
                for row in zip(*(self[n] for n in index)):
                    new_struct.append(list(row))
                return new_struct
            else:
                # This is not a recognized index.
                raise(self.BadIndex(f"The provided index, {index}, is not understood."))
        elif (type(index) == GENERATOR_TYPE):
            # Iterate over integer indices provided by generator
            new_struct = type(self)(names=self.names, types=self.types)
            for row in index:
                if (type(row) != int):
                    raise(self.BadIndex(f"Generators as indices are only allowed to produce integers, encountered {type(row)}."))
                new_struct.append(self[row][:])
            return new_struct
        elif ((type(index) == int) and 
              (index >= len(self)) and
              (index < len(self)*len(self[0]))):
            # This index is accessing this structure serially
            row = index // len(self[0])
            col = index % len(self[0])
            return self[row][col]
        else:
            # If this is a normal index, use the standard get-item
            return super().__getitem__(index)

    # Overwrite standard "[]" assignment operator to accept strings.
    def __setitem__(self, index, value):
        if (type(index) == str):
            # Get the column of information under the given name
            if type(self.names) == type(None):
                raise(self.NoNamesSpecified("This struct does not have assigned 'names'."))
            if (index not in self.names):
                raise(self.UnknownName(f"This struct does not have a column named '{index}'."))
            col_number = self.names.index(index)
            # The new type will be recovered from the assigned values
            new_type = None
            # Set the values in self
            for i,v in enumerate(value):
                if (i >= len(self)):
                    # If there are too many values, raise an error
                    raise(self.BadAssignment(f"Column assignment requires a column of length {len(self)}, received object of at least length {i+1}."))
                # Update the value of the next entry in the struct
                self[i][col_number] = v
                if (type(v) == type(None)):
                    # If this is a missing value, update internal records
                    if (id(self[i]) not in self.missing):
                        self.missing.add(id(self[i]))
                elif (type(new_type) == type(None)):
                    new_type = type(v)
                elif (type(v) != new_type):
                    # If this type is unexpected, raise an error
                    raise(self.BadAssignment(f"Assigned column values were of type {new_type}, but encountered second type {type(v)}."))
            else:
                # If there were not enough values, raise an error
                if (i != (len(self)-1)):
                    raise(self.BadAssignment(f"Column assignment requires a column of length {len(self)}, received object of only length {i+1}."))
            if (type(new_type) != type(None)):
                # Update the new type stored in this column of the struct
                self.types[col_number] = new_type
        elif (type(index) == tuple):
            # Otherwise, if the user provides a tuple of integers + slices
            if all(type(i) in {int,slice} for i in index):
                if (len(index) > 2):
                    raise(self.BadIndex(f"The provided index, {index}, is not understood."))
                # Make sure that both indices are slices
                if type(index[0]) == int:
                    index = (slice(index[0],index[0]+1), index[1])
                if type(index[1]) == int:
                    index = (index[0], slice(index[1],index[1]+1))
                # Make sure the value has a length
                singleton = True
                try: 
                    len(value)
                    singleton = False
                except TypeError: pass
                # Make sure values is either iterable or 
                # Iterate over both slices
                count = 0
                for row in range(len(self))[index[0]]:
                    for col in range(len(self[0]))[index[1]]:
                        if singleton:
                            new_val = value
                        else:
                            new_val = value[count]
                            count += 1
                        if ((type(new_val) == type(None)) and
                            (id(self[row]) not in self.missing)):
                            self.missing.add(id(self[row]))
                        elif (type(new_val) != self.types[col]):
                            raise(self.BadAssignment(f"Assigned value '{new_val}' of type {type(new_val)} does not match required type {self.types[col]}."))
                        self[row][col] = new_val
            else:
                raise(self.BadIndex(f"The provided index, {index}, is not understood."))
        else:
            # If this is a normal index, use the standard get-item
            return super().__setitem__(index, value)

    # Generate a descriptive string of this struct.
    # WARNING: This operation costs O( len(self.missing) * len(self) )
    def __str__(self):
        string = "\n"
        string += f"{type(self)}\n"
        string += f" names ({len(self.names)}):\n  {self.names}\n\n"
        string += f" types:\n  {self.types}\n\n"
        string += f" values ({len(self)}):\n"
        for i,row in enumerate(self):
            string += f"  {row}\n"
            if i >= self._max_display: break
        if len(self) > self._max_display:
            string += "   ...\n"
        if len(self.missing) > 0:
            # Get the indices of missing elements
            indices = []
            for i in range(len(self)):
                if (id(self[i]) in self.missing):
                    indices.append(i)
                if len(indices) > self._max_display: break
            # Print out the missing values
            string += f"\n missing ({len(self.missing)}): ["
            # Print out the indices of the rows with missing values
            for i in indices[:self._max_display]:
                if (string[-1] != "["): string += ", "
                string += f"{i}"
            if (len(indices) > self._max_display):
                string += ", ..."
            string += "]\n"
            # Print out the actual rows with missing values
            for i in indices[:self._max_display]:
                string += f"  {self[i]}\n"
            if len(indices) > self._max_display:
                string += "   ..."
        return string + "\n"

    # Define a convenience funciton for concatenating another
    # similarly typed and named struct to self.
    def __iadd__(self, struct):
        # Copy names *only* if this struct does not have any
        if type(self.names) == type(None):
            self.names = struct.names[:]
        # Load in the struct
        for row in struct:
            self.append(row)
        # Return self
        return self

    # Convenience method for loading from file.
    def load(self, path):
        import pickle, dill, gzip
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
                self += pickle.load(f)
        elif (ext == "dill"):
            with file_opener(path, "rb") as f:
                self += dill.load(f)
        elif (ext in {"csv", "txt", "tsv"}):
            mode = "r" + ("t" if compressed else "")
            with file_opener(path, mode) as f:
                self += read_struct(opened_file=f)
        else:
            raise(self.Unsupported(f"Cannot load file with extension '{ext}'."))
        return self

    # Convenience method for saving to a file
    def save(self, path):
        import pickle, dill, gzip
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
            with file_opener(path, "wb") as f:
                pickle.dump(self, f)
        elif (ext == "dill"):
            with file_opener(path, "wb") as f:
                dill.dump(self, f)
        elif (ext == "csv"):
            mode = "w" + ("t" if compressed else "")
            list_to_str = lambda l: ",".join(l) + "\n"
            with file_opener(path, mode) as f:
                print(",".join(self.names), file=f)
                for row in self:
                    print_row = [str(v) if type(v) != type(None) else "" for v in row]
                    print(",".join(print_row), file=f)
        else:
            raise(self.Unsupported(f"Cannot save {'compressed ' if compressed else ''}file with base extension '{ext}'."))
        return self


    #      Custom Methods     
    # ========================

    # Given a new list of types, re-cast all elements of columns with
    # changed types into the newly specified type.
    def retype(self, types):
        if type(self.types) == type(None): self.types = types
        if (len(self.types) != len(types)):
            raise(self.BadSpecifiedType(f"{type(self).retype} given {len(types)} types, this struct requires {len(self.types)}."))
        for c, (new_t, old_t) in enumerate(zip(types, self.types)):
            if (new_t == old_t): continue
            # Update the stored type
            self.types[c] = new_t
            # Retype all non-missing elements in that column (in place)
            for i in range(len(self)):
                if (type(self[i][c]) == type(None)): continue
                try:
                    self[i][c] = new_t(self[i][c])
                except ValueError:
                    raise(self.BadSpecifiedType(f"Type casting {new_t} for column {c} is not compatible with existing value '{self[i][c]}' on row {i}."))

    # Given a new column of data, add it to this Struct
    def add_column(self, column, name=None):
        # Verify the column length
        if (len(column) != len(self)):
            raise(self.BadData(f"Additional column has length {len(column)}, but this struct has {len(self)} rows."))
        # Verify the name
        if (type(name) == type(None)):
            num = 0
            while (str(num) in self.names): num += 1
            self.names.append(str(num))
        elif (type(name) != str):
            raise(self.BadSpecifiedName(f"Only string names are allowed. Received name '{name}' with type {type(name)}."))
        else:
            self.names.append(name)
        # Verify the column type dynamically. Add new values to all rows.
        new_type = None
        for i,val in enumerate(column):
            self[i].append(val)
            if (type(val) == type(None)):
                # Only add to missing values entry if it's not already there
                if (id(self[i]) not in self.missing):
                    self.missing.add(id(self[i]))
            elif (type(new_type) == type(None)):
                # Capture the new type (type must be right)
                new_type = type(val)
            elif (type(val) != new_type):
                # This is a new type, problem!
                raise(self.BadValue("Provided column has multiple types. Original type {new_type}, but '{val}' has type {type(val)}."))
        # Finally, record the new type
        self.types.append(new_type)

    # Generate a compact numpy structured array from the contents of self.
    def to_numpy_struct(self):
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

    # Generate a pair of functions. The first function will map rows
    # of Struct to a real-valued list. The second function will map
    # real-valued lists back to rows in the original space.
    def generate_mapping(self):
        # Make sure we can handle this type of data
        if not all(t in {int,float,str} for t in self.types):
            raise(self.Unsupported("A mapping is only designed for Structs composed of <float>, <int>, and <str> typed values."))
        # Get those categoricals that need to be converted to reals
        names = self.names[:]
        types = self.types[:]
        cat_names = [n for (n,t) in zip(names,types) if (t == str)]
        cat_values = {n:sorted(set(self[n])-{None}) for n in cat_names}
        cat_mappings = {n:regular_simplex(len(cat_values[n]))
                        for n in cat_names}
        # Calculate the dimension of the real space (and the column names)
        real_names = []
        for (n,t) in zip(names, types):
            if (t in {int, float}):
                real_names += [n]
            else:
                real_names += [f"{n}-{i+1}" for i in range(len(cat_values[n])-1)]
        real_dim = len(real_names)

        # Generate a function that takes a row in the original space
        # and generates a row in the associated real-valued space.
        def to_real(row):
            # Check for row length
            if (len(row) != len(names)):
                raise(self.BadElement(f"Provided row has {len(row)} elements, expected {len(names)}."))
            # Bulid the real-valued vector
            real_row = []
            for (n,t,v) in zip(names, types, row):
                # Check for the type of the elements of the row
                if (t != type(v)):
                    raise(self.BadValue(f"Value '{v}' in provided row has wrong type. Got {type(v)}, expected {t}."))
                if (t == str):
                    # WARNING: Not handling "unrecognized categories"
                    cat_index = cat_values[n].index(v)
                    mapped_vec = cat_mappings[n][cat_index]
                    real_row += list(mapped_vec)
                else:
                    real_row += [float(v)]
            # Return the real-valued vector
            return real_row

        # Generate a function that takes a row in the real-valued
        # space and generates a row in the original space.
        def real_to(real_row):
            import numpy as np
            real_row = list(real_row)
            # Verify length of real-vector
            if (len(real_row) != real_dim):
                raise(self.BadElement(f"Provided real vector has {len(real_row)} elements, expected {real_dim}."))
            # Build the original row
            row = []
            for (n,t) in zip(names, types):
                if (t == str):
                    num_categories = len(cat_values[n])
                    vec = np.array(real_row[:num_categories-1])
                    # Handle reverse-categorical mapping
                    closest_cat = np.argmin(np.sum((cat_mappings[n] - vec)**2,axis=1))
                    category = cat_values[n][closest_cat]
                    row += [category]
                    # Reduce the length of the real-valued row
                    real_row = real_row[num_categories-1:]
                elif (t == int):
                    row += [int(real_row.pop(0) + .5)]
                else:
                    row += [float(real_row.pop(0))]
            return row
        # Return the two mapping functions
        return to_real, real_to, real_names

    # Convert this Struct automatically to a real-valued array.
    # Return real-valued array, function for going to real from
    # elements of this Struct, function for going back to elements of
    # this struct from a real vector.
    def to_numpy_real(self):
        import numpy as np
        to_real, real_to, real_names = self.generate_mapping()
        data = []
        for row in self:
            if None in row: continue
            data.append( to_real(row) )
        array = np.array(data)
        # Container for important post-processing information. 
        # 
        # Info.to_real -- the function for processing regular
        #                 vectors into associated real vectors
        # Info.real_to -- the function for processing real vectors
        #                 back into regular vector
        # Info.names   -- meaningful column names for the real vectors
        class Info: pass
        # Generate an info object and return it with the array.
        info = Info()
        info.to_real = to_real
        info.real_to = real_to
        info.names = real_names
        # Return the numeric array
        return array, info

    # Give an overview (more detailed than just "str") of the contents
    # of this struct. Useful for quickly understanding how data looks.
    def summarize(self, max_display=None):
        # Set the "max_display" to the default value for this class
        if type(max_display) == type(None): max_display = self._max_display
        # Custom counter (dict with default value of 0)
        class Counter(dict):
            def __getitem__(self, val):
                if val in self: return super().__getitem__(val)
                else: return 0
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
        for n,t in zip(self.names, self.types):
            # Count the number of elements for each value
            counts = Counter()
            for val in self[n]:
                counts[val] = counts[val] + 1
            print(f"  {n:{name_len}s} {str(t):{type_len}s} ({len(counts)} unique value{'s' if len(counts) != 1 else ''})")
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
                        num = sum(counts[v] for v in counts if lower <= v < upper)
                        cap = "]"
                    else:
                        num = sum(counts[v] for v in counts if lower <= v <= upper)
                        cap = ")"
                    perc = 100. * (num / len(self))
                    print(f"    [{lower:.2e}, {upper:.2e}{cap} {num:{len(str(len(self)))}d} ({perc:5.1f}%) {'#'*round(perc/2)}")
            else:
                if t in {int, float}:
                    # Order the values by their inate ordering
                    ordered_vals = sorted(counts)
                else:
                    # Order the values by their frequency and print
                    ordered_vals = sorted(counts, key=lambda v: -counts[v])
                val_len = max(map(lambda v: len(str(v)), counts))
                val_len = max(val_len, len("none"))
                if (none_count > 0):
                    perc = 100. * (none_count / len(self))
                    print(f"    {'None':{val_len}s} {none_count:{len(str(len(self)))}d} ({perc:5.1f}%) {'#'*round(perc/2)}")
                for val in ordered_vals:
                    perc = 100. * (counts[val] / len(self))
                    print(f"    {str(val):{val_len}s} {counts[val]:{len(str(len(self)))}d} ({perc:5.1f}%) {'#'*round(perc/2)}")
            print()


# TODO: Rewrite "read_struct" function to use the Struct class. Do not
#       separately read into a python structure and then convert
#       that. Read straight into the Struct. Use type digression to
#       handle unexpected changes in types (when no types are
#       provided). When types are provided and conversion fails, fill
#       with "None" and print warning.

# Some tests for the struct
def test_struct():
    import os
    import numpy as np

    print("Testing Struct...", end=" ")
    a = Struct()

    # Verify append
    a.append([1,"a"])
    a.append([2,"b"])
    a.append([3,"c"])

    # Verify add_column
    a.add_column([1.2, 3.0, 2.4])

    # Verify in-place addition
    b = a[:]
    c = a[1:]
    b += c
    assert(tuple(b["0"]) == tuple([1,2,3,2,3]))

    # Verify slicing-based singleton value assignment
    b = a[:]
    b[:,0] = 1
    assert(tuple(b["0"]) == tuple([1,1,1]))

    # Verify slicing-based multiple value assignment
    b = a[:]
    b[:,0] = [3,1,4]
    assert(tuple(b["0"]) == tuple([3,1,4]))

    # Verify serial indexing
    assert(a[7] == "c")

    # Verify double indexing with integers
    assert(a[0,1] == "a")

    # Verify double indexing with slices
    assert(tuple(a[::-1,1]["1"]) == tuple(["c","b","a"]))

    # Verify that copies with "copy" are deep
    b = a.copy()
    b.retype([str,str,str])
    assert(a.types[0] != str)

    # Verify that copies with slicing are deep
    b = a[:]
    b.retype([str,str,str])
    assert(a.types[0] != str)

    # Verify standard index access
    assert(tuple(a[0]) == tuple([1,"a",1.2]))

    # Verify column-access and automatic column naming
    assert(tuple(a["0"]) == tuple([1,2,3]))

    # Verify slicing by index
    assert(tuple(a[:1][0]) == tuple(a[0]))

    # Verify slicing by names
    assert(tuple(a["0","2"][0]) == tuple([1,1.2]))

    # Verify automatic type generation AND automatic type-casting
    a.append([1,2,3])
    assert(tuple(a[-1]) == tuple([1,"2",3.0]))
    assert(tuple(map(type,a[-1])) == (int, str, float))

    # Verify caching of missing values
    a.append([-1,None,None])
    assert(len(a.missing) == 1)

    # Verify the ability to construct real-valued arrays (and go back)
    b, info = a.to_numpy_real()
    c = Struct(map(info.real_to, b))
    assert(tuple(a["1"][:-1]) == tuple(c["1"]))

    # Verify conversion into numpy structured array (without "None")
    b = a.to_numpy_struct()
    assert(type(b) == np.ndarray)
    assert(type(b.dtype.names) == tuple)
    assert(tuple(b['1']) == tuple(a['1'][:-1]))

    # Verify column item retrieval and assignment
    assert(a["0"][0] == a[0][0])
    a["0"][0] = -10
    assert(a[0][0] == -10)
    a["0"][0] = 1

    # Verify total column assignment
    first_col = list(a["0"])
    new_col = list(map(str,a["0"]))
    a["0"] = new_col
    assert(tuple(a["0"]) == tuple(new_col))
    a["0"] = first_col

    # Verify that missing values are handled correctly by add_column
    b = a[:]
    b.add_column([1,2,3,4,None])
    assert(id(b[-1]) in b.missing)
    assert(len(b.missing) == 1)

    # Verify copying a struct that has a 'None' in the first row
    b = Struct(names=["0","1"], types=[int,str])
    b.append([0,None])
    b.append([1,'b'])
    c = b[:]
    assert(tuple(b[0]) == tuple(c[0]))

    # Verify accessing multiple rows by list-access
    b = a[[0,2,3]]
    assert(tuple(b["0"]) == tuple([1,3,1]))

    # Verify load and save of a csv file
    a.save("a-test.csv")
    b = Struct().load("a-test.csv")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.csv")

    # Verify load and save of a pkl file
    a.save("a-test.pkl")
    b = Struct().load("a-test.pkl")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.pkl")

    # Verify load and save of a gzipped dill file
    a.save("a-test.dill.gz")
    b = Struct().load("a-test.dill.gz")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.dill.gz")

    # Verify load and save of a gzipped csv file
    a.save("a-test.csv.gz")
    b = Struct().load("a-test.csv.gz")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.csv.gz")

    # Verify column equivalence checker
    b = a[a["0"] == 1]
    assert(tuple(b["0"]) == tuple([1,1]))

    # Try providing a generator that does not give integers
    try:
        a[(v for v in ('a','b'))]
    except Struct.BadIndex: pass
    else: assert(False)

    # Try adding a too-short row
    try:
        b = a[:]
        b.append([9,"z"])
    except Struct.BadElement: pass
    else: assert(False)

    # Try a too-short column
    try:
        b = a[:]
        b.add_column([1])
    except Struct.BadData: pass
    else: assert(False)

    # Try adding a name that is not a string
    try:
        b = a[:]
        b.add_column([1,2,3,4,5], name=1)
    except Struct.BadSpecifiedName: pass
    else: assert(False)

    # Try adding a column with multiple types
    try:
        b = a[:]
        b.add_column([1,2,3,4,5.0])
    except Struct.BadValue: pass
    else: assert(False)

    # Try a mismatched-type slice assignment
    try:
        b = a[:]
        b[:,0] = 1.0
    except Struct.BadAssignment: pass
    else: assert(False)

    # Try a mismatched-type column assignment operation
    try:
        a["0"] = list(a["0"])[:-1] + ["0"]
    except Struct.BadAssignment: pass
    else: assert(False)

    # Try a too-short column assignment operation
    try:
        a["0"] = list(map(str,a["0"]))[:-1]
    except Struct.BadAssignment: pass
    else: assert(False)

    # Try a too-long column assignment operation
    try:
        a["0"] = list(map(str,a["0"])) + [1]
    except Struct.BadAssignment: pass
    else: assert(False)

    # Try a bad combination of names and types
    try:
        Struct(names=["a","b","c"], types=[str, str])
    except Struct.NameTypeMismatch: pass
    else: assert(False)

    # Try a bad value in "names"
    try:
        Struct(names=[1,2,3], types=["a","b","c"])
    except Struct.BadSpecifiedName: pass
    else: assert(False)

    # Try a bad value in "types"
    try:
        Struct(names=["a","b","c"], types=[1,2,3])
    except Struct.BadSpecifiedType: pass
    else: assert(False)

    # Try bad value
    try:
        a.append(["a","b","c"])
    except Struct.BadValue: pass
    else: assert(False)

    # Try bad element
    try:
        a.append([1,2])
    except Struct.BadElement: pass
    else: assert(False)

    # Try non-existing column index
    try:
        a["hello"]
    except Struct.UnknownName: pass
    else: assert(False)

    # Try column indexing an uninitialized struct
    try:
        Struct()["a"]
    except Struct.NoNamesSpecified: pass
    else: assert(False)

    # Done testing
    print("passed.")

# =============================================================
#      END OF Python 'Struct' with named and typed columns     
# =============================================================


# Run test cases
if __name__ == "__main__":
    test_struct()
