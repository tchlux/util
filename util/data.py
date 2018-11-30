from util.system import AtomicOpen

#      data.py     
# =================
COMMON_SEPERATORS = [" ", ",", ";", "	"]

UPDATE_FREQUENCY = .1  # How much time (in seconds) between updates.
MAX_ERROR_PRINTOUT = 10 # Only print out this many errors when processing data

import numpy as np
NP_TYPES = {str:(np.str_, 16),    # Which numpy types correspond to 
            int:(np.int64,),      # the python types of variables
            float:(np.float64,)}
PY_TYPES = {value[0]:key for (key,value) in NP_TYPES.items()}
GENERATOR_TYPE = type(_ for _ in ())


# coding: future_fstrings
# from future import print_function

# TODO:  Add tests for empty data object.
# TODO:  Add "Descriptor" class for Data.names and Data.types so that
#        they can be directly modified by users safely (with checks).
# TODO:  Add "Row" class that protects item assignment types over rows
#        and allows for item-wise equality operator.
# TODO:  Non-singleton comparison operator for Column Column comparison.
# TODO:  Add addition, subtraction, multiplication, and division
#        operators to column objects.
# TODO:  Make "predict" not require that target columns be convertible
#        to numeric form, only that they operate under a weighted sum.
# TODO:  Make printout shorter for data with lots of columns.
# TODO:  Check the __iadd__ for data object columns, breaking "missing".
# TODO:  Update __iadd__ to behave differently when adding data:
#        when adding data with subset of known columns, add rows and
#        columns (for unknown) and None for missing.
# TODO:  Make "to_matrix" automatically flatten columns that are typed
#        as numeric numpy arrays.

QUOTES = {'"'}
DEFAULT_DISPLAY_WAIT = 3.

# Define a class for reading a video file without actually storing it in memory.
class IndexedVideo:
    import cv2
    import numpy as np
    class FailedSet(Exception): pass
    class FailedRead(Exception): pass
    class BadIndex(Exception): pass
    class IndexTooLarge(Exception): pass
    # Initialize by getting the information about this video from file.
    def __init__(self, video_file_name, flatten=False, use_float=False):
        # Store the 'flatten' attirbute for flattening images.
        self.flatten = flatten
        self.float = use_float
        # Store the video using an OpenCV video capture object (with its attributes)
        self.vid = self.cv2.VideoCapture(video_file_name)
        self.frames = int(self.vid.get(self.cv2.CAP_PROP_FRAME_COUNT))
        self.width  = int(self.vid.get(self.cv2.CAP_PROP_FRAME_WIDTH))
        self.height = int(self.vid.get(self.cv2.CAP_PROP_FRAME_HEIGHT))
        self.fps    = int(self.vid.get(self.cv2.CAP_PROP_FPS))
        # Find the last available index in the video using the "set"
        # method. This must be done because of a bug in OpenCV.
        self.length = self.frames-1
        for i in range(-1, -self.frames, -1):
            # If the index is successfully accessed, break.
            try:   
                img = self[i]
                self.length = self.frames + i
                break
            # Otherwise, ignore the failed read operation.
            except self.FailedRead: pass
        # Use the first image to get the full shape of this video.
        self.shape = (self.length,) + self[0].shape
        self.file_name = video_file_name

    # Define the "length" of this video.
    def __len__(self): return self.length

    # Define the index operator for this video.
    def __getitem__(self, index):
        # Cover two possible (common) usage errors.
        if type(index) != int: raise(self.BadIndex("Only integer indices allowed."))
        if index > self.length: raise(self.IndexTooLarge(f"Index {index} greater than length of self, {self.length}."))
        # Convert negative indices into positive indices.
        if (index < 0):
            while (index < 0): index += self.length
        # Move the capture to the correct frame of the video.
        success = self.vid.set(self.cv2.CAP_PROP_POS_FRAMES, index)
        if not success: raise(self.FailedSet(f"Failed to set to video frame {index}."))
        # Get the image at that frame.
        success, image = self.vid.read()
        if not success: raise(self.FailedRead(f"Failed to read video frame {index}."))
        # Flatten / float the image if that setting is enabled.
        if self.flatten: image = image.flatten()
        if self.float:   image = self.np.asarray(image, dtype=float)
        # Return the image at that frame.
        return image


# Given "d" categories that need to be converted into real space,
# generate a regular simplex in (d-1)-dimensional space. (all points
# are equally spaced from each other and the origin) This process
# guarantees that all points are not placed on a sub-dimensional
# manifold (as opposed to "one-hot" encoding).
def regular_simplex(num_categories):
    import numpy as np
    class InvalidNumberOfCategories(Exception): pass
    # Special cases for one and two categories
    if num_categories < 1:
        raise(InvalidNumberOfCategories(
            "Number of categories must be an integer greater than 0."))
    elif num_categories == 1:
        return np.array([[]])
    elif num_categories == 2:
        return np.array([[0.],[1.]])
    # Standard case for >2 categories
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

# Return the given element without quotes (if it is a string)
def without_quotes(s): 
    is_str = (type(s) == str)
    if is_str:
        is_long_enough = (len(s) > 1)
        if is_long_enough:
            has_single_quotes = (s[0] == s[-1] == "'")
            has_double_quotes = (s[0] == s[-1] == '"')
            if (has_single_quotes or has_double_quotes):
                return s[1:-1]
    return s


# Remove nested seperators inside of quotes in the line
def remove_quoted_seperators(line, sep=None):
    cleaned_line = ""
    in_quotes = False
    for char in line:
        # Toggle 'inside' quotes.
        if (char in QUOTES):
            in_quotes = not in_quotes
            if (sep == None): char = ""
        # Add quoted bits without separ
        if (not in_quotes):
            cleaned_line += char
        elif (char != sep) and (sep != None):
            cleaned_line += char
    return cleaned_line


# Pre:  "filename" is the name of an existing file (including path)
# Post: "filename" is scanned for viable seperators, those are all
#       characters used the same number of times on each line. If at
#       the end of the file there are still multiple candidate
#       seperators, they are sorted by most-occuring first. If one of
#       COMMON_SEPERATORS exist in the set, it is used, otherwise the
#       most occurring seperator is used.
def detect_seperator(filename="<no_provided_file>", verbose=False, opened_file=None):
    import time
    lines = []
    not_viable = set()
    if not opened_file:
        f = AtomicOpen(filename).file
    else:
        f = opened_file
    # Identify the potential seperators from the first line
    first_line = remove_quoted_seperators(f.readline().strip())
    seperators = {char:first_line.count(char) for char in set(first_line)}
    # Update the user
    last_update = 0
    should_update = lambda: (verbose and ((time.time() - last_update) > UPDATE_FREQUENCY))
    raw_data = f.readlines()
    # Cycle the file (narrowing down the list of potential seperators
    for i,line in enumerate(raw_data):
        # Update the user on progress if appropriate
        if should_update():
            last_update = time.time()
            print(f"\r{(100.0*i/len(raw_data)):0.1f}% complete", flush=True, end="")
        # Clean the line of quoted strings.
        line = remove_quoted_seperators(line)
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
# Post: A data.Data class containing all data from the file.
def read_data(filename="<no_provided_file>", sep=None, types=None,
              verbose=False, opened_file=None):
    import time
    if sep == None: sep = detect_seperator(filename, verbose, opened_file)
    # Get the opened file object
    if (not opened_file):
        # Open the file to get the data
        f = AtomicOpen(filename).file
    else:
        f = opened_file
        file_name = f.name
    # Get the first line, assuming it's the header
    line = remove_quoted_seperators(f.readline(), sep)
    header = [without_quotes(n.strip()).strip() for n in line.strip().split(sep)]
    # Get the first line of data and learn the types of each column
    line = [without_quotes(v.strip()).strip() for v in 
            remove_quoted_seperators(f.readline().strip(), sep).split(sep)]
    # Automatically detect types (if necessary)
    type_digression = False
    if type(types) == type(None):
        types = list(map(get_type, line))
        type_digression = True
    # Initialize our data holder
    data = Data(names=header, types=types)
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
        # Clean up the line (in case there are nested seperators)
        line = remove_quoted_seperators(line, sep)
        # Update the user on progress if appropriate
        if should_update():
            last_update = time.time()
            print("\r%0.1f%% complete"% (100.0*i/len(raw_data)), flush=True, end="")
        # Read the line of data into a list format (replace empty string with None)
        list_line = line.strip().split(sep)
        # Remove leading and trailing quotes from strings if necessary
        list_line = [without_quotes(el).strip() for el in list_line]
        # Replace any empty strings with "None" for missing values
        while "" in list_line: list_line[list_line.index("")] = None
        # Try and add the line of data
        try:
            data.append(list_line)
        except Data.BadValue:
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
    # Return the data
    return data

# ======================================================
#      Python 'Data' with named and typed columns     
# ======================================================

# Define a python data matrix structure (named and typed columns)
class Data(list):
    # Default values for 'self.names' and 'self.types'
    names = None
    types = None
    # Default values for "numeric" representation of self.
    numeric = None
    # Default maximum printed rows of this data.
    _max_display = 10
    # Default maximum string length of printed column values.
    _max_str_len = 50

    # Local descriptive exceptions
    class NameTypeMismatch(Exception): pass
    class BadSpecifiedType(Exception): pass
    class BadSpecifiedName(Exception): pass
    class NoNamesSpecified(Exception): pass
    class ImproperUsage(Exception): pass
    class UnknownName(Exception): pass
    class Unsupported(Exception): pass
    class BadAssignment(Exception): pass
    class BadElement(Exception): pass
    class BadTarget(Exception): pass
    class BadValue(Exception): pass
    class BadIndex(Exception): pass
    class BadData(Exception): pass
    class Empty(Exception): pass
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
                index = range(len(self.indices))[index]
            elif (hasattr(index, "__iter__")):
                pass
            else:
                raise(self.data.BadIndex(f"Data column does not support indexing with type {type(index)}."))
            # Convert provided index into a list of integers.
            indices = []
            for i in index:
                if (type(i) != int):
                    raise(self.data.BadIndex(f"Data column index iterable contained type {type(i)}, expected {int}."))                    
                indices.append( self.indices[i] )
            # Return a new column with modified indices.
            return self.data.Column(self.data, self.column, indices)
        def __len__(self): return len(self.indices)
        # Define this as an interator
        def __iter__(self):
            self.index = 0
            return self
        def __next__(self):
            if (self.index >= len(self.indices)):
                raise(StopIteration)
            self.index += 1
            return self.data[self.indices[self.index-1], self.column]
        # Use the iteration to generate a string of this column
        def __str__(self):
            return str(list(self))
        # Inequality operator
        def __ne__(self, group):
            return self.__eq__(group, equality=False)
        # Iterate over the indices that exist in a set of values
        def __eq__(self, group, equality=True):
            # Try to convert given group to an iterable type for comparison.
            if ((type(group) == self.data.types[self.column]) or (type(group) == type(None))):
                group = [group]
            elif (type(group) in {set, list}): 
                pass
            else:
                raise(self.data.BadData(f"There is no defined equality operator for the provided type {type(group)}, when it does not match expected type {self.data.types[self.column]}."))
            # Iterate over rows, generating indices that match the equality.
            for i,v in enumerate(self):
                if ((equality) and (v in group)):
                    yield i
                elif ((not equality) and (v not in group)):
                    yield i
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
    class Row(Column):
        pass

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
            if any(type(t) != type for t in types):
                raise(self.BadSpecifiedType(f"An entry of provided types was not of {type(type)} class."))
            self.types = types[:]
        self.missing = set()
        # If data was provided, put it into this data
        if type(data) != type(None):
            for row in data:
                self.append(row)

    # Redefined append operator that 
    def append(self, element):
        # Make sure the element has a "len" property
        try:    len(element)
        except: raise(self.BadValue(f"Invalid appended element of type {type(element)}, which does not have 'len'."))
        # Convert the element into a python list
        if type(element) != list:
            try:    element = list(element)
            except: raise(self.BadValue(f"Invalid appended element, failed conversion to type {list}."))
        missing_values = False
        # Check length in case names already exist
        if (type(self.names) != type(None)):
            if (len(element) != len(self.names)):
                raise(self.BadElement(f"Only elements of length {len(self.names)} can be added to this data."))
        # Check length in case types already exist
        if (type(self.types) != type(None)):
            if len(element) != len(self.types):
                raise(self.BadElement(f"Only elements of length {len(self.types)} can be added to this data."))                
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
        super().append(list(element))
        # Add to the list of missing values if necessary
        if (missing_values): self.missing.add(id(self[-1]))
        # If we are storing a numeric representation and there are no missing values
        elif (type(self.numeric) != type(None)):
            import numpy as np
            # Get the real representation of the new row
            new_row = np.array(self.numeric.to_real(self[-1]))
            # Append the new row to the matrix of data.
            self.numeric.data = np.concatenate((self.numeric.data,
                                                new_row[None,:]))

    # Overwrite the standard "[]" get operator to accept strings as well.
    def __getitem__(self, index):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            raise(self.Empty("Cannot get item from empty data."))
        if ((type(index) == tuple) and (len(index) == 2) and 
            (all(type(i) == int for i in index))):
            # This index is accessing a (row,col) entry
            return super().__getitem__(index[0])[index[1]]
        elif ((type(index) == tuple) and (len(index) == 2) and 
              (type(index[0]) == int and type(index[1]) == str)):
            if (index[1] not in self.names):
                raise(self.BadIndex(f"Column index {index[1]} is not a recognized column name."))
            # This index is accessing a (row,col) entry
            return super().__getitem__(index[0])[self.names.index(index[1])]
        elif (type(index) == int):
            # This index is accessing a Data row
            return super().__getitem__(index)
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
            names = [self.names[col] for col in cols]
            types = [self.types[col] for col in cols]
            new_data = Data(names=names, types=types)
            for row in rows:
                if type(row) != int:
                    raise(self.BadIndex(f"The provided row index of type {type(row)}, {row}, is not understood."))
                new_data.append( [self[row,col] for col in cols] )
            return new_data

    # Overwrite the standard "[]" get operator to accept strings as well.
    def __setitem__(self, index, value):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            raise(self.Empty("Cannot set item for empty data."))
        self.numeric = None # Reset the numeric representation (if it existed)
        # Special case assignment of new column
        if ((type(index) == str) and (index not in self.names)):
            return self.add_column(value, name=index)
        # Get the list of rows and list of columns being assigned.
        rows, cols = self._index_to_rows_cols(index)
        # Assume if it has an iterator that is not a singleton.
        # Assume that if the value type matches all assigned types it is a singleton.
        singleton = (not hasattr(value, '__iter__')) or \
                    (all(self.types[c] in {type(value), type(None)} for c in cols))
        # If it has a "length" then use that to verify singleton status.
        if (singleton and hasattr(value, "__len__")):
            singleton = len(value) != (len(rows)*len(cols))
        # Assume singleton status is correct.
        if not singleton: value_iter = value.__iter__()
        # Iterate over rows and perform assignment
        step = 0
        new_type = type(None)
        for row in rows:
            if type(row) != int:
                raise(self.BadIndex(f"The provided row index of type {type(row)}, {row}, is not understood."))
            for col in cols:
                # Retreive next value from iterator if necessary.
                if not singleton:
                    try:
                        value = value_iter.__next__()
                    except StopIteration:
                        raise(self.BadValue(f"Provided iterable only contained {step} elements, expected more."))
                # Handle type of new value appropriately
                if (len(rows) == len(self)):
                    # Check to find the type of the new column (if appropriate)
                    if (new_type == type(None)):
                        new_type = type(value)
                    # Check the type of the new value being assigned
                    elif (type(value) != new_type):
                        raise(self.BadAssignment(f"Provided value {value} of type {type(value)} does not match expected type {new_type}."))
                # The existing column doesn't have a type, assign it.
                elif (self.types[col] == type(None)):
                    self.types[col] = type(value)
                # Mark this row as contianing a missing value if assigned 'None'.
                elif (type(value) == type(None)):
                    if (id(self[row]) not in self.missing):
                        self.missing.add(id(self[row]))
                # Check against the existing column type (for verification).
                elif (type(value) != self.types[col]):
                    raise(self.BadAssignment(f"Provided value {value} of type {type(value)} does not match expected type {self.types[col]}."))
                # Make the assignment
                super().__getitem__(row)[col] = value
                # Remove entry from missing if appropriate.
                if (type(value) != type(None)):
                    if (id(self[row]) in self.missing):
                        if (None not in self[row]):
                            self.missing.remove(id(self[row]))
                step += 1
        # Update column type if it was totally reassigned.
        if (len(rows) == len(self)):
            self.types[cols[0]] = new_type
        # Verify that the iterable has been exhausted.
        if (not singleton):
            try:
                value = value_iter.__next__()
                raise(self.BadAssignment(f"Column assignment requires a column of length {len(self)}, provided values iterable was too long."))
            except StopIteration: pass

    # Custom function for mapping values to strings in the printout of self.
    def _val_to_str(self, v): 
        # Get the string version of the value.
        if not ((type(v) == type(self)) or (issubclass(type(v),type(self)))):
            string = str(v).replace("\n"," ")
        else:
            # Custom string for a 'data' type object.
            string = v.__str__(short=True)
        # Shorten the string if it needs to be done.
        if len(string) > self._max_str_len: 
            string = string[:self._max_str_len]+".."
        return string

    # Printout a brief table-format summary of this data.
    def __str__(self, short=False):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            return "This Data has no contents.\n"
        # Custom short string for a 'data' type object.
        if short: return f"Data ({len(self.names)} x {len(self)}) -- {self.names}"
        # Make a pretty table formatted output
        rows = []
        rows += [ self.names ]
        rows += [ list(map(lambda t: str(t)[8:-2],self.types)) ]
        for i,row in enumerate(self):
            rows += [ row ]
            if i >= self._max_display: break
        rows = [list(map(self._val_to_str,r)) for r in rows]
        lens = [max(len(r[c]) for r in rows)
                for c in range(len(self.names))]
        rows = [[v + " "*(lens[i]-len(v)) for i,v in enumerate(row)]
                for row in rows]
        rows = [" " + (" | ".join(row)) + "\n" for row in rows]
        string = "\n" + "="*len(rows[0]) + "\n"
        string += f"Size: ({len(self)} x {len(self.names)})\n\n"
        string += rows[0]
        string += rows[1]
        string += "-"*len(rows[0]) + "\n"
        for row in rows[2:]:                string += row
        if (len(self) > self._max_display): string += "  ...\n"
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
                if len(indices) > self._max_display:
                    break
            # Print out the indices of the rows with missing values
            print_indices = ', '.join(map(str, indices[:self._max_display]))
            string += f" missing values ({len(self.missing)}) at rows: \n  [{print_indices}"
            # Add elipses if there are a lot of missing values.
            if (len(indices) > self._max_display): string += ", ..."
            string += "]\n"
        string += "="*len(rows[0]) + "\n"
        return string

    # Return True if empty.
    def empty(self):
        return ((type(self.names) == type(None)) or 
                (type(self.types) == type(None)))

    # Define a convenience funciton for concatenating another
    # similarly typed and named data to self.
    def __iadd__(self, data):
        self.numeric = None # Reset the numeric representation (if it existed)
        # Check for improper usage
        if type(data) != Data:
            raise(self.Unsupported(f"In-place addition only supports type {Data}, but {type(data)} was given."))
        # Special case for being empty.
        if data.empty(): 
            return self
        # Special case for being empty.
        elif self.empty():
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

    # Make sure the copy of this data is a deep copy.
    def copy(self): return self[:]

    # Overwrite the 'pop' method to track missing values.
    def pop(self, index=-1):
        # Popping of a column
        if type(index) == int:
            if (id(self[index]) in self.missing):
                self.missing.remove( id(self[index]) )
            return super().pop(index)
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
        if (id(index) in self.missing):
            self.missing.remove( id(index) )
        return super().remove( index )


    # ========================
    #      Custom Methods     
    # ========================

    # Generate hash values using a near-zero probability of conflich
    # hashing mechanism that can handle non-hashable types. Does this
    # by hashing the raw bytes (from pickl) with sha256.
    def _hash(self, something):
        import hashlib, pickle
        return hashlib.sha256(pickle.dumps(something)).hexdigest()

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
        elif (type(index) == str):
            rows = list(range(len(self)))
            if index not in self.names:
                raise(self.UnknownName(f"This data does not have a column named '{index}'."))
            cols = [self.names.index(index)]
        elif ((type(index) == tuple) and (len(index) == 2) and 
            (all(type(i) == int for i in index))):
            # This index is accessing a (row,col) entry
            rows = [index[0]]
            cols = [index[1]]
        elif (type(index) == slice):
            # Construct a new "data" with just the specified rows (deep copy)
            rows = list(range(len(self))[index])
            cols = list(range(len(self.names)))
        elif (hasattr(index, "__iter__")):
            index = tuple(index)
            # Special case for when a list of ints is used to access rows
            if all(type(i)==int for i in index):
                rows = index
                cols = list(range(len(self.names)))
            # Special case for when a list of strings is used to access columns
            elif all(type(i)==str for i in index):
                rows = list(range(len(self)))
                cols = []
                for i in index:
                    if i not in self.names:
                        raise(self.UnknownName(f"This data does not have a column named '{i}'."))
                    cols.append( self.names.index(i) )
            else:
                # Otherwise this should be a two-index for [row,col] access.
                if (len(index) != 2):
                    raise(self.BadIndex(f"The provided index, {index}, is not understood."))
                # Create row index
                if type(index[0]) == int:
                    rows = [index[0]]
                elif type(index[0]) == slice:
                    rows = list(range(len(self))[index[0]])
                elif hasattr(index[0], "__iter__"):
                    rows = list(index[0])
                else:
                    raise(self.BadIndex(f"The provided index, {index}, is not understood."))
                # Create column index
                if type(index[1]) == int:
                    cols = [index[1]]
                elif type(index[1]) == str:
                    if (index[1] not in self.names):
                        raise(self.UnknownName(f"This data does not have a column named '{index[1]}'."))
                    cols = [self.names.index(index[1])]
                elif hasattr(index[1], "__iter__"):
                    cols = []
                    for i in index[1]:
                        if (type(i) not in {str,int}):
                            raise(self.BadIndex(f"The provided column index of type {type(i)}, {i}, is not understood."))
                        elif (type(i) == str):
                            if (i not in self.names):
                                raise(self.UnknownName(f"This data does not have a column named '{i}'."))
                            i = self.names.index(i)
                        cols.append( i )
                elif type(index[1]) == slice:
                    cols = list(range(len(self.names))[index[1]])
                else:
                    raise(self.BadIndex(f"The provided index, {index}, is not understood."))
        else:
            raise(self.BadIndex(f"The provided index, {index}, is not understood."))
        return rows, cols

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
                read_data_kwargs["opened_file"] = f
                self += read_data(**read_data_kwargs)
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

    # Given a new list of types, re-cast all elements of columns with
    # changed types into the newly specified type.
    def retype(self, types, columns=None):
        # Special case for being empty.
        if ((type(self.names) == type(None)) or 
            (type(self.types) == type(None))):
            raise(self.Empty("Cannot retype empty data."))
        self.numeric = None # Reset the numeric representation (if it existed)
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
    def add_column(self, column, name=None):
        self.numeric = None # Reset the numeric representation (if it existed)
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
        new_type = type(None);
        for i,val in enumerate(column):
            # Verify valid index first..
            if (i >= len(self)):
                # Remove the added elements if the length was not right
                for j in range(len(self)): self[j].pop(-1)
                # Raise error for too long of a column
                raise(self.BadData(f"Provided column has more elements than {len(self)}, the length of this data."))
            # Append the value to this row..
            self[i].append(val)
            if (type(val) == type(None)):
                # Only add to missing values entry if it's not already there
                if (id(self[i]) not in self.missing):
                    self.missing.add(id(self[i]))
            elif (new_type == type(None)):
                # Capture the new type (type must be right)
                new_type = type(val)
            elif (type(val) != new_type):
                # This is a new type, problem!
                raise(self.BadValue(f"Provided column has multiple types. Original type {new_type}, but '{val}' has type {type(val)}."))
        # Verify the column length
        if (i < len(self)-1):
            # Remove the added elements if the length was not right
            for j in range(i+1): self[j].pop(-1)
            # Raise error for too short of a column
            raise(self.BadData(f"Provided column has length less than {len(self)}, the length of this data."))
        # Finally, record the new type
        self.types.append(new_type)

    # Generate a pair of functions. The first function will map rows
    # of Data to a real-valued list. The second function will map
    # real-valued lists back to rows in the original space.
    def generate_mapping(self):
        if self.empty(): raise(self.ImproperUsage("Cannot map empty Data."))
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

    # Convert this Data automatically to a real-valued array.
    # Return real-valued array, function for going to real from
    # elements of this Data, function for going back to elements of
    # this data from a real vector.
    def to_matrix(self, print_column_width=70):
        import numpy as np
        to_real, real_to, real_names, num_inds, cat_inds = self.generate_mapping()
        data = []
        for row in self:
            # How to handle rows with missing values?
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

    # Collect the dictionaries of unique values (with counts) for each column.
    def counts(self):
        column_info = {n:{} for n in self.names}
        for row in self:
            for n,val in zip(self.names, row):
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
        unique = Data(names=self.names, types=self.types)
        found = set()
        start = time.time()
        for i,row in enumerate(self):
            # Update user on progress if too much time has elapsed..
            if (time.time() - start) > display_wait_sec:
                print(f" {100.*i/len(self):.2f}%", end="\r", flush=True)
                start = time.time()
            # Get the hashed value, check to see if its been found.
            hashed = self._hash(row)
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
            hash_value = self._hash([row[c] for c in self_columns])
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
            hash_value = self._hash([row[c] for c in other_columns])
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
        random.seed(seed)
        indices = list(range(len(self)))
        random.shuffle(indices)
        # Re-seed the random number generator as not to muck with
        # other processes that might be using it.
        random.seed()
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
        # Get the numeric representation of self.
        if type(self.numeric) == type(None):
            self.numeric = self.to_matrix()
        # Pick the filling model based on data size.
        if (type(model) == type(None)):
            # Use Voronoi if there are fewer than 10000 points.
            if (len(self) <= 10000):
                from util.algorithms import Voronoi, uniqueify
                model = uniqueify(Voronoi)()
            # Otherwise use nearest neighbor (can't get points + weights from NN)
            else:
                from util.algorithms import NearestNeighbor, uniqueify
                model = uniqueify(NearestNeighbor)()
        # Set the weights accordingly
        if (type(weights) == type(None)):
            weights = np.ones(self.numeric.data.shape[1])
        elif (len(weights) != len(self.names)):
            raise(self.BadValue("Weights vector is length {len(weights)}, expected length {len(self.names)}."))
        else:
            # Convert the provided weight to the correct shape.
            col_weights, weights = weights, []
            col_inds = list(range(len(self.numeric.names)))
            while (len(col_inds) > 0):
                if (col_inds[0] not in self.numeric.cats):
                    col_inds.pop(0)
                    weights.append( col_weights.pop(0) )
                else:
                    for i in range(1,len(col_inds)):
                        not_categorical = (col_inds[i] not in self.numeric.cats)
                        beginning_of_next = ("1" == 
                            self.numeric.names[col_inds[i]].split('-')[-1])
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
        numeric_row = self.numeric.to_real(row, fill=True)
        # Identify the indices of known values and unknown values.
        known_values = np.array([i for i in range(len(numeric_row))
                                 if (type(numeric_row[i]) != type(None))])
        unknown_values = np.array([i for i in range(len(numeric_row))
                                   if (type(numeric_row[i]) == type(None))])
        # Get the part of the numeric row
        no_none_numeric_row = np.array([numeric_row[i] for i in known_values])
        # Normalize the numeric row (to debias numeric scale differences)
        no_none_numeric_row = ((no_none_numeric_row - self.numeric.shift[known_values])
                               / self.numeric.scale[known_values]) * weights[known_values]
        # Get the normalized training point locations.
        normalized_numeric_data = ((self.numeric.data - self.numeric.shift)
                                   / self.numeric.scale) * weights
        # Fit the model and make a prediction for this row.
        if hasattr(model, "WeightedApproximator"):
            # Fit the model (just to the points)
            model.fit(normalized_numeric_data[:,known_values])
            del(normalized_numeric_data)
            # Get the source points and source weights for the prediction
            pts, wts = model.predict(np.array(no_none_numeric_row))
            # Compute the fill row by performing the weighted sum.
            fill_row = np.sum((self.numeric.data[pts,:][:,unknown_values].T * wts), axis=1)
            # Convert points and weights into python list (in case they're not)
            pts = list(pts)
            wts = list(wts)
        else:
            # The model makes predictions using all points.
            # Fit the model to the input points and output points.
            model.fit(normalized_numeric_data[:,known_values], 
                      normalized_numeric_data[:,unknown_values])
            del(normalized_numeric_data)
            pts = list(range(self.numeric.data.shape[0]))
            wts = [1 / len(pts)] * len(pts)
            # Use the model to make a prediction.
            fill_row = model.predict(no_none_numeric_row)
            # De-normalize the predicted output.
            fill_row = ((fill_row / weights[unknown_values]) * 
                        self.numeric.scale[unknown_values] - 
                        self.numeric.shift[unknown_values])
        # Pad fill row with 0's so that it is the correct length for conversion.
        fill_row = [0.] * len(known_values) + list(fill_row)
        # python3 -m pytest data.py
        # Transfer the filled numeric row into the original row
        fill_row = self.numeric.real_to(fill_row, bias=bias)
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

    # Given a single column of this data to try and predict, run a
    # k-fold cross validation over predictions using 'fill'.
    def predict(self, target_columns, model=None, weights=None, bias={}, k=10):
        import numpy as np
        # Get the numeric representation of self.
        if (type(self.numeric) == type(None)):
            self.numeric = self.to_matrix()
        # Pick the filling model based on data size.
        if (type(model) == type(None)):
            # Use Voronoi if there are fewer than 10000 points.
            if (len(self) <= 10000):
                from util.algorithms import Voronoi, uniqueify
                model = uniqueify(Voronoi)()
            # Otherwise use nearest neighbor (can't get points + weights from NN)
            else:
                from util.algorithms import NearestNeighbor, uniqueify
                model = uniqueify(NearestNeighbor)()
        # Set the weights accordingly
        if (type(weights) == type(None)):
            weights = np.ones(self.numeric.data.shape[1])
        elif (len(weights) != len(self.names)):
            raise(self.BadValue("Weights vector is length {len(weights)}, expected length {len(self.names)}."))
        else:
            # Convert the provided weight to the correct shape.
            col_weights, weights = weights, []
            col_inds = list(range(len(self.numeric.names)))
            while (len(col_inds) > 0):
                if (col_inds[0] not in self.numeric.cats):
                    col_inds.pop(0)
                    weights.append( col_weights.pop(0) )
                else:
                    for i in range(1,len(col_inds)):
                        not_categorical = (col_inds[i] not in self.numeric.cats)
                        beginning_of_next = ("1" == 
                            self.numeric.names[col_inds[i]].split('-')[-1])
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
        normalized_data = weights * (self.numeric.data - self.numeric.shift) / self.numeric.scale
        # Identify those columns that are being predicted
        target_column_indices = [self.names.index(c) for c in target_columns]
        results = Data(names=[c + " Predicted" for c in target_columns]+["Prediction Index"],
                       types=[self.types[i] for i in target_column_indices]+[int])
        # Identify the numeric ranges of those columns
        indices = set(target_column_indices)
        sample = self.numeric.to_real([
            (v if i in indices else None)
            for i,v in enumerate(self[0])])
        del(indices)
        train_cols = np.array([i for (i,v) in enumerate(sample) if type(v) == type(None)])
        test_cols =  np.array([i for (i,v) in enumerate(sample) if type(v) != type(None)])
        del(sample)
        # Cycle through and make predictions.
        for i,(train_rows, test_rows) in enumerate(
                Data.k_fold(range(self.numeric.data.shape[0]), k=k)):
            print(f"  Fold {i+1} of {k}.. ", end="", flush=True)
            train_data = normalized_data[train_rows]
            print("fitting.. ", end="", flush=True)
            model.fit( normalized_data[train_rows,:][:,train_cols], 
                       normalized_data[train_rows,:][:,test_cols]   )
            print("predicting.. ", end="", flush=True)
            predictions = model.predict( normalized_data[test_rows,:][:,train_cols] )
            # De-normalize predictions
            predictions = ( (predictions / weights[test_cols]) * 
                            self.numeric.scale[test_cols] + 
                            self.numeric.shift[test_cols] )
            print("storing.. ", end="\r", flush=True)
            # Store the prediction results.
            for j,i in enumerate(test_rows):
                row = list(self.numeric.data[i,train_cols]) + list(predictions[j,:])
                row = self.numeric.real_to(row, bias=bias)
                results.append( [row[i] for i in target_column_indices] + [i] )
            print(" "*70,end="\r",flush=True)
        # Sort results by the original index
        results.sort(key=lambda r: r[-1])
        results.pop("Prediction Index")
        # Return the results in a new Data object.
        return results

    # Give an overview (more detailed than just "str") of the contents
    # of this data. Useful for quickly understanding how data looks.
    def summarize(self, max_display=None):
        # Special case for an empty data
        if len(self) == 0:
            print(self)
            return
        # Set the "max_display" to the default value for this class
        if type(max_display) == type(None): max_display = self._max_display
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
            for val in self[n]:
                if type(val) in {list, dict, set}: val = str(val)
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
                val_len = min(val_len, self. _max_str_len)
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


# Some tests for the data
def test_data():
    import os
    import numpy as np

    print("Testing Data...", end=" ")
    a = Data()

    # Verify append
    a.append([1,"a"])
    a.append([2,"b"])
    a.append([3,"c"])

    # Verify add_column

    # Verify add column (with edge case, none type)
    a.add_column([None,None,None])
    assert(a.types[-1] == type(None))
    # Reassign the missing value column to floats, verify type update
    a[a.names[-1]] = [1.2, 3.0, 2.4]
    assert(a.types[-1] == float)
    # WARNING: Need to check doing (<int>, <str>) assignment

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

    # Verify double indexing with integers
    assert(a[0,1] == "a")

    # Verify double indexing with slices
    assert(tuple(a[::-1,1]["1"]) == tuple(["c","b","a"]))

    # Verify standard index access
    assert(tuple(a[0]) == tuple([1,"a",1.2]))

    # Verify column-access and automatic column naming
    assert(tuple(a["0"]) == tuple([1,2,3]))

    # Verify slicing by index
    assert(tuple(a[:1][0]) == tuple(a[0]))

    # Verify slicing by names
    assert(tuple(a["0","2"][0]) == tuple([1,1.2]))

    # Verify that copies with "copy" are deep
    b = a.copy()
    b.retype([str,str,str])
    assert(a.types[0] != str)

    # Verify that copies with slicing are deep
    b = a[:]
    b.retype([str,str,str])
    assert(a.types[0] != str)

    # Verify automatic type generation AND automatic type-casting
    a.append([1,2,3])
    assert(tuple(a[-1]) == tuple([1,"2",3.0]))
    assert(tuple(map(type,a[-1])) == (int, str, float))

    # Verify caching of missing values
    a.append([-1,None,None])
    assert(len(a.missing) == 1)

    # Verify the ability to construct real-valued arrays (and go back)
    numeric = a.to_matrix()
    c = Data(map(numeric.real_to, numeric.data))
    assert(tuple(a[:-1]["1"]) == tuple(c["1"]))
    
    # Verify conversion into numpy structured array (without "None")
    b = a.to_struct()
    assert(type(b) == np.ndarray)
    assert(type(b.dtype.names) == tuple)
    assert(tuple(b['1']) == tuple(a['1'][:-1]))

    # Verify column item retrieval and assignment
    assert(a["0"][0] == a[0][0])
    a["0"][0] = -10
    a["0"][0] = 1

    # Verify total column assignment
    first_col = list(a["0"])
    new_col = list(map(str,a["0"]))
    a["0"] = new_col
    assert(tuple(a["0"]) == tuple(new_col))
    a["0"] = first_col

    # Verify new column assignment
    b = a[:]
    first_col = list(b["0"])
    new_col = list(map(str,b["0"]))
    b["0-str"] = new_col
    assert(tuple(map(str,b["0"])) == tuple(b["0-str"]))

    # Verify that missing values are handled correctly by add_column
    b = a[:]
    b.add_column([1,2,3,4,None])
    assert(id(b[-1]) in b.missing)
    assert(len(b.missing) == 1)

    # Verify copying a data object that has a 'None' in the first row
    b = Data(names=["0","1"], types=[int,str])
    b.append([0,None])
    b.append([1,'b'])
    c = b[:]
    assert(tuple(b[0]) == tuple(c[0]))

    # Verify accessing multiple rows by list-access
    b = a[[0,2,3]]
    assert(tuple(b["0"]) == tuple([1,3,1]))

    # Verify load and save of a csv file
    a.save("a-test.csv")
    b = Data().load("a-test.csv")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.csv")

    # Verify load and save of a pkl file
    a.save("a-test.pkl")
    b = Data.load("a-test.pkl")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.pkl")

    # Verify load and save of a gzipped dill file
    a.save("a-test.dill.gz")
    b = Data().load("a-test.dill.gz")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.dill.gz")

    # Verify load and save of a gzipped csv file
    a.save("a-test.csv.gz")
    b = Data().load("a-test.csv.gz")
    assert(tuple(a["0"]) == tuple(b["0"]))
    os.remove("a-test.csv.gz")

    # Verify column equivalence / equality checker by element
    b = a[a["0"] == 1]
    assert(tuple(b["0"]) == tuple([1,1]))

    # Verify one of the inequality operators
    b = a[a["0"] < 2]
    assert(tuple(b["0"]) == tuple([1,1,-1]))

    # Verify column membership checker by set
    b = a[a["1"] == {'a','b','2'}]
    assert(tuple(b["1"]) == tuple(['a','b','2']))

    # WARNING: Not individually verifying *all* comparison operators.
    # WARNING: No tests for list-based index assignment
    # WARNING: No tests for generator-based index assignment

    # Verify the generation of a random k-fold cross validation.
    b = Data()
    for i in range(37):
        b.append([i])
    # Collect the list of all the "testing" rows seen
    all_test = []
    for (train, test) in b.k_fold(k=11, seed=1):
        all_test += list(test["0"])
    assert(sorted(all_test) == list(range(len(b))))

    # Verify the identification of unique rows and collection process.
    b = a[:]
    b += a
    cols = ['0','1']
    b = b[cols].unique().collect(b)
    assert(tuple(b[0,-1]) == tuple(2*[a[0,-1]]))

    # Test the 'fill' method.
    b = a[:]
    b[-2,1] = 'd'
    b.append( b[0] )
    assert(str(b.fill(b[-2]).data) == str([-1,'a',1.7999999999999998]))
    assert(str(a.fill(a[-1]).data) == str([-1,'2',2.1]))

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

    # Attempt to index an empty data
    b = Data()
    try:   b[0]
    except Data.Empty: pass
    else: assert(False)

    # Try providing a generator that does not give strictly integers
    try:   a[(v for v in ('a',1,'b'))]
    except Data.BadIndex: pass
    else: assert(False)

    # Try adding a too-short row
    try:
        b = a[:]
        b.append([9,"z"])
    except Data.BadElement: pass
    else: assert(False)

    # Try a too-short column
    try:
        b = a[:]
        b.add_column([1])
    except Data.BadData: pass
    else: assert(False)

    # Try a too-long column
    try:
        b = a[:]
        b.add_column(map(float,range(1000)))
    except Data.BadData: pass
    else: assert(False)

    # Try adding a name that is not a string
    try:
        b = a[:]
        b.add_column([1,2,3,4,5], name=1)
    except Data.BadSpecifiedName: pass
    else: assert(False)

    # Try adding a column with multiple types
    try:
        b = a[:]
        b.add_column([1,2,3,4,5.0])
    except Data.BadValue: pass
    else: assert(False)

    # Try a mismatched-type slice assignment
    try:
        b = a[:]
        b[:,0] = [1.0, 1, 1, 1, 1]
    except Data.BadAssignment: pass
    else: assert(False)

    # Try a mismatched-type column assignment operation
    try:   a["0"] = list(a["0"])[:-1] + ["0"]
    except Data.BadAssignment: pass
    else: assert(False)

    # Try a too-short column assignment operation
    try:   a["0"] = list(map(str,a["0"]))[:-1]
    except Data.BadValue: pass
    else: assert(False)

    # Try a too-long column assignment operation
    try:   a["0"] = list(map(str,a["0"])) + [1]
    except Data.BadAssignment: pass
    else: assert(False)

    # Try a bad combination of names and types
    try:   Data(names=["a","b","c"], types=[str, str])
    except Data.NameTypeMismatch: pass
    else: assert(False)

    # Try a bad value in "names"
    try:   Data(names=[1,2,3], types=["a","b","c"])
    except Data.BadSpecifiedName: pass
    else: assert(False)

    # Try a bad value in "types"
    try:   Data(names=["a","b","c"], types=[1,2,3])
    except Data.BadSpecifiedType: pass
    else: assert(False)

    # Try bad value
    try:   a.append(["a","b","c"])
    except Data.BadValue: pass
    else: assert(False)

    # Try bad element
    try:   a.append([1,2])
    except Data.BadElement: pass
    else: assert(False)

    # Try non-existing column index
    try:   a["hello"]
    except Data.UnknownName: pass
    else: assert(False)

    # Try column indexing an uninitialized data
    try:   Data()["a"]
    except Data.Empty: pass
    else: assert(False)

    # Done testing
    print("passed.")

# =============================================================
#      END OF Python 'Data' with named and typed columns     
# =============================================================


# Run test cases
if __name__ == "__main__":
    test_data()
