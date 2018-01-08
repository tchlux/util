import numpy as np
from util import COMMON_SEPERATORS, UPDATE_RATE, \
    MAX_ERROR_PRINTOUT, NP_TYPES, PY_TYPES
from util.classes import AtomicOpen

class NoFix(Exception): pass
class FailedFix(Exception): pass

# Return the type of a string (try converting to int and float)
def get_type(string):
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
def detect_seperator(filename, verbose=False):
    lines = []
    not_viable = set()
    with AtomicOpen(filename) as f:
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
        sep = max(seperators, key=lambda k: seperators[k])
        print("WARNING: Using '%s' as a seperator, even though it is unexpected."%(sep))
    # Return the automatically detected seperator
    return sep

# Pre:  "filename" is the name of a txt file that has a header
#       "sep" is the seperator used in the data file
#       "display" is True if the user wants the detected types printed
# Post: The header of the file as a list of strings,
#       the types of each header as a list,
#       the raw data as a [list of [lists of rows]]
def read_data_and_type(filename, sep=None, update_rate=UPDATE_RATE):
    if sep == None: sep = detect_seperator(filename, update_rate!=0)
    # Open the file to get the data
    with AtomicOpen(filename) as f:
        # Get the first line, assuming it's the header
        header = f.readline().strip().split(sep)
        # Get the first line of data and learn the types of each column
        line = f.readline().strip().split(sep)
        types = list(map(get_type, line))
        # Add the first line to the data store
        row_1 = tuple(t(v) for t,v in zip(types,line))
        # Get the rest of the raw data from the file and close it
        raw_data = f.readlines()
    # Print out a nicely formatted description of the detected types
    if update_rate != 0:
        print()
        print("Automatically detected types:")
        first_line  = ""
        second_line = ""
        for h,t in zip(header,types):
            first_line += "  "
            second_line += "  "
            length = max(len(h), len(t.__name__))
            first_line += h + " "*(length-len(h))
            second_line += t.__name__ + " "*(length-len(t.__name__))
        print(first_line)
        print(second_line)
    # Now process the rest of the data, dynamically overwriting old data
    to_remove = []
    for i,line in enumerate(raw_data):
        if (update_rate != 0) and not (i % update_rate):
            print("\r%0.1f%% complete"% (100.0*i/len(raw_data)), 
                  flush=True, end="")

        try:
            # Holder for the cast row of data
            raw_data[i] = tuple()
            # Cycle through the row, casting the types
            for c, (cast, value) in enumerate(zip(types, line.strip().split(sep))):
                try:
                    cast_value = cast(value)
                except:
                    # If the number was thought to be an int, try float instead
                    if (cast == int):
                        try:
                            cast_value = float(value)
                            types[c] = float
                        # If trying float failed, then there's a problem
                        except ValueError:
                            raise(FailedFix())
                    # If it's not an int to begin with, there's a problem
                    else:
                        raise(NoFix())
                raw_data[i] += (cast_value,)
        except (ValueError, FailedFix, NoFix):
            to_remove.append(i)
            if update_rate != 0:
                if len(to_remove) > MAX_ERROR_PRINTOUT:
                    if len(to_remove) == MAX_ERROR_PRINTOUT+1:
                        print("Suppressing output for errors... See all "+
                              "erroneous lines in seperate printout.")
                else:
                    print("\nERROR ON LINE %i, REMOVING FROM DATA:\n  %s\n"%(i+3,line))
    if update_rate != 0: print()
    # Remove lines that caused trouble when casting types
    for i in to_remove[::-1]:
        raw_data.pop(i)
    if (len(to_remove) > 0):
        print("WARNING: Encountered potentially erroneous lines while\n"+
              "         reading '%s'.\n"%(filename)+
              "         Consider setting 'verbose=True' and checking errors.")
    if (update_rate != 0) and (len(to_remove) > MAX_ERROR_PRINTOUT):
        print("ERRONEOUS LINES FROM DATA:\n  %s"%to_remove)
    # Return the header, types, and raw data
    return header, types, [row_1] + raw_data

# Given a file name, automatically detect the seperator and the types
# of each column, return all of this in a numpy structured array.
def read_struct(file_name, sep=None, verbose=False):
    print_rate = 0 if (not verbose) else UPDATE_RATE
    if verbose:
        print("Opening data file and reading into Python List...")
    try:
        header, types, data = read_data_and_type(file_name, sep, print_rate)
        # Sort the data by column precedent (standard python sort)
        # data.sort()
    except TypeError:
        msg = "Current version cannot read data with missing values\n"+\
              "            or values that change type throughout column. Types\n"+\
              "            are automatically detected based on first data row."
        raise(Exception("DATA ERROR: "+msg))
    # Convert the data into numpy form
    if verbose: print("Converting the data into numpy format...")
    dtypes = []
    for i,(h,t) in enumerate(zip(header,types)):
        if (t != str):
            pair = (h,) + NP_TYPES[t]
        else:
            max_len = max(len(row[i]) for row in data)
            pair = (h,) + NP_TYPES[t][:1]+(max_len,)
        dtypes.append(pair)
    dtypes = np.dtype(dtypes)
    if verbose: print("Data types:",dtypes)
    data = np.array(data, dtype=dtypes)
    return data


