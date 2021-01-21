import inspect
import numpy as np

# This subclass of a NumPy matrix has relevant information to
# the process used when transforming a Data object into a
# purely numerical matrix of values for prediction purposes.
# 
# Numeric.to_real   -- the function for processing regular
#                      vectors into associated real vectors
# Numeric.from_real -- the function for processing real vectors
#                      back into regular vector
# Numeric.names     -- meaningful column names for the real vectors
# Numeric.nums      -- standard numerical indices
# Numeric.cats      -- translated categorical indices
# Numeric.data      -- The numeric version of the data.
# Numeric.shift     -- How to shift this data to have min 0.
# Numeric.scale     -- How to scale this data (after shift) to 
#                      have a domain of [0., 1.]                    
class Numeric(np.ndarray):
    print_column_width = None
    to_real = None
    from_real = None
    names = None
    nums = None
    cats = None
    shift = None
    scale = None

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

# Assign the comments to the official documentation spot.
Numeric.__doc__ = inspect.getcomments(Numeric).replace("# ","")
