from util.math import transpose
from util.decorators import type_check

TABLE_SEPERATOR = "\\\\\n    "
# TABLE_TEXT%(<column_layout>, <num_columns>, <contents>)
TABLE_TEXT = """\\begin{table}
  \centering
  \\begin{tabular}{%s}
    \hline
    %s
    \hline
  \end{tabular}
  \caption{}
  \label{table:}
\end{table}"""
TABLE_MULTI_COL = "\multicolumn{%i}{%s}{%s}"
N_COLS = lambda n, side: side+"|".join(["c"]*n)+side


# Given a float 'value' and a format 'fmt' (with 'e' in it), generate
# a scientific notation string, but simplify numbers with exponents
# whose absolute value is less than or equal to 'min_exp'.
def simplify_scientific(value, fmt='.2e', min_exp=2):
    v, e = f"{value:{fmt}}".split("e")
    if "." in v: digits = len(v) - v.index(".") - 1
    else:        digits = v
    e = int(e)
    # If we are simplifying some sig figs, fix them.
    if (abs(e) <= min_exp):
        v = list(v)
        # (numerically stable) expansion of exponential form
        if   (e < 0): v = ['0',v.pop(1)] + ['0']*(-e-1) + v
        elif (e > 0): v = [v[0]]+v[2:2+e] + ['0']*(e-digits) + ['.'] + v[2+e:]
        v, e = ''.join(v + (['0'] if v[-1] == '.' else [])), 0
    string = f"{v} \\times 10^{{{e}}}"
    if (min_exp > 0): string = string.replace(" \\times 10^{0}", "")\
                                     .replace("10^{1}","10")
    return string

# Generate the text for a latex table block from a two dimensional
# python list. Emphasize numeric column values based on their sorted
# index and provided "wrap_" keyword arguments.
# 
#   "py_list"        -> A list of lists containing table rows.
# 
#   "header"         -> A row that will be placed above the table with
#                       an additional \hline separating it from contents.
# 
#   "fmt"            -> A format string for the numeric values.
# 
#   "wrap_nums"      -> A tuple of (before,after) strings to wrap
#                       around numbers, defaults to math mode dollar.
# 
#   "side_bars"      -> True if table should have outer-edge side bars.
# 
#   "simplify_sig"   -> The maximum absolute value exponent that will
#                       be written out fully for 'e' numeric formats.
# 
#   "decimal_align"  -> If provided, numbers will be converted into a
#                       2-column format aligned by the decimal rounded
#                       to "simplify_sig" nonzero digits.
# 
#   "wrap_i_j=(a,b)" -> Numerical elements of column 'i' with rank 'j'
#                       are converted to "a" + value + "b".
@type_check([list])
def latex_table(py_list, header=[], fmt=None, wrap_nums=("$","$"),
                side_bars=True, simplify_sig=-1, decimal_align=False, **kwargs):
    # Check row sizes to make sure they're all the same
    num_columns = max({len(l) for l in py_list})
    side = "|" if side_bars else ""
    table_rows = []
    # Get the numeric columns, translate keyword arguments into local
    # variables, and get the unique numeric column values sorted.
    is_numeric = lambda v: isinstance(v,int) or isinstance(v,float)
    to_wrap = [tuple(map(int,k[len('wrap_'):].split("_")))
               for k in kwargs if 'wrap_' == k[:len('wrap_')]]
    col_vals = [sorted(set(v for v in col if is_numeric(v))) for col in
                transpose(r for r in py_list if len(r) == num_columns)]
    doubled_columns = []
    # Create all of the rows of the table.
    for row in py_list:
        str_row = []
        # Fix the numeric values in "row" accordingly.
        for i in range(len(row)):
            str_row.append( str(row[i]) )
            # Skip if not formatting this column.
            if (not fmt) or (not is_numeric(row[i])):
                # Make the value take up 2 cells if it should.
                if i in doubled_columns:
                    str_row[-1] = TABLE_MULTI_COL%(2,"c",str_row[-1])
                continue
            # Get the format of this column based on usage.
            if type(fmt) == str: fmt_i = fmt
            else:                fmt_i = fmt[i]
            # Get the value and format.
            value = row[i]
            if ("e" in fmt_i): str_row[-1] = simplify_scientific(value, fmt_i, simplify_sig)
            else:              str_row[-1] = f"{value:{fmt_i}}"
            # If we have a full length row, wrap "indices" with text.
            before = []
            after = []
            if (len(row) == num_columns):
                for (col_idx, rank) in to_wrap:
                    if (i != range(num_columns)[col_idx]): continue
                    # Format the "col_vals" if they are not already formatted.
                    if is_numeric(col_vals[i][rank]):
                        col_vals[i] = [f"{v:{fmt_i}}" for v in col_vals[i]]
                        col_vals[i] = sorted(set(col_vals[i]), key=lambda v:col_vals[i].index(v))
                    # Check to see if this value matches the rank.
                    if (f"{value:{fmt_i}}" != col_vals[i][rank]): continue
                    # If it matches, then add the formatting.
                    start, end = kwargs[f'wrap_{col_idx}_{rank}']
                    before.insert(0,start)
                    after.append(end)
            # Wrap numbers with the number wrapping characters.
            if (len(wrap_nums) == 2) and all(type(w) == str for w in wrap_nums):
                start, end = wrap_nums
                # Add the numeric wrappers
                before.insert(0,start)
                after.append(end)
            elif (len(wrap_nums) == len(row)):
                start, end = wrap_nums[i]
                # Add the numeric wrappers
                before.insert(0,start)
                after.append(end)
            else:
                raise(Exception(f"'wrap_nums' must be 2-tuple or have same length as rows. Received '{wrap_nums}'."))
            # If this row is decimal aligned, then split it into two.
            if decimal_align:
                # Add this index to the set of doubled columns if it's not there.
                if i not in doubled_columns:
                    doubled_columns.append(i)
                # WARNING, this will not fix rows above it that are
                #          not properly doubled up.
                v = str_row.pop(-1)
                if ("e" in v):
                    # Convert a "pretty" (<base> \times 10^{<exp>}) into "<base>e<exp>".
                    if "times" in v:
                        n,_,e = v.strip().split()
                        e = e.split("{")[1].rstrip("}")
                        v = n+"e"+e
                    # Get the digits to the left and right of the decimal.
                    lr = ( v.split("e")[0] ).split(".")
                    if len(lr) == 1: left,right = lr[0], ""
                    else:            left,right = lr
                    # Undo the exponentiation by padding with 0's.
                    shift = int(v.split("e")[-1])
                    for i in range(abs(shift)):
                        if (shift > 0):
                            right = right + "0"
                            left  = left + right[0]
                            right = right[1:]
                        else:
                            left  = "0" + left
                            right = left[-1] + right
                            left  = left[:-1]
                else:
                    # Get the digits to the left and right of the decimal.
                    lr = v.split(".")
                    if len(lr) == 1: left,right = lr, ""
                    else:            left,right = lr
                # Strip off excessive 0's
                left = left.lstrip("0")
                right = right.rstrip("0")
                # Make sure there is some string there.
                if (left == ""):  left = "0"
                if (right == ""): right = "0"
                # Wrap both the left and right appropraitely.
                left_after = [c for c in "".join(after) if c in {"$","}"}]
                left = "".join(before+ [left] +left_after)
                right = "".join(before+ [right] +after)
                str_row += [left + "&" + right]
            else:
                str_row[-1] = "".join(before+ [str_row[-1]] +after)
        # Add the textual row to the table.
        table_rows.append( str_row )

    # Join all of the string table rows together with "&" (col delimeter).
    for i in range(len(table_rows)):
        table_rows[i] = " & ".join(list(map(str,table_rows[i])))
    
    table_rows[-1] += "\\\\"
    # Add the header if it was provided.
    if (len(header) > 0):
        if (len(header) < num_columns):
            row = TABLE_MULTI_COL%(num_columns, N_COLS(len(header),side), row)
        elif decimal_align and (len(doubled_columns) > 0):
            for i in doubled_columns:
                edge = "|" if ((i+1 == num_columns) or side) else ""
                header[i] = TABLE_MULTI_COL%(2,"c"+edge,str(header[i]))
            row = " & ".join(list(map(str,header)))
        else:
            row = " & ".join(list(map(str,header)))
        table_rows = [row] + table_rows
    # Add a horizontal line after the header.
    if len(header) > 0:
        table_rows[1] = "\hline"+TABLE_SEPERATOR[2:]+table_rows[1]
    # Set the column layout for the table.
    column_layout = N_COLS(num_columns,side)
    # Update the column layout for doubled-up columns that include decimals.
    if decimal_align:
        row = ["c"]*num_columns
        for i in doubled_columns: row[i] = "r@{.}l"
        column_layout = side+"|".join(row)+side
    # Return the latex table string
    return TABLE_TEXT%(column_layout, TABLE_SEPERATOR.join(table_rows))
