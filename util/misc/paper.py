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
    digits = len(v) - v.index(".") - 1
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
#   "wrap_i_j=(a,b)" -> Numerical elements of column 'i' with rank 'j'
#                       are converted to "a" + value + "b".
# 
@type_check([list])
def latex_table(py_list, header=[], fmt=None, wrap_nums=("$","$"),
                side_bars=True, simplify_sig=2, **kwargs):
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
    # Add the header if it was provided.
    if (len(header) > 0):
        row = " & ".join(list(map(str,header)))
        if (len(header) < num_columns):
            row = TABLE_MULTI_COL%(num_columns, N_COLS(len(header),side), row)
        table_rows += [row]
    # Create all of the rows of the table.
    for row in py_list:
        # Fix the numeric values in "row" accordingly.
        for i in range(len(row)):
            # Skip if not formatting this column.
            if (not fmt) or (not is_numeric(row[i])): continue
            # Get the format of this column based on usage.
            if type(fmt) == str: fmt_i = fmt
            else:                fmt_i = fmt[i]
            # Get the value and format.
            value = row[i]
            if ("e" in fmt_i): row[i] = simplify_scientific(value, fmt_i, simplify_sig)
            else:              row[i] = f"{value:{fmt_i}}"
            # If we have a full length row, wrap "indices" with text.
            if (len(row) == num_columns):
                for (col_idx, rank) in to_wrap:
                    if (i != range(num_columns)[col_idx]): continue
                    if (value != col_vals[i][rank]): continue
                    start, end = kwargs[f'wrap_{col_idx}_{rank}']
                    row[i] = start + row[i] + end
            # Wrap numbers with the number wrapping characters.
            if (len(wrap_nums) == 2) and all(type(w) == str for w in wrap_nums):
                start, end = wrap_nums
            elif (len(wrap_nums) == len(row)):
                start, end = wrap_nums[i]
            else:
                raise(Exception("'wrap_nums' must be 2-tuple or have same length as rows."))
            row[i] = start + row[i] + end
        # Add the textual row to the table.
        table_rows.append( " & ".join(list(map(str,row))) )
    table_rows[-1] += "\\\\"
    # Add a horizontal line after the header.
    if len(header) > 0:
        table_rows[1] = "\hline"+TABLE_SEPERATOR[2:]+table_rows[1]
    # Return the latex table string
    return TABLE_TEXT%(N_COLS(num_columns,side), TABLE_SEPERATOR.join(table_rows))
