from util.decorators import type_check

TABLE_SEPERATOR = "\\\\\n    "
# TABLE_TEXT%(<column_layout>, <num_columns>, <contents>)
TABLE_TEXT = """
\\begin{table}
  \centering
  \\begin{tabular}{%s}
    \hline
    %s
    \hline
  \end{tabular}
  \caption{}
  \label{table:%s}
\end{table}"""
TABLE_MULTI_COL = "\multicolumn{%i}{%s}{%s}"
N_COLS = lambda n: "|"+"|".join(["c"]*n)+"|"

# Generate the text for a latex table block from a two dimensional
# python list.
@type_check([list])
def latex_table(py_list, header=""):
    # Check row sizes to make sure they're all the same
    num_columns = max({len(l) for l in py_list})
    table_rows = [TABLE_MULTI_COL%(num_columns, N_COLS(1), header)]
    for row in py_list:
        table_rows.append( " & ".join(list(map(str,row))) )
    table_rows[-1] += "\\\\"
    if len(table_rows) >= 2:
        table_rows[1] = "\hline"+TABLE_SEPERATOR[2:]+table_rows[1]
    # Print out the latex table
    print(TABLE_TEXT%(N_COLS(num_columns), TABLE_SEPERATOR.join(table_rows), header))
    
