import var_elim.scripts.config as config
import pandas as pd
import os
import math


FORMAT = {
    "model":  lambda item: item.ljust(9),
    "method": lambda item: item.ljust(9),

    "nvar":               lambda item: str(int(item)).rjust(6),
    "ncon":               lambda item: str(int(item)).rjust(6),
    "n-elim":             lambda item: str(int(item)).rjust(5),
    "n-elim-bound":       lambda item: str(int(item)).rjust(5) if not math.isnan(item) else "--".rjust(5),
    "nnz":                lambda item: str(int(item)).rjust(6),
    "nnz-linear":         lambda item: str(int(item)).rjust(6),
    "nnz-hessian":        lambda item: str(int(item)).rjust(6),
    "nnode-pyomo":        lambda item: str(int(item)).rjust(7),
    "nnode-nl-linear":    lambda item: str(int(item)).rjust(7),
    "nnode-nl-nonlinear": lambda item: str(int(item)).rjust(7),

    "nnz/con":        lambda item: "%3.2f" % item,
    "nnz-linear/con": lambda item: "%3.2f" % item,

    "success":         lambda item: str(item).ljust(5),
    "feasible":        lambda item: str(item).ljust(5),
    "elim-time":       lambda item: "%5.2f" % item,
    "solve-time":      lambda item: "%5.2f" % item,
    "function-time":   lambda item: "%5.2f" % item,
    "jacobian-time":   lambda item: "%5.2f" % item,
    "hessian-time":    lambda item: "%5.2f" % item,
    "n-iter":          lambda item: str(int(item)).rjust(3),
    "ave-ls-trials":   lambda item: "%3.1f" % item,
    "function-per100": lambda item: "%5.2f" % item,
    "jacobian-per100": lambda item: "%5.2f" % item,
    "hessian-per100":  lambda item: "%5.2f" % item,
    "other-per100":    lambda item: "%5.2f" % item,
}


# Derived quantities that are more reader-friendly to display
CALCULATE = {
    "nnz/con":        lambda row: row["nnz"] / row["ncon"],
    "nnz-linear/con": lambda row: row["nnz-linear"] / row["ncon"],
}


# Map from names used in the code (e.g. column headers) to names more suitable
# for the paper
RENAME = {
    # TODO: How to avoid hardcoding columns widths?
    "model": "Model".ljust(9),
    "method": "Method".ljust(9),
    "nvar": "Var.".ljust(6),
    "ncon": "Con.".ljust(6),
    "n-elim": "Elim.".ljust(5),
    "nnz/con": "NNZ/Con.",
    "nnz-linear/con": "Lin. NZ/Con.",
    "nnz-hessian": "Hess. NNZ",
    "nnode-nl-linear": "Lin. Nodes",
    "nnode-nl-nonlinear": "Nonlin. Nodes",

    "distill": "DIST",
    "mb-steady": "MB",
    "opf": "OPF",
    "pipeline": "PIPE",

    "no-elim": "--",
    "d1": "LD1",
    "trivial": "JP",
    "linear-d2": "LD2",
    "d2": "D2",
    "ampl": "GR",
    "matching": "LM",
}


def dataframe_to_latex(df, columns=None):
    if columns is None:
        # If not specified, we will use all columns
        columns = list(df.columns)
        # Pandas's convention for unnamed columns (including row indices
        # written from a dataframe) is to label them "Unnamed: N". We want
        # to ignore these.
        columns = [c for c in columns if not c.startswith("Unnamed")]

    lines = []
    for i, row in df.iterrows():
        line = []
        for c in columns:
            if c not in row:
                line.append(CALCULATE[c](row))
            elif c == "method" or c == "model":
                line.append(RENAME[row[c]])
            else:
                line.append(row[c])
        lines.append(line)
        #    [row[c] if c in row else CALCULATE[c](row) for c in columns]
        #)

    # TODO: Any formatting for column headers?
    paper_columns = [RENAME[c] if c in RENAME else c for c in columns]
    header_line = " & ".join(paper_columns) + "\\\\\n"
    latex_lines = [header_line, "\\hline\n"]
    last_model = None
    for line in lines:
        if last_model is not None and line[0] != last_model:
            # TODO: Multirow for each model?
            latex_lines.append("\\hline\n")
        last_model = line[0]

        formatted_line = [
            # Dispatch to a custom formatter for each column. If no formatter is
            # defined, we just convert to str
            FORMAT.get(col, lambda c: str(c))(item) for item, col in zip(line, columns)
        ]
        latex_lines.append(" & ".join(formatted_line) + "\\\\\n")
    latex_lines.append("\\hline\n")
    latex_table = "".join(latex_lines)
    return latex_table


def _generate_structure_table(df):
    columns = [
        "model",
        "method",
        "nvar",
        "ncon",
        "n-elim",
        #"n-elim-bound",
        "nnz/con",
        "nnz-linear/con",
        "nnz-hessian",
        #"nnode-pyomo",
        #"nnode-nl-linear",
        #"nnode-nl-nonlinear",
    ]
    return dataframe_to_latex(df, columns=columns)


def _generate_solvetime_table(df):
    columns = [
        "model",
        "method",
        "elim-time",
        "solve-time",
        #"function-time",
        #"jacobian-time",
        #"hessian-time",
        "n-iter",
        #"ave-ls-trials",
        "function-per100",
        "jacobian-per100",
        "hessian-per100",
        "other-per100",
    ]
    return dataframe_to_latex(df, columns=columns)


def generate_table(df, which):
    dispatcher = {
        "structure": _generate_structure_table,
        "solvetime": _generate_solvetime_table,
    }
    return dispatcher[which](df)


def main(args):
    df = pd.read_csv(args.input_fpath)

    if args.which is not None:
        which = args.which
    if "structure" in args.input_fpath and "solvetime" not in args.input_fpath:
        which = "structure"
    elif "solvetime" in args.input_fpath and "structure" not in args.input_fpath:
        which = "solvetime"
    else:
        raise ValueError(
            "--which (table type) not provided and cannot be inferred from"
            " input file name"
        )
    table_str = generate_table(df, which)

    if not args.no_save:
        if args.fname is not None:
            output_fname = args.fname
        else:
            input_basename = os.path.basename(args.input_fpath)
            # exclude the .csv or .CSV extension
            input_noext = input_basename[:-4]
            output_fname = input_noext + ".txt"

        output_fpath = os.path.join(args.results_dir, output_fname)
        print(f"Writing output to {output_fpath}")
        with open(output_fpath, "w") as f:
            f.write(table_str)

    print(table_str)


if __name__ == "__main__":
    argparser = config.get_argparser()
    argparser.add_argument(
        "input_fpath",
        help="CSV file containing results to write to a table",
    )
    argparser.add_argument("--fname", default=None, help="Basename of output file.")
    argparser.add_argument(
        "--which",
        default=None,
        help="Type of table. Options: 'structure', 'solvetime'",
    )
    args = argparser.parse_args()
    if args.model is not None:
        raise ValueError("--model argument is not supported")
    if not args.input_fpath.endswith(".csv") and not args.input_fpath.endswith(".CSV"):
        raise ValueError("Provided file must have extension .csv or .CSV")
    main(args)
