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
    "nnode-pyomo":        lambda item: str(int(item)).rjust(7),
    "nnode-nl-linear":    lambda item: str(int(item)).rjust(7),
    "nnode-nl-nonlinear": lambda item: str(int(item)).rjust(7),

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
        lines.append([row[c] for c in columns])

    # TODO: Any formatting for column headers?
    header_line = " & ".join(columns) + "\\\\\n"
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
    return dataframe_to_latex(df)


def _generate_solvetime_table(df):
    return dataframe_to_latex(df)


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
    if "structure" in args.input_fpath:
        which = "structure"
    elif "solvetime" in args.input_fpath:
        which = "solvetime"
    else:
        raise ValueError(
            "--which (table type) not provided and cannot be inferred from"
            " input file name"
        )
    table_str = generate_table(df, which)

    if args.fname is not None:
        output_fpath = os.path.join(args.results_dir, args.fname)
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
    argparser.add_argument(
        "--fname",
        default=None,
        help="Basename of output file. Default (None) is to not write a file.",
    )
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
