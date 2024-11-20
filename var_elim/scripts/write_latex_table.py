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
    "n-elim-ub":          lambda item: str(int(item)).rjust(5) if not math.isnan(item) else "--".rjust(5),
    "n-elim-lb":          lambda item: str(int(item)).rjust(5) if not math.isnan(item) else "--".rjust(5),
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
    "elim-time":       lambda item: "%5.1f" % item if item >= 0.005 else "   --",
    "solve-time":      lambda item: "%5.1f" % item if item >= 0.05 else "%5.2f" % item,
    "build-time":      lambda item: "%5.1f" % item if item >= 0.05 else "%5.2f" % item,
    "init-time":       lambda item: "%5.1f" % item if item >= 0.05 else "%5.2f" % item,
    "function-time":   lambda item: "%5.1f" % item if item >= 0.05 else "%5.2f" % item,
    "jacobian-time":   lambda item: "%5.1f" % item if item >= 0.05 else "%5.2f" % item,
    "hessian-time":    lambda item: "%5.1f" % item if item >= 0.05 else "%5.2f" % item,
    "n-iter":          lambda item: str(int(item)).rjust(3),
    "ave-ls-trials":   lambda item: "%3.1f" % item,
    "function-per100": lambda item: "%5.2f" % item,
    "jacobian-per100": lambda item: "%5.2f" % item,
    "hessian-per100":  lambda item: "%5.2f" % item,
    "other-per100":    lambda item: "%5.2f" % item,

    r"Function (\%)":    lambda item: "%3i" % round(item),
    r"Jacobian (\%)":    lambda item: "%3i" % round(item),
    r"Hessian (\%)":     lambda item: "%3i" % round(item),
    r"Other (\%)":       lambda item: "%3i" % round(item),

    # At this point, we only use these keys for percent converged in the sweep summary
    # table. But if we end up using these in other tables, we should update the column
    # headers to something more specific.
    "distill":   lambda item: "%3i" % round(item),
    "mb-steady": lambda item: "%3i" % round(item),
    "pipeline":  lambda item: "%3i" % round(item),
    "total":     lambda item: "%3i" % round(item),
}


# Derived quantities that are more reader-friendly to display
CALCULATE = {
    "nnz/con":        lambda row: row["nnz"] / row["ncon"],
    "nnz-linear/con": lambda row: row["nnz-linear"] / row["ncon"],
    r"Function (\%)":   lambda row: 100 * row["function-per100"] / (row["function-per100"] + row["jacobian-per100"] + row["hessian-per100"] + row["other-per100"]),
    r"Jacobian (\%)":   lambda row: 100 * row["jacobian-per100"] / (row["function-per100"] + row["jacobian-per100"] + row["hessian-per100"] + row["other-per100"]),
    r"Hessian (\%)":    lambda row: 100 * row["hessian-per100"] / (row["function-per100"] + row["jacobian-per100"] + row["hessian-per100"] + row["other-per100"]),
    r"Other (\%)":      lambda row: 100 * row["other-per100"] / (row["function-per100"] + row["jacobian-per100"] + row["hessian-per100"] + row["other-per100"]),
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

    "build-time": r"$t_{\rm build}$",
    "elim-time": r"$t_{\rm elim}$",
    "init-time": r"$t_{\rm init}$",
    "solve-time": r"$t_{\rm solve}$",
    "n-iter": "Iter.",
    r"Function (\%)": r"Func.",
    r"Jacobian (\%)": r"Jac.",
    r"Hessian (\%)": r"Hess.",
    r"Other (\%)": r"Other",

    "distill": "DIST",
    "mb-steady": "MB",
    "opf": "OPF",
    "pipeline": "PIPE",

    "no-elim": "--",
    "d1": "LD1",
    "trivial": "ECD2",
    "linear-d2": "LD2",
    "d2": "D2",
    "ampl": "GR",
    "matching": "LM",

    "total": "Total",
}


def dataframe_to_latex(df, columns=None, methods=None):
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
        if methods is not None and row["method"] not in methods:
            # If methods were specified, skip any methods that were omitted.
            # This is for generating the matching-specific results table.
            continue
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

    if "model" in columns:
        model_idx = columns.index("model")
    else:
        model_idx = None
    if "method" in columns:
        method_idx = columns.index("method")
    else:
        method_idx = None

    last_model = None
    last_method = None
    for line in lines:
        # branching on methods is None is kind of a hack here. This might
        # not be what we want for all tables
        #if methods is None and (
        #    last_model is not None and line[0] != last_model
        #):
        #
        # Is this logic better?
        if (
            last_model is not None and model_idx is not None and line[model_idx] != last_model
            and last_method is not None and method_idx is not None and line[method_idx] != last_method
        ):
            # TODO: Multirow for each model?
            latex_lines.append("\\hline\n")

        if model_idx is not None:
            last_model = line[model_idx]
        if method_idx is not None:
            last_method = line[method_idx]

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
        "build-time",
        "elim-time",
        "init-time",
        "solve-time",
        #"function-time",
        #"jacobian-time",
        #"hessian-time",
        "n-iter",
        #"ave-ls-trials",
        #"function-per100",
        #"jacobian-per100",
        #"hessian-per100",
        #"other-per100",
        r"Function (\%)",
        r"Jacobian (\%)",
        r"Hessian (\%)",
        r"Other (\%)",
    ]
    return dataframe_to_latex(df, columns=columns)


def _generate_convergence_summary_table(df):
    columns = ["method", "distill", "mb-steady", "pipeline", "total"]
    return dataframe_to_latex(df, columns=columns)


def _generate_matching_table(df):
    columns = ["model", "method", "n-elim-lb", "n-elim", "n-elim-ub"]
    methods = ("matching",)
    return dataframe_to_latex(df, columns=columns, methods=methods)


def generate_table(df, which):
    dispatcher = {
        "structure": _generate_structure_table,
        "solvetime": _generate_solvetime_table,
        "matching-bounds": _generate_matching_table,
        "convergence-summary": _generate_convergence_summary_table,
    }
    return dispatcher[which](df)


def main(args):
    df = pd.read_csv(args.input_fpath)

    if args.which is not None:
        which = args.which
    elif "structure" in args.input_fpath and "solvetime" not in args.input_fpath:
        which = "structure"
    elif "solvetime" in args.input_fpath and "structure" not in args.input_fpath:
        which = "solvetime"
    elif "sweep-summary" in args.input_fpath:
        which = "convergence-summary"
    else:
        raise ValueError(
            "--which (table type) not provided and cannot be inferred from"
            " input file name"
        )

    print(f"Generating {which} table")
    table_str = generate_table(df, which)

    if not args.no_save:
        if args.fname is not None:
            output_fname = args.fname
        else:
            input_basename = os.path.basename(args.input_fpath)
            # exclude the .csv or .CSV extension
            input_noext = input_basename[:-4]
            if which == "matching-bounds":
                output_fname = input_noext + "-matching.txt"
            else:
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
