import var_elim.scripts.config as config


def main(args):
    model = config.CONSTRUCTOR_LOOKUP[args.model]()

    has_testproblem = args.model in config.TESTPROBLEM_LOOKUP
    if has_testproblem:
        testproblem = config.TESTPROBLEM_LOOKUP[args.model]
        params = testproblem.parameters
        ranges = [testproblem.parameter_ranges[p] for p in params]
        for p, (lb, ub) in zip(params, ranges):
            print(p)
            model_p = model.find_component(p)
            print(f"LB = {lb}")
            print(f"UB = {ub}")
            model_p.pprint()

    config.ELIM_LOOKUP[args.method](model)
    solver = config.get_optimization_solver()
    solver.solve(model, tee=True)


if __name__ == "__main__":
    argparser = config.get_argparser()
    args = argparser.parse_args()
    if args.model is None:
        raise ValueError("--model must be provided.")
    if args.method is None:
        raise ValueError("--method must be provided.")
    main(args)
