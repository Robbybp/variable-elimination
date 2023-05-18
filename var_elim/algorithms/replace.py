def define_variable_from_constraint(variable, constraint):
    """Get the expression that defines the variable according to the
    constraint.

    This only works if the variable participates linearly in the constraint.
    This should be used to generate substitution maps for
    ``replace_expressions``.

    Returns
    -------
    Numeric expression
        Defines variable using the expression in the constraint

    """
    # TODO:
    # - Convert constraint to standard repn
    # - Get linear coefficient of variable (if variable is nonlinear, throw
    #   error)
    # - Re-construct expression without varible's term
    pass
