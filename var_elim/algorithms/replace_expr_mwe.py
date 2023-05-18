import pyomo.environ as pyo
from pyomo.core.expr.visitor import replace_expressions


def main():
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3])
    expr = m.x[1] + 2*m.x[2] + 3*m.x[3]

    # replace_expressions accepts a map from expression (including just a
    # variable) ids to the "target" expression that the original should
    # be replaced with.
    substitution_map = {
        id(m.x[2]): m.x[1]**2 + 1/m.x[3]
    }

    new_expr = replace_expressions(expr, substitution_map)

    print("Before replacement") 
    print("------------------") 
    print(expr)

    print("After replacement") 
    print("-----------------") 
    print(new_expr)


if __name__ == "__main__":
    main()
