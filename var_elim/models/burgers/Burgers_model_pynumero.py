import pyomo.environ as pe
from pyomo import dae

def build_burgers_model(nfe_x=50, nfe_t=50, start_t=0, end_t=1, add_init_conditions=True):
    dt = (end_t - start_t) / float(nfe_t)

    start_x = 0
    end_x = 1
    dx = (end_x - start_x) / float(nfe_x)

    m = pe.Block(concrete=True)
    m.omega = pe.Param(initialize=0.02)
    m.v = pe.Param(initialize=0.01)
    m.r = pe.Param(initialize=0)

    m.x = dae.ContinuousSet(bounds=(start_x, end_x))
    m.t = dae.ContinuousSet(bounds=(start_t, end_t))

    m.y = pe.Var(m.x, m.t)
    m.dydt = dae.DerivativeVar(m.y, wrt=m.t)
    m.dydx = dae.DerivativeVar(m.y, wrt=m.x)
    m.dydx2 = dae.DerivativeVar(m.y, wrt=(m.x, m.x))

    m.u = pe.Var(m.x, m.t)

    def _y_init_rule(m, x):
        if x <= 0.5 * end_x:
            return 1
        return 0

    m.y0 = pe.Param(m.x, default=_y_init_rule)

    def _upper_x_bound(m, t):
        return m.y[end_x, t] == 0

    m.upper_x_bound = pe.Constraint(m.t, rule=_upper_x_bound)

    def _lower_x_bound(m, t):
        return m.y[start_x, t] == 0

    m.lower_x_bound = pe.Constraint(m.t, rule=_lower_x_bound)

    def _upper_x_ubound(m, t):
        return m.u[end_x, t] == 0

    m.upper_x_ubound = pe.Constraint(m.t, rule=_upper_x_ubound)

    def _lower_x_ubound(m, t):
        return m.u[start_x, t] == 0

    m.lower_x_ubound = pe.Constraint(m.t, rule=_lower_x_ubound)

    def _lower_t_bound(m, x):
        if x == start_x or x == end_x:
            return pe.Constraint.Skip
        return m.y[x, start_t] == m.y0[x]

    def _lower_t_ubound(m, x):
        if x == start_x or x == end_x:
            return pe.Constraint.Skip
        return m.u[x, start_t] == 0

    if add_init_conditions:
        m.lower_t_bound = pe.Constraint(m.x, rule=_lower_t_bound)
        m.lower_t_ubound = pe.Constraint(m.x, rule=_lower_t_ubound)

    # PDE
    def _pde(m, x, t):
        if t == start_t or x == end_x or x == start_x:
            e = pe.Constraint.Skip
        else:
            last_t = None
            # print(foo.last_t, t-dt, abs(foo.last_t - (t-dt)))
            # assert math.isclose(foo.last_t, t - dt, abs_tol=1e-6)
            e = m.dydt[x, t] - m.v * m.dydx2[x, t] + m.dydx[x, t] * m.y[x, t] == m.r + m.u[x, last_t]
        
        return e

    m.pde = pe.Constraint(m.x, m.t, rule=_pde)

    # Discretize Model
    disc = pe.TransformationFactory('dae.finite_difference')
    disc.apply_to(m, nfe=nfe_t, wrt=m.t, scheme='BACKWARD')
    disc.apply_to(m, nfe=nfe_x, wrt=m.x, scheme='CENTRAL')

    # Solve control problem using Pyomo.DAE Integrals
    def _intX(m, x, t):
        return (m.y[x, t] - m.y0[x]) ** 2 + m.omega * m.u[x, t] ** 2

    m.intX = dae.Integral(m.x, m.t, wrt=m.x, rule=_intX)

    def _intT(m, t):
        return m.intX[t]

    m.intT = dae.Integral(m.t, wrt=m.t, rule=_intT)

    def _obj(m):
        e = 0.5 * m.intT
        for x in sorted(m.x):
            if x == start_x or x == end_x:
                pass
            else:
                e += 0.5 * 0.5 * dx * dt * m.omega * m.u[x, start_t] ** 2
        return e

    m.obj = pe.Objective(rule=_obj)

    return m

def main():
    m = build_burgers_model(nfe_x = 100, nfe_t = 2000, end_t = 10)
    ipopt =pe.SolverFactory('ipopt')
    ipopt.solve(m, tee = True)
    
if __name__ == "__main__": 
    main()
