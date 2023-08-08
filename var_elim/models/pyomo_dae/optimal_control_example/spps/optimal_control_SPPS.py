#Code to replicate a SPPS optimal control problem
#Written by: Sakshi Naik and Megan Walsh 
#Course: Optimal control and reinforcement learning at Carnegie Mellon University

import pyomo.environ as pyo
import pyomo.dae as dae
import matplotlib.pyplot as plt

def make_model(nfe = 50, ncp = 3):
    m = pyo.ConcreteModel()
    #Fixed concentration of A1 in the injection streams
    CA1_in = 0.1
    CA2_in = 0.1
    V = 1
    
    
    #Parameters
    m.tf = pyo.Param(initialize = 100)
    m.K1 = pyo.Param(initialize = 10)
    m.K2 = pyo.Param(initialize = 1)
    m.K3 = pyo.Param(initialize = 10)
    m.K4 = pyo.Param(initialize = 1)
    m.K5 = pyo.Param(initialize = 1)
    m.K6 = pyo.Param(initialize = 1)
    m.K7 = pyo.Param(initialize = 1) 
    m.K8 = pyo.Param(initialize = 1)
    m.pen1 = pyo.Param(initialize = 1e-2, mutable = True)
    m.pen2 = pyo.Param(initialize = 1e-2, mutable = True)
    
    # pen1_list = [1e-3, 1e-2, 1e-1, 1]
    # pen2_list = [1e-3, 1e-2, 1e-1, 1]
    pen1_list = [0.001]
    pen2_list = [0.1]
    
    #Set
    m.t = dae.ContinuousSet(bounds = (0,m.tf))
    
    
    #Variables for concentrations
    m.S = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.A1 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.SA1 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.SA1A1 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.SA1A1A2 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.SA1A1A2A2 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.SA1A2 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.SA1A2A2 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.SA2 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.SA2A2 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.A2 = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    
    #Variables for flow
    m.FA1 = pyo.Var(m.t, bounds = (0,10))
    m.FA2 = pyo.Var(m.t, bounds = (0,10))
    m.Fout = pyo.Var(m.t, domain = pyo.NonNegativeReals)
    m.du1 = pyo.Var(m.t)
    m.du2 = pyo.Var(m.t)
    
    
    #Derivative Variables
    m.dS = dae.DerivativeVar(m.S, wrt = m.t)
    m.dA1 = dae.DerivativeVar(m.A1, wrt = m.t)
    m.dSA1 = dae.DerivativeVar(m.SA1, wrt = m.t)
    m.dSA1A1 = dae.DerivativeVar(m.SA1A1, wrt = m.t)
    m.dSA1A1A2 = dae.DerivativeVar(m.SA1A1A2, wrt = m.t)
    m.dSA1A1A2A2 = dae.DerivativeVar(m.SA1A1A2A2, wrt = m.t)
    m.dSA1A2 = dae.DerivativeVar(m.SA1A2, wrt = m.t)
    m.dSA1A2A2 = dae.DerivativeVar(m.SA1A2A2, wrt = m.t)
    m.dSA2 = dae.DerivativeVar(m.SA2, wrt = m.t)
    m.dSA2A2 = dae.DerivativeVar(m.SA2A2, wrt = m.t)
    m.dA2 = dae.DerivativeVar(m.A2, wrt = m.t)
    
    
    #Constraints
    def _cS(m, t):
        return m.dS[t] ==  -m.K1*m.S[t]*m.A1[t] - m.K5*m.S[t]*m.A2[t]
    m.cS = pyo.Constraint(m.t,rule = _cS)
    
    def _cA1(m, t):
        return m.dA1[t] ==  -m.K1*m.S[t]*m.A1[t] - m.K2*m.A1[t]*m.SA1[t]  +\
            m.FA1[t]/V*(CA1_in) - m.Fout[t]/V*m.A1[t]
    m.cA1 = pyo.Constraint(m.t,rule = _cA1)
    
    def _cSA1(m, t):
        return m.dSA1[t] ==  m.K1*m.S[t]*(m.A1[t]) - \
            m.K2*m.SA1[t]*(m.A1[t]) - m.K3*m.SA1[t]*m.A2[t]
    m.cSA1 = pyo.Constraint(m.t,rule = _cSA1)
    
    def _cSA1A1(m, t):
        return m.dSA1A1[t] ==  m.K2*m.SA1[t]*(m.A1[t])  - m.K7*m.SA1A1[t]*m.A2[t]
    m.cSA1A1 = pyo.Constraint(m.t,rule = _cSA1A1)
    
    def _cSA1A2(m, t):
        return m.dSA1A2[t] == m.K3*m.SA1[t]*m.A2[t] - m.K4*m.SA1A2[t]*m.A2[t]
    m.cSA1A2 = pyo.Constraint(m.t,rule = _cSA1A2)
    
    def _cSA1A2A2(m, t):
        return m.dSA1A2A2[t] == m.K4*m.SA1A2[t]*m.A2[t]
    m.cSA1A2A2 = pyo.Constraint(m.t,rule = _cSA1A2A2)
    
    def _cSA2(m, t):
        return m.dSA2[t] == m.K5*m.S[t]*m.A2[t] - m.K6*m.SA2[t]*m.A2[t]
    m.cSA2 = pyo.Constraint(m.t,rule = _cSA2)
    
    def _cSA2A2(m, t):
        return m.dSA2A2[t] ==  m.K6*m.SA2[t]*m.A2[t]
    m.cSA2A2 = pyo.Constraint(m.t,rule = _cSA2A2)
    
    def _cSA1A1A2(m, t):
        return m.dSA1A1A2[t] == m.K7*m.SA1A1[t]*m.A2[t] - m.K8*m.SA1A1A2[t]*m.A2[t]
    m.cSA1A1A2 = pyo.Constraint(m.t,rule = _cSA1A1A2)
    
    def _cSA1A1A2A2(m,t):
        return m.dSA1A1A2A2[t] == m.K8*m.SA1A1A2A2[t]
    m.cSA1A1A2A2 = pyo.Constraint(m.t,rule = _cSA1A1A2A2)
    
    def _cA2(m,t):
        return m.dA2[t] == -m.K3*m.SA1[t]*m.A2[t] - m.K4*m.SA1A2[t]*m.A2[t] -\
            m.K5*m.S[t]*m.A2[t] - m.K6*m.SA2[t]*m.A2[t] - m.K7*m.SA1A1[t]*m.A2[t]-\
            m.K8*m.SA1A1A2[t]*m.A2[t] + m.FA2[t]/V*(CA2_in) - m.Fout[t]/V*m.A2[t]     
    m.cA2 = pyo.Constraint(m.t,rule = _cA2)
    
    def _Flowrateout(m,t):
        return m.Fout[t] == m.FA1[t] + m.FA2[t] 
    m.Flow_bal = pyo.Constraint(m.t,rule = _Flowrateout)
    
    def _cdu1(m, ti):
        if ti == 0:
            return m.du1[ti] == 0
        return m.du1[ti] == m.FA1[ti] - m.FA1[m.t.prev(ti)]
    m.cdu1 = pyo.Constraint(m.t,rule = _cdu1)
    
    def _cdu2(m, ti):
        if ti == 0:
            return m.du2[ti] == 0
        return m.du2[ti] == m.FA2[ti] - m.FA2[m.t.prev(ti)]
    m.cdu2 = pyo.Constraint(m.t,rule = _cdu2)
    
    def _init_conditions(model):
        yield m.A1[0] == 0
        yield m.SA1[0] == 0
        yield m.SA1A1[0] == 0
        yield m.S[0] == 1
        yield m.SA1A1A2[0] == 0
        yield m.SA1A1A2A2[0] == 0
        yield m.SA1A2[0] == 0
        yield m.SA1A2A2[0] ==0
        yield m.SA2[0] == 0
        yield m.SA2A2[0] == 0
        yield m.A2[0] == 0
        yield m.FA1[0] == 0
        yield m.FA2[0] == 0
        
        
    m.init_conditions = pyo.ConstraintList(rule=_init_conditions)
    
    
    # Discretize model using Orthogonal Collocation
    discretizer = pyo.TransformationFactory('dae.collocation')
    
    discretizer.apply_to(m,nfe=nfe, ncp = ncp, scheme = 'LAGRANGE-RADAU')
    m = discretizer.reduce_collocation_points(m,var=m.FA1,
                                                  ncp=1,
                                                  contset=m.t)
    m = discretizer.reduce_collocation_points(m,var=m.FA2,
                                                  ncp=1,
                                                  contset=m.t)
    
    #objective 
    def _obj(model):
    	return -m.SA1A2[m.tf] +m.pen1*sum(m.du1[t]**2 for t in m.t) +m.pen2*sum(m.du2[t]**2 for t in m.t) 
    m.obj = pyo.Objective(rule=_obj)
    
    solver=pyo.SolverFactory('ipopt')

    return m

m = make_model(nfe = 2, ncp = 1)
solver=pyo.SolverFactory('ipopt')
solver.solve(m, tee= True)

# plotter(131, m.t ,m.SA1, m.SA1A2, title='Differential Variables') 
# plotter(132, m.t, m.FA1, title='Control Variable', points='o') 
# plotter(133, m.t, m.FA2, title='Control Variable', points='o') 
# plt.show()

