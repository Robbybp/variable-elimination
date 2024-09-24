import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent
import pickle
import os
def build_qcp(data, init_custom = False, init_custom_param = None):
    #Load kUSD as units in the pint_registry if the unit is absent
    try:
        getattr(pyunits.pint_registry, 'kUSD')
    except:
        pyunits.load_definitions_from_strings(['kUSD = [currency]'])
    
    model = pyo.ConcreteModel()

    # -------------------SETS-----------------------------
    # Nodes
    model.s_N = pyo.Set(initialize=data['s_N'], doc='Nodes')
    model.s_NMS = pyo.Set(initialize=data['s_NMS'], doc='Mixers/Splitters')
    model.s_NP = pyo.Set(initialize=data['s_NP'], doc='Production pads')
    model.s_NC = pyo.Set(initialize=data['s_NC'], doc='Completion pads')
    model.s_NS = pyo.Set(initialize=data['s_NS'], doc='Storage')
    model.s_ND = pyo.Set(initialize=data['s_ND'], doc='Disposal')
    model.s_NW = pyo.Set(initialize=data['s_NW'], doc='Fresh water source')
    model.s_NTIN = pyo.Set(initialize = data['s_NTIN'], doc = 'Inlet to Treatment Facility')
    model.s_NTTW = pyo.Set(initialize = data['s_NTTW'], doc = 'Treatmed Water of Treatment Facility')
    model.s_NTCW = pyo.Set(initialize = data['s_NTCW'], doc = 'Concentrated Water of Treatment Facility') 
    
    # Arcs
    model.s_A = pyo.Set(initialize = data['s_A'], doc = "Arcs")
    model.s_Ain = pyo.Param(model.s_N, initialize = data['s_Ain'], doc = 'Arcs going into node')
    model.s_Aout = pyo.Param(model.s_N, initialize = data['s_Aout'], doc = "Arcs coming out of node")

    # Components
    model.s_Q = pyo.Set(initialize = data['s_Q'], doc = "Components")
    
    # TODO: Can't do this because then you will have to separate sets for each node
    model.s_Qalpha = pyo.Param(model.s_NTIN, initialize = data['s_Qalpha'], doc = "Components being treated at Treated Water Node")

    # Time Periods
    model.s_T = pyo.Set(initialize = data['s_T'], doc = "Time Periods")

    # ----------------------PARAMETERS-------------------------------------
    # Generation and Consumption Parameters
    model.p_FGen = pyo.Param(model.s_NMS|model.s_NP|model.s_NC, model.s_T, initialize = data['p_FGen'], doc = "Flow Generated", units = pyunits.bbl/pyunits.day)
    model.p_CGen = pyo.Param(model.s_NMS|model.s_NP|model.s_NC, model.s_Q, model.s_T, initialize = data['p_CGen'], doc = "Concentration in flow generated", units = pyunits.g/pyunits.liter)
    model.p_SGen = pyo.Param(model.s_NMS|model.s_NP|model.s_NC, model.s_Q, model.s_T, initialize = data['p_SGen'], doc = "Solids generated", units = pyunits.bbl/pyunits.day*pyunits.g/pyunits.liter)
    model.p_FCons = pyo.Param(model.s_NMS|model.s_NP|model.s_NC, model.s_T, initialize = data['p_FCons'], doc = "Flow consumed", units = pyunits.bbl/pyunits.day)

    # Inventory Parameters
    model.p_I0 = pyo.Param(model.s_NS, initialize = data['p_I0'], doc = "Initial Inventory", units = pyunits.bbl)
    model.p_IS0 = pyo.Param(model.s_NS, model.s_Q, initialize = data['p_IS0'], doc = "Initial solids in inventory", units = pyunits.bbl*pyunits.g/pyunits.liter)
    
    # Treatment Efficiencies
    #Should be set to 1 in the case where we use MVC as the treatment unit. 
    model.p_alpha = pyo.Param(model.s_NTIN, model.s_Q, initialize = data['p_alpha'], doc = "Treatment Efficiency")

    # Capacity limits for sinks
    model.p_Cap = pyo.Param(model.s_ND | model.s_NW | model.s_NTIN, initialize = data['p_Cap'], doc = "Capacity of sources and sinks", units = pyunits.bbl/pyunits.day)
    model.p_Cap_treat_min = pyo.Param(model.s_NTIN, initialize = data['p_Cap_treat_min'], doc = "Minimum flow to the treatment node", units = pyunits.bbl/pyunits.day)
    # Time discretization
    model.p_dt = pyo.Param(initialize = data['p_dt'], doc = "Time discretized for inventory", units = pyunits.day)

    # Minimum treatment concentration required
    model.p_Cmin = pyo.Param(model.s_NTCW, model.s_Q, initialize = data['p_Cmin'], mutable = True, doc = "Minimum concentration required at the concentrated water node", units = pyunits.g/pyunits.liter)
    model.p_Cmax = pyo.Param(initialize = 300, doc = "Maximum concentration in the network", units = pyunits.g/pyunits.liter)
    # Cost Parameters
    model.p_betaArc = pyo.Param(model.s_A, initialize = data['p_betaArc'], doc = "Operational costs of arcs" ,units = pyunits.kUSD/pyunits.bbl)
    model.p_betaD = pyo.Param(model.s_ND, initialize = data['p_betaD'], mutable = True, doc = "Operational cost of disposal site",units = pyunits.kUSD/pyunits.bbl)
    model.p_betaW = pyo.Param(model.s_NW, initialize = data['p_betaW'], doc = "Costs of sourcing freshwater",units = pyunits.kUSD/pyunits.bbl)
    model.p_betaT = pyo.Param(model.s_NTIN, initialize = data['p_betaT'], mutable = True, doc = "Operating cost of treatment facility",units = pyunits.kUSD/pyunits.bbl)
    model.p_betaS = pyo.Param(model.s_NS, initialize = data['p_betaS'], mutable = True, doc = "Costs of storing produced water",units = pyunits.kUSD/pyunits.bbl)
    
    #Differentiating the costs to send water to inventory slightly to avoid degeneracy
    data_inv_cost = {}
    counter = 0
    for n in model.s_NS:
        delta = 0.00001
        for t in model.s_T:
            data_inv_cost[n, t] = pyo.value(model.p_betaS[n]+ len(model.s_T)*delta - delta*counter)
            counter = counter + 1
            
   
    model.p_betaSt = pyo.Param(model.s_NS, model.s_T, initialize = data_inv_cost, mutable = True, doc = "Costs of storing produced water",units = pyunits.kUSD/pyunits.bbl)
    
    model.p_gammaS = pyo.Param(model.s_NS, initialize = data['p_gammaS'], mutable = True, doc = "Earnings from retrieving water from storage tanks",units = pyunits.kUSD/pyunits.bbl)
    model.p_gammaT = pyo.Param(model.s_NTTW, initialize = data['p_gammaT'], mutable = True, doc = "Earnings from treating water",units = pyunits.kUSD/pyunits.bbl)
    
    # Functions
    model.p_nodeUp = pyo.Param(model.s_A, initialize = data['p_nodeUp'], doc = "Upstream node of the arc")
    model.p_nodeDown = pyo.Param(model.s_A, initialize = data['p_nodeDown'], doc = "Downstream node of the arc")
    model.p_tp = pyo.Param(model.s_A, initialize = data['p_tp'], doc = "Type of arc")
    model.p_treatedIN = pyo.Param(model.s_NTTW | model.s_NTCW, initialize = data['p_treatedIN'], doc = "Inlet node for each treatment facility")
    
    # --------------------VARIABLES-----------------------------
    model.v_F = pyo.Var(model.s_A, model.s_T, bounds = data['p_Fbounds'], doc = "Flow in arc", initialize = 5, units = pyunits.liter/pyunits.s)
    model.v_S = pyo.Var(model.s_A, model.s_Q, model.s_T, bounds = (0,30000),doc = "Solid flow in arc", initialize = 30, units = pyunits.g/pyunits.s)
    model.v_I = pyo.Var(model.s_NS, model.s_T, bounds=data['p_Ibounds'], doc = "Inventory in node", initialize = 10, units = pyunits.liter)
    model.v_IS = pyo.Var(model.s_NS, model.s_Q, model.s_T, bounds = (0, 30000), doc ="Solids in the inventory node", initialize = 100, units = pyunits.g)
    model.v_C = pyo.Var(model.s_NTIN |model.s_NTTW | model.s_NTCW | model.s_NS , model.s_Q, model.s_T, bounds = (0,300), doc = "Concentration of water at nodes", initialize = 100, units = pyunits.g/pyunits.liter)
    model.v_alphaW = pyo.Var(model.s_NTIN, model.s_T, bounds = (0.1,0.9), doc = "Treatment Efficiency", initialize = 0.5)
    model.v_Ctreatin = pyo.Var(model.s_NTIN, model.s_T, doc = "Concentration of salts going to treatment", bounds = (0, 210), initialize = 100, units = pyunits.g/pyunits.liter)
    #Custom initialization
    if init_custom:
        model.v_I[:, :] = init_custom_param['I_init']
        model.v_IS[:, :, :] = init_custom_param['IS_init']
   
    #Node R01 is just pretreatment. So fix it's alpha 
    #Scaling factor 
    scale_fac = 1/100
    # ----------------------OBJECTIVE---------------------------
    
    @model.Expression()
    def arc_cost(m):
        return sum(sum(pyunits.convert(m.p_betaArc[a], pyunits.kUSD/pyunits.liter)*m.v_F[a,t] for a in m.s_A) for t in m.s_T)*pyunits.convert(m.p_dt, pyunits.s)
    
    @model.Expression()
    def disp_cost(m):
        #This disposal cost function was used for the PWS meeting results
        #To avoid disposal cost when we are disposing brine and fresh water 
        #return sum(sum(pyunits.convert(m.p_betaD[n], pyunits.kUSD/pyunits.liter)*sum(m.v_F[a,t] for a in m.s_Ain[n]) for n in m.s_ND if n not in ['D02', 'D03']) for t in m.s_T)*pyunits.convert(m.p_dt, pyunits.s)
        return sum(sum(pyunits.convert(m.p_betaD[n], pyunits.kUSD/pyunits.liter)*sum(m.v_F[a,t] for a in m.s_Ain[n]) for n in m.s_ND) for t in m.s_T)*pyunits.convert(m.p_dt, pyunits.s)
    
    @model.Expression()
    def fresh_cost(m):
        return sum(sum(pyunits.convert(m.p_betaW[n], pyunits.kUSD/pyunits.liter)*sum(m.v_F[a,t] for a in m.s_Aout[n]) for n in m.s_NW) for t in m.s_T)*pyunits.convert(m.p_dt, pyunits.s)
    
    @model.Expression()
    def stor_cost(m):
        return sum(sum(pyunits.convert(m.p_betaSt[n,t], pyunits.kUSD/pyunits.liter)*sum(m.v_F[a,t] for a in m.s_Ain[n]) for n in m.s_NS) for t in m.s_T)*pyunits.convert(m.p_dt, pyunits.s)
    
    @model.Expression()
    def stor_rev(m):
        return sum(sum(pyunits.convert(m.p_gammaS[n], pyunits.kUSD/pyunits.liter)*sum(m.v_F[a,t] for a in m.s_Aout[n]) for n in m.s_NS) for t in m.s_T)*pyunits.convert(m.p_dt, pyunits.s)
    
    @model.Expression()
    def treatment_rev(m):
        return sum(sum(pyunits.convert(m.p_gammaT[n], pyunits.kUSD/pyunits.liter)*sum(m.v_F[a,t] for a in m.s_Aout[n]) for n in m.s_NTTW) for t in m.s_T)*pyunits.convert(m.p_dt, pyunits.s)
    
    # #Constraint to avoid degeneracy in the inventory variable
    # @model.Expression()
    # def inv_penalty(m):
    #     return 1e-3*sum(sum((m.v_I[n, t] - m.v_I[n, m.s_T.prev(t)])**2 for n in m.s_NS)for t in m.s_T if t != m.s_T.first())
    # @model.Expression()
    # def treat_cost(m):
    #     return sum(sum(pyunits.convert(m.p_betaT[n], pyunits.kUSD/pyunits.liter)*sum(m.v_F[a,t] for a in m.s_Ain[n]) for n in m.s_NTIN) for t in m.s_T)*pyunits.convert(m.p_dt, pyunits.s)
    
    @model.Objective(doc = "Cost Minimization Objective")
    def obj(m):
        return (m.arc_cost + m.disp_cost + m.fresh_cost + m.stor_cost - m.stor_rev - m.treatment_rev)

    # ---------------------------CONSTRAINTS-----------------------------------    
    # flow for non-inventory terms
    @model.Constraint(model.s_NP | model.s_NMS, model.s_T, doc = "General flow equation for non inventory nodes")
    def noinvflow(m,n,t):
        return sum(m.v_F[a,t] for a in m.s_Ain[n]) + pyunits.convert(m.p_FGen[n,t], pyunits.liter/pyunits.sec)== sum(m.v_F[a,t] for a in m.s_Aout[n]) + pyunits.convert(m.p_FCons[n,t], pyunits.liter/pyunits.s)
    assert_units_consistent(model.noinvflow)
    
    # solids flow for non-inventory terms 
    @model.Constraint(model.s_NP | model.s_NMS, model.s_Q, model.s_T, doc = "General solids flow equation for non inventory nodes")
    def noinvsolidsflow(m,n,q,t):
        return scale_fac*(sum(m.v_S[a,q,t] for a in m.s_Ain[n]) + pyunits.convert(m.p_SGen[n,q,t], pyunits.g/pyunits.s)) == scale_fac*(sum(m.v_S[a,q,t] for a in m.s_Aout[n]))
    assert_units_consistent(model.noinvsolidsflow)
    
    # Splitter and Mixers outlet concentration constraints
    model.splitter_conc = pyo.ConstraintList()
    for n in model.s_NMS:
        if len(model.s_Aout[n]) > 1:
            for q in model.s_Q:
                for t in model.s_T:
                    for a in model.s_Aout[n]:
                        model.splitter_conc.add(scale_fac*(model.v_S[a, q, t]*sum(model.v_F[a, t] for a in model.s_Ain[n])) == scale_fac*(sum(model.v_S[a, q, t] for a in model.s_Ain[n])*model.v_F[a, t]))
    assert_units_consistent(model.splitter_conc)     
        
    # Completion Pads
    # flow and concentration constraints are done together below
    model.Cflow = pyo.ConstraintList(doc = "Completions pad flow")
    model.Csolids = pyo.ConstraintList(doc = "Completions pad solids flow")
    for n in model.s_NC:
        for t in model.s_T:
            #condition where completion pad is producing 
            if pyo.value(model.p_FGen[n,t])>0 and pyo.value(model.p_FCons[n,t]) == 0:
                model.Cflow.add(pyunits.convert(model.p_FGen[n,t], pyunits.liter/pyunits.s) == sum(model.v_F[a,t] for a in model.s_Aout[n]))  # flow constraint
                model.Cflow.add(0 == sum(model.v_F[a,t] for a in model.s_Ain[n]))  # flow constraint
                for q in model.s_Q:
                    #Redundant constraint since the inequality bound takes care of this
                    model.Csolids.add(sum(model.v_S[a, q, t] for a in model.s_Ain[n]) == 0)
                    model.Csolids.add(pyunits.convert(model.p_SGen[n,q,t], pyunits.g/pyunits.s) == sum(model.v_S[a, q, t] for a in model.s_Aout[n]))    # solid flow constraint
                    
                    
            #condition where completion pad is consuming
            elif pyo.value(model.p_FGen[n,t]) == 0 and pyo.value(model.p_FCons[n,t]) > 0:
                model.Cflow.add(pyunits.convert(model.p_FCons[n,t], pyunits.liter/pyunits.s) == sum(model.v_F[a,t] for a in model.s_Ain[n]))  # flow constraint
                model.Cflow.add(0 == sum(model.v_F[a,t] for a in model.s_Aout[n]))  # flow constraint
                # Solids are not tracked in this case
                #The sum of solids out of the completions pad is 0. Since solids have a LB of 0
                #All of them have to be 0
                #Redundant constraint since the inequality bound takes care of this
                for q in model.s_Q:
                    model.Csolids.add(sum(model.v_S[a,q,t] for a in model.s_Aout[n]) == 0)
            
            #condition where completion pad is neither producing nor consuming
            elif pyo.value(model.p_FGen[n,t]) == 0 and pyo.value(model.p_FCons[n,t]) == 0:
                model.Cflow.add(0 == sum(model.v_F[a,t] for a in model.s_Aout[n]))  # flow constraint
                model.Cflow.add(0 == sum(model.v_F[a,t] for a in model.s_Ain[n]))  # flow constraint
                for q in model.s_Q:
                    model.Csolids.add(sum(model.v_S[a, q, t] for a in model.s_Aout[n]) == 0)  # solid flow constraint
                    model.Csolids.add(sum(model.v_S[a, q, t] for a in model.s_Ain[n]) == 0)  # solid flow constraint
            
            #condition where completion pad is both producing and consuming
            elif pyo.value(model.p_FGen[n,t]) > 0 and pyo.value(model.p_FCons[n,t]) > 0:
                model.Cflow.add(pyunits.convert(model.p_FGen[n,t], pyunits.liter/pyunits.s) == sum(model.v_F[a,t] for a in model.s_Aout[n]))  # flow constraint
                model.Cflow.add(model.p_FCons[n,t] == sum(model.v_F[a,t] for a in model.s_Ain[n]))  # flow constraint
                for q in model.s_Q:
                    model.Cconc.add(pyunits.convert(model.p_SGen[n,q,t], pyunits.g/pyunits.s)== sum(model.S[a,q,t] for a in model.s_Aout[n]))    # solids flow constraint
    
    assert_units_consistent(model.Cflow)     
    assert_units_consistent(model.Csolids)     
        
    # Storage
    @model.Constraint(model.s_NS, model.s_T, doc = "Storage inventory balance")
    def Sinv(m,n,t):
        if t == m.s_T.first():
            return scale_fac*(m.v_I[n,t])== scale_fac*(pyunits.convert(m.p_I0[n], pyunits.liter) + sum(m.v_F[a,t] for a in m.s_Ain[n])*pyunits.convert(m.p_dt, pyunits.s) - sum(m.v_F[a,t] for a in m.s_Aout[n])*pyunits.convert(m.p_dt, pyunits.s))
        else:
            return scale_fac*(m.v_I[n,t]) == scale_fac*(m.v_I[n,m.s_T.prev(t)] + sum(m.v_F[a,t] for a in m.s_Ain[n])*pyunits.convert(m.p_dt, pyunits.s) - sum(m.v_F[a,t] for a in m.s_Aout[n])*pyunits.convert(m.p_dt, pyunits.s))
    assert_units_consistent(model.Sinv)
    
    @model.Constraint(model.s_NS, model.s_Q, model.s_T, doc = "Storage solids balance")
    def Sconc(m,n,q,t):
        if t == m.s_T.first():
            return scale_fac*(m.v_IS[n, q, t]) == scale_fac*(pyunits.convert(m.p_IS0[n, q], pyunits.g)+ 
                    sum(m.v_S[a,q,t] for a in m.s_Ain[n])*pyunits.convert(m.p_dt, pyunits.s) - sum(m.v_S[a,q,t] for a in m.s_Aout[n])*pyunits.convert(m.p_dt, pyunits.s))
        else:
            return scale_fac*(m.v_IS[n, q, t]) == scale_fac*(m.v_IS[n, q, m.s_T.prev(t)] + 
                    sum(m.v_S[a, q, t] for a in m.s_Ain[n])*pyunits.convert(m.p_dt, pyunits.s) - sum(m.v_S[a, q, t] for a in m.s_Aout[n])*pyunits.convert(m.p_dt, pyunits.s))
    assert_units_consistent(model.Sconc)
    
    #We need to maintain a concentration variable in the inventory unit
    #Since it is a storage tank and the flows aren't just in and out for a simple flow balance
    @model.Constraint(model.s_NS, model.s_Q, model.s_T)
    def Cinvcal(m,n,q,t):
        return m.v_C[n,q,t]*m.v_I[n,t] == m.v_IS[n, q, t]
    assert_units_consistent(model.Cinvcal)
    
    # #This constraint is essential if there is a pretreatment unit before storage. 
    # #Helps a little with situations where flow is 0 to storage and concentration can take any value since 
    # #v_IS is also 0 
    @model.Constraint(model.s_NS, model.s_Q, model.s_T)
    def Cinvzero_pretreatment(m,n,q,t):
        for a in m.s_Ain[n]:
            if a[0].endswith('TW'):
                return m.v_C[n,q,t] == 0
            else:
                return pyo.Constraint.Skip
    assert_units_consistent(model.Cinvzero_pretreatment)
    
    @model.Constraint(model.s_NS, model.s_Q, model.s_T)
    def Srelation(m,n,q,t):
        return scale_fac*sum(m.v_S[a, q, t] for a in m.s_Aout[n]) == scale_fac*m.v_C[n,q,t]*sum(m.v_F[a, t] for a in m.s_Aout[n])
    assert_units_consistent(model.Srelation)
    
    # Disposal
    @model.Constraint(model.s_ND, model.s_T, doc = "Disposal flow restriction")
    def Dflow(m,n,t):
        return sum(m.v_F[a,t] for a in m.s_Ain[n]) <= pyunits.convert(m.p_Cap[n], pyunits.liter/pyunits.s)
    assert_units_consistent(model.Dflow)
    
    # Freshwater 
    @model.Constraint(model.s_NW, model.s_T, doc = "Freshwater flow restriction")
    def Wflow(m,n,t):
        return sum(m.v_F[a,t] for a in m.s_Aout[n]) <= pyunits.convert(m.p_Cap[n], pyunits.liter/pyunits.s)
    assert_units_consistent(model.Wflow)
    
    @model.Constraint(model.s_NW, model.s_Q, model.s_T, doc = "Freshwater solids flow")
    def WSflow(m, n, q, t): 
        return sum(m.v_S[a, q, t] for a in m.s_Aout[n]) == 0
    assert_units_consistent(model.WSflow)
    
    # Treatment Inlet Node
    @model.Constraint(model.s_NTIN, model.s_T, doc = "Treatment inlet flow restriction")
    def TINflowmax(m,n,t):
        return sum(m.v_F[a,t] for a in m.s_Ain[n]) <= 20*pyunits.liter/pyunits.s
    assert_units_consistent(model.TINflowmax)
    
    # Treatment inlet node minimum flow restriction
    @model.Constraint(model.s_NTIN, model.s_T, doc = "Treatment inlet flow minimum")
    def TINflowmin(m,n,t):
        return sum(m.v_F[a,t] for a in m.s_Ain[n]) >= 3*pyunits.liter/pyunits.s
    assert_units_consistent(model.TINflowmin)
    
    # Relate concentration to solids flow into the treatment unit
    @model.Constraint(model.s_NTIN, model.s_Q, model.s_T, doc = "Linking solids to treatment inlet to inlet concentration")
    def TINconc_rel(m, n, q, t):
        return scale_fac*sum(model.v_F[a, t] for a in model.s_Ain[n])*m.v_C[n, q, t] == scale_fac*sum(model.v_S[a, q, t] for a in m.s_Ain[n])
    assert_units_consistent(model.TINconc_rel)
    
    #Make a total concentration variable for the treatment unit since it can only handle one component 
    @model.Constraint(model.s_NTIN, model.s_T, doc = "Linking total conc to component conc")
    def TINtotalconc(m, n, t):
        return scale_fac*sum(model.v_F[a, t] for a in model.s_Ain[n])*m.v_Ctreatin[n, t] == scale_fac*sum(sum(model.v_S[a, q, t] for a in m.s_Ain[n]) for q in m.s_Q)
    assert_units_consistent(model.TINtotalconc)
    
    #Bound on minimum total concentration
    @model.Constraint(model.s_NTIN, model.s_T, doc = "Minimum bound on inlet concentration of the treatment unit")
    def CINmin(m, n, t):
        return model.v_Ctreatin[n, t] >= 70*pyunits.g/pyunits.liter
    assert_units_consistent(model.CINmin)
    
    # Treatment Treated Water Node
    @model.Constraint(model.s_NTTW, model.s_T, doc = "Flow of Treated Water Stream")
    def NTTWflow(m,n,t):
        assert len(m.s_Aout[n]) == 1
        n1 = m.p_treatedIN[n]
        return m.v_F[m.s_Aout[n][0],t] == m.v_alphaW[n1, t]*(m.v_F[m.s_Ain[n1][0],t])
    assert_units_consistent(model.NTTWflow)
    
    model.NTTWconc = pyo.ConstraintList(doc = "Concentration of components being treated at Treated Water Node")
    for n in model.s_NTTW:
        assert len(model.s_Aout[n]) == 1
        n1 = model.p_treatedIN[n]
        for t in model.s_T:
            for q in model.s_Qalpha[n1]:
                model.NTTWconc.add(model.v_C[n,q,t] == 0)
    assert_units_consistent(model.NTTWconc)
    #We need a constraint saying salt out ot TTW is 0. Link conc with salt here too.
    
    # Treatment Concentrated Water Node
    @model.Constraint(model.s_NTCW, model.s_T, doc = "Flow of Concentrated Water Stream")
    def NTCWflow(m,n,t):
        assert len(m.s_Aout[n]) == 1
        n1 = m.p_treatedIN[n]
        return m.v_F[m.s_Aout[n][0],t] == (1-m.v_alphaW[n1, t])*(m.v_F[m.s_Ain[n1][0],t])
    assert_units_consistent(model.NTCWflow)

    model.NTCWconc = pyo.ConstraintList(doc = "Concentration of components being treated at Concentrated Water Node")
    for n in model.s_NTCW:
        assert len(model.s_Aout[n]) == 1
        n1 = model.p_treatedIN[n]
        for t in model.s_T:
            for q in model.s_Qalpha[n1]:
                model.NTCWconc.add(scale_fac*model.v_C[n,q,t]*model.v_F[model.s_Aout[n][0], t] == scale_fac*model.v_C[n1,q,t]*model.v_F[model.s_Ain[n1][0], t])
    assert_units_consistent(model.NTCWconc)
    
    #Relating outlet solids flow to concentration at the node
    @model.Constraint(model.s_NTCW|model.s_NTTW, model.s_Q, model.s_T, doc = "Outlet solids flow")
    def solids_flow(m, n, q, t):
        return scale_fac*sum(m.v_S[a, q, t] for a in m.s_Aout[n]) == scale_fac*m.v_C[n,q,t]*sum(m.v_F[a,t] for a in m.s_Aout[n])
    assert_units_consistent(model.solids_flow)
    
    # Treatment: Cases where treatment is acting as a splitter
    model.NTSrcconc = pyo.ConstraintList(doc = "Concentration of components not being treated")
    for n in model.s_NTTW | model.s_NTCW:
        assert len(model.s_Aout[n]) == 1
        n1  = model.p_treatedIN[n]
        for t in model.s_T:
            for q in (model.s_Q - model.s_Qalpha[n1]):
                model.NTSrcconc.add(model.v_C[n,q,t] == model.v_C[n1,q,t])
    assert_units_consistent(model.NTSrcconc)
            
    # Restricting the minimum concentration at the concentrated water node
    @model.Constraint(model.s_NTCW, model.s_Q, model.s_T, doc = "Minimum concentration out of concentrated water node")
    def minconccon(m,n,q,t):
        return m.v_C[n,q,t] >= m.p_Cmin[n,q]
    assert_units_consistent(model.minconccon)
    
    # Bounds on the salt in the lines
    model.saltbound = pyo.ConstraintList(doc ="Bounding salt flow in inlet of nodes")
    for a in model.s_A:
        for q in model.s_Q:
            for t in model.s_T:
                model.saltbound.add(model.v_S[a, q, t] <= model.p_Cmax*model.v_F[a, t])
    assert_units_consistent(model.saltbound)
    
    return model

def mee_svr_model(I=2, flag = 1):
    #Create model
    m = pyo.ConcreteModel()
    
    "Parameter definitions"
    #Number of evaporator effects
    N_evap = I
    
    #Number of compression stages
    N_compr = 1
    
    #Number of components
    N_components = 1
    
    if flag == 1:
        #Flowrate of feed stream
        m.flow_feed = pyo.Param(initialize = 10,
                              units = pyo.units.kg/pyo.units.s,
                              mutable = True)
        
        #Salt in the feed stream
        m.salt_feed = pyo.Param(initialize =  100,
                              units = pyo.units.g/pyo.units.kg,
                              mutable = True)
    else:
        #Flowrate of feed stream
        m.flow_feed = pyo.Var(initialize = 20,
                              units = pyo.units.kg/pyo.units.s,
                              bounds = (0,100))
        
        #Salt in the feed stream
        m.salt_feed = pyo.Var(initialize = 80,
                              units = pyo.units.g/pyo.units.kg,
                              bounds = (70,210))
        
    
    #Temperature of the feed stream
    m.feed_temperature = pyo.Param(initialize = 25,
                                   units = pyo.units.C)
    
    #Salt in brine concentrate
    m.salt_outlet_spec = pyo.Param(initialize = 250,
                                    units = pyo.units.g/pyo.units.kg, mutable = True)
    #Water recovery
    m.water_recovery_fraction = pyo.Var(initialize = 0.3,
                                        bounds = (0.1,0.9))
    
    #Temperature constraint parameters
    m.DT_min = pyo.Param(initialize = 2, units = pyo.units.C)
    m.DT_min_1 = pyo.Param(initialize = 2, units = pyo.units.C)
    m.DT_min_2 = pyo.Param(initialize = 2, units = pyo.units.C)
    m.DT_min_stage = pyo.Param(initialize = 2, units = pyo.units.C)
   
    #Antoine coefficients
    m.a = pyo.Param(initialize = 12.98437)
    m.b = pyo.Param(initialize = -2001.77468)
    m.c = pyo.Param(initialize = 139.61335)
    
    #Specific heat capacity of vapor
    m.cp_vapor = pyo.Param(initialize = 1.873)
    
    #Minimum pressure delta between evaporators
    m.DP_min = pyo.Param(initialize = 0.1) #==============>NEW
    
    #Overall heat transfer coefficient (Known parameter)
    m.overall_heat_transfer_coef =  pyo.Param(initialize = 100)
    
    #Heat capacity ratio
    m.gamma = pyo.Param(initialize = 1.33)
    
    #Maximum compression ratio
    m.CR_max = pyo.Param(initialize = 4)
    
    #Efficiency of compressors (isentropic efficiency)
    m.eta = pyo.Param(initialize = 0.75, mutable = True)
    
    #Costing parameters
    m.r = pyo.Param(initialize = 0.1)
    m.years = pyo.Param(initialize = 10)
    m.cepci_2022 = pyo.Param(initialize = 813)
    m.cepci_2003 = pyo.Param(initialize = 402)
    
    "Set definitions"
    #Set of Evaporator effects
    m.i = pyo.Set(initialize = range(N_evap))
    i_first = m.i.first()
    i_last = m.i.last()
    
    #Set of all evaporators except the first 
    m.i_except_1 = pyo.Set(initialize = range(1,N_evap))
    
    #Set of compression stages
    m.j = pyo.Set(initialize = range(N_compr))
    j_first = m.j.first()
    j_last = m.j.last()
    
    m.j_except_1 = pyo.Set(initialize = range(1,N_compr))
    
    #Set of components in water
    m.components = pyo.Set(initialize = range(N_components))
    
    "Variable definitions"
    #====================================================================
                            #All flow variables
    #Flow of brine evaporator
    m.flow_brine = pyo.Var(m.i, 
                           bounds = (1e-20, 100), 
                           initialize = [pyo.value(m.flow_feed)*pyo.value(1-m.water_recovery_fraction)]*I,
                           units = pyo.units.kg/pyo.units.s)
    #[pyo.value(m.flow_feed)*pyo.value(1-m.water_recovery_fraction)]*I
    #Flow of vapor evaporator
    m.flow_vapor_evaporator = pyo.Var(m.i, 
                           bounds = (1e-20, 100), 
                           initialize = [pyo.value(m.flow_feed - m.flow_brine[i_last])]*I,
                           units = pyo.units.kg/pyo.units.s)
    #[pyo.value(m.flow_feed - m.flow_brine[i_last])]*I
    #Flow of super heated vapor
    m.flow_super_heated_vapor = pyo.Var(bounds = (1e-20, 100),
                                        initialize = pyo.value(m.flow_vapor_evaporator[i_last]),
                                        units = pyo.units.kg/pyo.units.s)
    
    #Flow of treated water
    #Just wanted to make it easy to find the treated water variable in the network
    #Writing an equality constraint between flow_vapor from last evaporator and 
    #treated water variable
    m.flow_treated_water = pyo.Var( bounds = (1e-20, 100), 
                            initialize = pyo.value(m.flow_vapor_evaporator[i_last]),
                            units = pyo.units.kg/pyo.units.s)
    
    m.flow_brine_last = pyo.Var(bounds = (1e-20, 100),
                                initialize = 1,
                                units = pyo.units.kg/pyo.units.s)
    #======================================================================
                           #All concentration variables
    #Flow of salt/TDS (For now we are considering only one component in water)
    m.salt = pyo.Var(m.i, 
                     bounds = (1e-20, 300), 
                     initialize = [100]*I,
                     units = pyo.units.g/pyo.units.kg)
    
    #Salt mass fraction (XS in paper)
    m.salt_mass_frac= pyo.Var(m.i,
                              bounds = (1e-20,1),
                              initialize = 0.5)
    #[pyo.value(m.salt[i])/1000 for i in range(0,I)]
    #Salt mass fraction in the feed 
    #Question
    #maybe shouldn't be a parameter. There must be some relationship between 
    #salt_mass_frac_feed amd salt_feed
    m.salt_mass_frac_feed= pyo.Var(bounds= (1e-10,1),
                                   initialize = pyo.value(m.salt_feed)/1000)
    #pyo.value(m.salt_feed)/1000
    #======================================================================
                          #All pressure variables
    #Vapor pressure in evaporator effects
    m.evaporator_vapor_pressure = pyo.Var(m.i,
                                          bounds = (1, 200),
                                          initialize = [70]*I,
                                          units = pyo.units.kPa)
    m.super_heated_vapor_pressure = pyo.Var(m.j,
                                            bounds = (1,200),
                                            initialize = 60.540,
                                            units = pyo.units.kPa)
    #====================Not sure if we need this (Can this be replaced by evaporator vapor pressure?)
    #This is used to say that the vapor of the i-1 th evaporator should be hotter than the vapor of the ith 
    #evaporator for it to be able to evaporate the feed in the ith evaporator
    
    m.saturated_vapor_pressure= pyo.Var(m.i,
                                        bounds = (1, 300),
                                        initialize = 101.00,
                                        units = pyo.units.kPa) 
        
    #=======================================================================
                          #All temperature variables
    #Actal temperature of feed entering the evaporator after preheating
    m.evaporator_feed_temperature = pyo.Var(bounds = (1, 200),
                                            initialize = 25,
                                            units = pyo.units.C)
    
    m.evaporator_ideal_temperature = pyo.Var(m.i,
                                             bounds = (1, 200),
                                             initialize = 25,
                                             units = pyo.units.C)
    m.evaporator_brine_temperature = pyo.Var(m.i,
                                             bounds = (1, 200),
                                             initialize = 35,
                                             units = pyo.units.C)
    m.evaporator_condensate_temperature = pyo.Var(m.i,
                                                  bounds = (1, 200),
                                                  initialize  = 30,
                                                  units = pyo.units.C)
    m.super_heated_vapor_temperature = pyo.Var(m.j,
                                               bounds = (1,300),
                                               initialize = 50,
                                               units = pyo.units.C)
    
    #This is linked with saturated vapor pressure. Not sure if we need this either
    m.evaporator_saturated_vapor_temperature = pyo.Var(m.i,
                                                       bounds = (1, 200),
                                                       initialize = 30,
                                                       units = pyo.units.C)
    
    m.fresh_water_temperature = pyo.Var(bounds = (1, 200),
                                        initialize = 25,
                                        units = pyo.units.C)
    
    m.LMTD = pyo.Var(m.i,
                     initialize = 1,
                     units = pyo.units.C)
    
    m.theta_1 = pyo.Var(m.i,
                        bounds = (1, 200),
                        initialize = 1,
                        units = pyo.units.C)
    
    m.theta_2 = pyo.Var(m.i,
                        bounds = (1, 200),
                        initialize = 1,
                        units = pyo.units.C)
    
    m.preheater_LMTD = pyo.Var(initialize = 1,
                               units = pyo.units.C)
    
    m.preheater_theta_1 = pyo.Var(bounds = (1, 200),
                                  initialize = 1,
                                  units = pyo.units.C)
    
    m.preheater_theta_2 = pyo.Var(bounds = (1, 200),
                                  initialize = 1,
                                  units = pyo.units.C)
    
    m.isentropic_temperature = pyo.Var(m.j,
                                       bounds = (1e-20, 200),
                                       initialize = 25,
                                       units = pyo.units.C)
    m.mixer_temperature = pyo.Var(bounds = (1e-20, 200),
                                  initialize = 25,
                                  units = pyo.units.C)
    
    #=======================================================================
                        #All enthalpy variables
    #Specific Enthalpies of brine and vapor in the evaporator
    m.evaporator_brine_enthalpy = pyo.Var(m.i,
                                          domain = pyo.Reals,
                                          initialize = 400)
    
    m.evaporator_vapor_enthalpy = pyo.Var(m.i,
                                          domain = pyo.Reals,
                                          initialize = 400)
    #Defined only for first evaporator
    m.evaporator_condensate_vapor_enthalpy = pyo.Var(domain = pyo.Reals,
                                                     initialize = 100)
    
    m.evaporator_condensate_enthalpy = pyo.Var(domain = pyo.Reals,
                                               initialize = 100)
    
    m.super_heated_vapor_enthalpy = pyo.Var(m.j,
                                      domain = pyo.Reals,
                                      initialize = 100)
    #Enthalpy of the feed stream
    m.enthalpy_feed = pyo.Var(domain = pyo.Reals,
                              initialize = 200)
    
    #=======================================================================
                            #All area variables
    m.evaporator_total_area = pyo.Var(bounds = (14, 1500),
                                               initialize = 60,
                                               units = pyo.units.m**2)
    m.each_evaporator_area = pyo.Var(m.i,
                                     bounds = (14,372),
                                     initialize = 30,
                                     units = pyo.units.m**2)
    
    #=======================================================================
                            #Other Variables
    #Evaporator heat flow - Q
    m.evaporator_heat_flow = pyo.Var(m.i, 
                                     bounds = (1e-20, None),
                                     initialize = 1)
    #Boiling point elevation
    m.bpe = pyo.Var(m.i,
                    bounds = (1e-20, 200),
                    initialize = [15]*I)
    
    #Latent heat
    m.latent_heat = pyo.Var(m.i_except_1,
                            bounds = (1e-20, None),
                            initialize = 1)
    
    #Heat transfer coefficient
    m.heat_transfer_coef = pyo.Var(m.i,
                    bounds = (1e-20, None),
                    initialize = 1)
    
    #Specific heat capacity of condensate
    m.cp_condensate = pyo.Var(bounds = (1e-20, None),
                              initialize = 1)
    
    #Specific heat capacity of feed 
    m.cp_feed =  pyo.Var(bounds = (1e-20, None),
                         initialize = 1)
    
    #Preheater area
    m.preheater_area = pyo.Var(bounds = (1, 372),
                               initialize = 20)
    
    #preheater heat transfer coefficient
    m.preheater_heat_transfer_coef = pyo.Var(bounds = (1e-20, None),
                                             initialize = 0.1)
    
    #Each compressor work
    m.compressor_work = pyo.Var(m.j,
                                bounds = (1e-20, 30000),
                                initialize = 1)
    
    m.total_compressor_work = pyo.Var(bounds = (1e-20, 30000),
                                      initialize = 1)
    #=======================================================================
                        #Variables for costing
    m.compressor_capacity = pyo.Var(m.j,
                                    bounds = (1e-20,30000),
                                    initialize = 150)
    
    m.annual_fac = pyo.Var(bounds = (1e-20, None),
                           initialize = 1)
    
    m.cepci_ratio = pyo.Var(bounds = (1e-20, None),
                           initialize = 1)
    
    m.compressor_capex = pyo.Var(bounds = (1e-20, None),
                                 initialize = 1)
    
    m.evaporator_capex = pyo.Var(bounds = (1e-20, None),
                                 initialize = 1)
    
    m.preheater_capex = pyo.Var(bounds = (1e-20, None),
                                initialize = 1)
    
    m.CAPEX = pyo.Var(bounds = (1e-20, None),
                      initialize=1)
    
    m.OPEX = pyo.Var(bounds = (1e-20, None),
                      initialize=1)
     
    #=======================================================================
    "Model Constraints"
    
    #=======================================================================
                            #Evaporator Constraints
    #Flow balance across the evaporator (Equation 1, 3)
    def _evaporator_flow_balance(m, i):
        if i != i_last:
            return m.flow_brine[i+1] - m.flow_brine[i] - m.flow_vapor_evaporator[i] == 0
        else: 
            return m.flow_feed - m.flow_brine[i] -m.flow_vapor_evaporator[i] == 0 
    
    m.evaporator_flow_balance = pyo.Constraint(m.i, rule = _evaporator_flow_balance)
    
    #Linking treated water variable with vapor flow from last evaporator
    def _treated_water_linking_cons(m):
        return m.flow_treated_water == sum(m.flow_vapor_evaporator[i] for i in m.i)
    m.treated_water_linking_cons = pyo.Constraint(rule = _treated_water_linking_cons)
    
    #Linking brine water from the last evaporator to the flow_brine_last var
    def _brine_water_linking_cons(m):
        return m.flow_brine_last == m.flow_brine[i_last]
    m.brine_water_linking_cons = pyo.Constraint(rule = _brine_water_linking_cons)
    
    #Solid balance in the evaporator (Equation 3, 4)
    def _evaporator_salt_balance(m, i):
        if i != i_last:
            return m.flow_brine[i+1]*m.salt[i+1] - m.flow_brine[i]*m.salt[i] == 0
        else:
            return m.flow_feed*m.salt_feed - m.flow_brine[i]*m.salt[i] == 0
    
    m.evaporator_salt_balance = pyo.Constraint(m.i, rule = _evaporator_salt_balance)
    
    #Estimating the ideal temperature in an evaporator (Equation 5)
    def _ideal_evaporator_temperature(m,i):
        return pyo.log(m.evaporator_vapor_pressure[i]) ==\
            m.a + m.b/(m.evaporator_ideal_temperature[i] +m.c) 
    
    m.evaporator_ideal_temp_con = pyo.Constraint(m.i, rule = _ideal_evaporator_temperature)
    
    #Boiling point elevation (Equation 6)
    def _bpe_con(m, i):
        return m.bpe[i] == 0.1581 \
            + 2.769*m.salt_mass_frac[i]\
            - 0.002676*m.evaporator_ideal_temperature[i]\
            + 41.78*m.salt_mass_frac[i]**2 \
            + 0.134*m.salt_mass_frac[i]*m.evaporator_ideal_temperature[i]
    
    m.bpe_con = pyo.Constraint(m.i, rule = _bpe_con)
    
    #Relating mass fraction of salt to brine salinity (Equation 7)
    def _match_mass_frac_to_salinity(m, i):
        return m.salt_mass_frac[i] - 0.001*m.salt[i] == 0
    m.match_mass_frac_to_salinity = pyo.Constraint(m.i, rule = _match_mass_frac_to_salinity)
    
    #Relating salt mass frac of feed to salinity of feed
    def _match_mass_frac_to_salinity_feed(m):
        return m.salt_mass_frac_feed == 0.001*m.salt_feed
    m.match_mass_frac_to_salinity_feed = pyo.Constraint(rule = _match_mass_frac_to_salinity_feed)
    
    #Relate brine temperature to ideal temperature and bpe (Equation 8)
    def _brine_temp_con(m, i):
        return m.evaporator_brine_temperature[i] - m.evaporator_ideal_temperature[i]\
            - m.bpe[i] == 0
    
    m.brine_temp_con = pyo.Constraint(m.i, rule = _brine_temp_con)
    
    #Energy balance in the evaporator (Equations 9, 10)
    def _evaporator_energy_balance(m, i):
        if i != i_last:
            return 1/1000*(m.evaporator_heat_flow[i]\
                + m.flow_brine[i+1]*m.evaporator_brine_enthalpy[i+1]\
                - m.flow_brine[i]*m.evaporator_brine_enthalpy[i]\
                - m.flow_vapor_evaporator[i]*m.evaporator_vapor_enthalpy[i]) == 0
        else:
            return 1/1000*(m.evaporator_heat_flow[i]\
                + m.flow_feed*m.enthalpy_feed\
                - m.flow_brine[i]*m.evaporator_brine_enthalpy[i]\
                - m.flow_vapor_evaporator[i]*m.evaporator_vapor_enthalpy[i]) == 0
        
    m.evaporator_energy_balance = pyo.Constraint(m.i, rule = _evaporator_energy_balance)
    
    #Estimating the enthalpies (Equations 11, 12, 13)
    def _enthalpy_vapor_estimate(m, i):
        return 1/1000*(m.evaporator_vapor_enthalpy[i]) == 1/1000*(-13470 + 1.84*m.evaporator_brine_temperature[i])
    m.enthalpy_vapor_estimate = pyo.Constraint(m.i, rule = _enthalpy_vapor_estimate)
    
    def _enthalpy_brine_estimate(m, i):
        return 1/1000*(m.evaporator_brine_enthalpy[i]) == 1/1000*(-15940 + 8787*m.salt_mass_frac[i] + 3.557*m.evaporator_brine_temperature[i])
    m.enthalpy_brine_estimate = pyo.Constraint(m.i, rule = _enthalpy_brine_estimate)
    
    def _enthalpy_feed_estimate(m):
        return 1/1000*m.enthalpy_feed == 1/1000*(-15940 + 8787*m.salt_mass_frac_feed + 3.557*m.evaporator_feed_temperature)
    m._enthalpy_feed_estimate = pyo.Constraint(rule = _enthalpy_feed_estimate)
    
    #==============================================================================
    #Estimating the condensate temperature in the first evaporator from the outlet compressor pressure
    #Not sure if this constraint is needed. Should be calculated by energy balance
    def _condensate_temperature_estimate(m, i):
        if i == i_first:
            return pyo.log(m.super_heated_vapor_pressure[j_last])\
                - m.a - m.b/(m.evaporator_condensate_temperature[i] +m.c) == 0
        else:
            return pyo.log(m.evaporator_vapor_pressure[i - 1])\
                - m.a - m.b/(m.evaporator_condensate_temperature[i] +m.c) == 0
            
    m.condensate_temperature_estimate = pyo.Constraint(m.i, rule = _condensate_temperature_estimate)
    
    #Estimating enthalpy of the condensate vapor and condensate 
    def _enthalpy_condensate_vapor_estimate(m):
        return 1/1000*(m.evaporator_condensate_vapor_enthalpy) == 1/1000*(-13470 + 1.84*m.evaporator_condensate_temperature[i_first])
    m.enthalpy_condensate_vapor_estimate = pyo.Constraint(rule = _enthalpy_condensate_vapor_estimate)
    
    def _enthalpy_condensate_estimate(m):
        return 1/1000*(m.evaporator_condensate_enthalpy) == 1/1000*(-15940 + 3.557*m.evaporator_condensate_temperature[i_first])
    m.enthalpy_condensate_estimate = pyo.Constraint(rule = _enthalpy_condensate_estimate)
    
    #Flow balance for super heated vapor (Equation 15)
    m.flow_balance_super_heated_vapor = pyo.Constraint(expr = m.flow_super_heated_vapor
                                                       == m.flow_vapor_evaporator[i_last])
    
    
    #Heat requirements in evaporators
    def _evaporator_heat_balance(m,i):
        if i == i_first:
            #Heat requirements in the 1st evaporator (Equation 14)
            return 1/1000*(m.evaporator_heat_flow[i_first]) ==\
                1/1000*(m.flow_super_heated_vapor*m.cp_vapor*\
                    (m.super_heated_vapor_temperature[j_last] - m.evaporator_condensate_temperature[i_first])+\
                        m.flow_super_heated_vapor*\
                    (m.evaporator_condensate_vapor_enthalpy- m.evaporator_condensate_enthalpy))
    
        else:
            #Heat requirements in other evaporators (Equation 16)
            return 1/1000*(m.evaporator_heat_flow[i]) == 1/1000*(m.flow_vapor_evaporator[i-1]*m.latent_heat[i])
    m.evaporator_heat_balance = pyo.Constraint(m.i, rule = _evaporator_heat_balance)
    
    #Calculating latent heat from temperature (Equation 17)
    def _latent_heat_estimation(m, i):
        if i == i_first:
            return pyo.Constraint.Skip
        return m.latent_heat[i] == 2502.5 - 2.3648*m.evaporator_saturated_vapor_temperature[i] + \
            1.840*(m.evaporator_saturated_vapor_temperature[i-1] - m.evaporator_saturated_vapor_temperature[i])
    m.latent_heat_estimation = pyo.Constraint(m.i, rule = _latent_heat_estimation)
    
    #Calculating saturated vapor temperature from saturated vapor pressure
    def _saturated_vapor_temp_estimate(m, i):
        if I == 1:
            return pyo.Constraint.Skip
        return pyo.log(m.saturated_vapor_pressure[i]) \
            - m.a - m.b/(m.evaporator_saturated_vapor_temperature[i] +m.c) == 0
    m.saturated_vapor_temp_estimate = pyo.Constraint(m.i, rule = _saturated_vapor_temp_estimate)
    
    #Pressure and temperature feasiility constraints (Equation 18)
    def _pressure_gradient_constraint(m, i):
        if i != i_last:
            return m.evaporator_vapor_pressure[i] >= m.evaporator_vapor_pressure[i+1] + m.DP_min
        else:
            return pyo.Constraint.Skip
    m.pressure_gradient_constraint = pyo.Constraint(m.i, rule = _pressure_gradient_constraint)
    
    #Relating evaporator vapor pressure with saturated vapor pressure (Equation 19)
    def _relating_pressure_in_evaporator(m, i):
        if i != i_first:
            return m.saturated_vapor_pressure[i] == m.evaporator_vapor_pressure[i-1]
        else:
            return m.saturated_vapor_pressure[i] == m.super_heated_vapor_pressure[j_last]
    m.relating_pressure_in_evaporator = pyo.Constraint(m.i, rule = _relating_pressure_in_evaporator)
    
    #Calculating the heat transfer coefficient (Equation 20)
    def _heat_transfer_coef_calculation(m,i):
        return m.heat_transfer_coef[i] == 0.001*(1939.4 + 1.40562*m.evaporator_brine_temperature[i]
                                                 - 0.00207525*m.evaporator_brine_temperature[i]**2
                                                 + 0.0023186*m.evaporator_brine_temperature[i]**3)
    m.heat_transfer_coef_calculation = pyo.Constraint(m.i, rule = _heat_transfer_coef_calculation)
    
    #Evaporator heat transfer area calculation (Equation 21)
    def _total_evaporator_heat_transfer_area(m):
        return m.evaporator_total_area == sum(m.each_evaporator_area[i] for i in m.i)
    m.total_evaporator_heat_transfer_area = pyo.Constraint(rule = _total_evaporator_heat_transfer_area)
    
    #Area of the first evaporator
    def _first_evaporator_area_calculation(m):
        return m.each_evaporator_area[i_first] == m.flow_super_heated_vapor*m.cp_vapor*\
            (m.super_heated_vapor_temperature[j_last] - m.evaporator_condensate_temperature[i_first])\
            /(m.overall_heat_transfer_coef*m.LMTD[i_first])\
            + m.flow_super_heated_vapor*(m.evaporator_condensate_vapor_enthalpy - m.evaporator_condensate_enthalpy)\
            /(m.heat_transfer_coef[i_first]*(m.evaporator_condensate_temperature[i_first] - m.evaporator_brine_temperature[i_first])+1e-6)
            
    m.first_evaporator_area_calculation = pyo.Constraint(rule = _first_evaporator_area_calculation)
    
    def _evaporator_total_area_from_heat_calculation(m):
        if I == 1:
            return pyo.Constraint.Skip
        else:
            return m.evaporator_total_area == sum((m.evaporator_heat_flow[i]/(m.heat_transfer_coef[i]*m.LMTD[i])) for i in m.i)
    m.evaporator_total_area_from_heat_calculation = pyo.Constraint(rule = _evaporator_total_area_from_heat_calculation)
    
    #Chen approximation for LMTD (Equation 22-26)
    #For the other evaporators should theta 1 use T_vapor_evaporator/saturated temperatures
    #Paper uses saturated temperatures
    def _theta_1_calculation(m, i):
        if i==i_first:
            return m.theta_1[i] == m.super_heated_vapor_temperature[j_last] - m.evaporator_brine_temperature[i]
        else:
            return m.theta_1[i] == m.evaporator_saturated_vapor_temperature[i] - m.evaporator_brine_temperature[i]
    m.theta_1_calculation = pyo.Constraint(m.i, rule = _theta_1_calculation)
    
    def _theta_2_calculation(m, i):
        if I == 1:
            return m.theta_2[i] == m.evaporator_condensate_temperature[i] -  m.evaporator_feed_temperature     
        else:
            if i == i_first:
                return m.theta_2[i] == m.evaporator_condensate_temperature[i] - m.evaporator_brine_temperature[i+1]
            elif i == i_last:
                return m.theta_2[i] == m.evaporator_condensate_temperature[i] - m.evaporator_feed_temperature
            else:
                return m.theta_2[i] == m.evaporator_condensate_temperature[i] -  m.evaporator_brine_temperature[i+1]     
    m.theta_2_calculation = pyo.Constraint(m.i, rule = _theta_2_calculation)
            
    def _LMTD_calculation(m, i):
        alpha = m.theta_1[i]/m.theta_2[i]
        eps = 1e-10
        return m.LMTD[i] == m.theta_2[i]*((alpha -1)**2 + eps)**0.5/(pyo.log(alpha)**2 + eps)**0.5
        #return m.LMTD[i] == (0.5*m.theta_1[i]*m.theta_2[i]*(m.theta_1[i]+m.theta_2[i]))**(1/3)
    m.LMTD_calculation = pyo.Constraint(m.i, rule = _LMTD_calculation)
    
    #Restrictions on area for uniform distribution (Equation 27, 28)
    def _area_restriction_con_1(m, i):
        if I == 1:
            return pyo.Constraint.Skip
        else:
            return m.each_evaporator_area[i] <= 3*m.each_evaporator_area[i-1]
    m.area_restricion_con_1 = pyo.Constraint(m.i_except_1, rule = _area_restriction_con_1)
    
    def _area_restriction_con_2(m, i):
        if I == 1:
            return pyo.Constraint.Skip
        else:
            return m.each_evaporator_area[i] >= 1*m.each_evaporator_area[i-1]
    m.area_restricion_con_2 = pyo.Constraint(m.i_except_1, rule = _area_restriction_con_2)
    
    #Temperature constraints to avoid temperature crossovers in evaporator effects (Equation 29-36)
    def _temp_con_1(m):
        return m.super_heated_vapor_temperature[j_last] >= m.evaporator_condensate_temperature[i_first] + m.DT_min_1
    m.temp_con_1 = pyo.Constraint(rule = _temp_con_1)
    
    def _temp_con_2(m, i):
        if I == 1:
            return pyo.Constraint.Skip
        else:
            return m.evaporator_brine_temperature[i-1] >= m.evaporator_condensate_temperature[i] + m.DT_min_1
    m.temp_con_2 = pyo.Constraint(m.i_except_1, rule = _temp_con_2)
    
    def _temp_con_3(m, i):
        if i == i_last:
            return pyo.Constraint.Skip
        else:        
            return m.evaporator_brine_temperature[i] >= m.evaporator_brine_temperature[i+1] + m.DT_min_stage
    m.temp_con_3 = pyo.Constraint(m.i, rule = _temp_con_3)
    
    def _temp_con_4(m):
        return m.evaporator_brine_temperature[i_last] >= m.evaporator_feed_temperature + m.DT_min_2
    m.temp_con_4 = pyo.Constraint(rule = _temp_con_4)
    
    def _temp_con_5(m, i):
        if i == i_last:
            return pyo.Constraint.Skip
        else:
            return m.evaporator_condensate_temperature[i] >= m.evaporator_brine_temperature[i+1] + m.DT_min
    m.temp_con_5 = pyo.Constraint(m.i, rule = _temp_con_5)
    
    def _temp_con_6(m):
        return m.evaporator_condensate_temperature[i_last] >= m.evaporator_feed_temperature + m.DT_min
    m.temp_con_6 = pyo.Constraint(rule = _temp_con_6)
    
    def _temp_con_7(m,i):
        return m.evaporator_condensate_temperature[i] >= m.evaporator_brine_temperature[i] + m.DT_min
    m.temp_con_7 = pyo.Constraint(m.i, rule = _temp_con_7)
    
    def _temp_con_8(m,i):
        if I == 1:
            return pyo.Constraint.Skip
        else:
            return m.evaporator_saturated_vapor_temperature[i] >= m.evaporator_brine_temperature[i] + m.DT_min
    m.temp_con_8 = pyo.Constraint(m.i, rule = _temp_con_8)
    
    #=====================================================================================
                                        #Preheater constraints
    #Write mixer constraints to determine the temperature of the inlet stream to the preheater
    #Energy balance in the mixer (Used a similar form of Equation 50)
    #Check if this is the correct way to determine the mixer temperature
    #This is just taking the weighted average to calculate the temperature of the mixed stream
    def _mixer_energy_balance(m):
        return m.mixer_temperature == sum(m.flow_vapor_evaporator[i]*
                                         m.evaporator_condensate_temperature[i] for i in m.i)/sum(m.flow_vapor_evaporator[i] for i in m.i)
        # return m.flow_vapor_evaporator[i_first]*(m.mixer_temperature - m.evaporator_condensate_temperature[i_first]) ==\
        #     m.flow_vapor_evaporator[i_last]*(m.evaporator_condensate_temperature[i_last] - m.mixer_temperature)
    m.mixer_energy_balance = pyo.Constraint(rule = _mixer_energy_balance)
    
    #=================================================Checked until here
    #Why do we have evaporator ideal temperature in the energy balance for flash
    #Energy balance in the preheater (Equation 43)
    def _energy_balance_pre_heater(m):
        return sum (m.flow_vapor_evaporator[i] for i in m.i)*m.cp_condensate*(m.mixer_temperature- m.fresh_water_temperature)\
            == m.flow_feed*m.cp_feed*(m.evaporator_feed_temperature - m.feed_temperature)
    m.energy_balance_pre_heater = pyo.Constraint(rule = _energy_balance_pre_heater)
    
    #Cp_feed calculation (Equation 44)
    def _cp_feed_estimate(m):
        return m.cp_feed == 0.001*(4206.8-6.6197*m.salt_feed + 1.2288e-2*m.salt_feed**2 +\
                                    (-1.1262 + 5.418e-2*m.salt_feed)*m.feed_temperature)
    m.cp_feed_estimate = pyo.Constraint(rule = _cp_feed_estimate)
    
    #Cp_condensate calculation (Equation 45)
    #Why are we using evaporator's ideal temperature to calculate this quantity
    #I believe we should use the temperature of the stream that is going into the 
    #preheater. But I am not sure if we can use this cp correlation then
    
    def _cp_condensate_estimate(m):
        return m.cp_condensate == 0.001*(4206.8 - 1.1262*m.mixer_temperature)
    m.cp_condensate_estimate = pyo.Constraint(rule = _cp_condensate_estimate)
    
    #Preheater Heat tranfer coefficient calculation
    def _preheater_heat_transfer_coef_con(m):
        return m.preheater_heat_transfer_coef == 0.001*(1939.4 + 1.40562*m.mixer_temperature
                                                 - 0.00207525*m.mixer_temperature**2
                                                 + 0.0023186*m.mixer_temperature**3)
    m.preheater_heat_transfer_coef_con = pyo.Constraint(rule = _preheater_heat_transfer_coef_con)
    
    #Heat transfer area of feed preheater(Equation 46)
    def _preheater_area_calculation(m):
        return m.preheater_area == sum(m.flow_vapor_evaporator[i] for i in m.i)*m.cp_condensate*\
            (m.mixer_temperature - m.fresh_water_temperature)/(m.preheater_heat_transfer_coef*m.preheater_LMTD + 1e-20)
    m.preheater_area_calculation = pyo.Constraint(rule = _preheater_area_calculation)
    
    #Preheater LMTD calculation (Equation 47)
    def _preheater_LMTD_calculation(m):
        alpha = m.preheater_theta_1/m.preheater_theta_2
        eps = 1e-10
        return m.preheater_LMTD== m.preheater_theta_2*((alpha -1)**2 + eps)**0.5/(pyo.log(alpha)**2 + eps)**0.5
        #return m.preheater_LMTD == (0.5*m.preheater_theta_1*m.preheater_theta_2*(m.preheater_theta_1+m.preheater_theta_2))**(1/3)
    m.preheater_LMTD_calculation = pyo.Constraint(rule = _preheater_LMTD_calculation)
    
    def _preheater_theta_1_calculation(m):
        return m.preheater_theta_1 == m.mixer_temperature-m.evaporator_feed_temperature
    m.preheater_theta_1_calculation = pyo.Constraint(rule = _preheater_theta_1_calculation)
    
    def _preheater_theta_2_calculation(m):
        return m.preheater_theta_2 == m.fresh_water_temperature - m.feed_temperature
    m.preheater_theta_2_calculation = pyo.Constraint(rule = _preheater_theta_2_calculation)
    
    #=====================================================================================
                                    #Single stage Compressor
    #Isentropic temperature constraints (Equation 51, 52)
    def _isentropic_temp_calculation(m, j):
        if j == j_first:
            return m.isentropic_temperature[j] == (m.evaporator_brine_temperature[i_last] + 273.15)*\
                (m.super_heated_vapor_pressure[j]/m.evaporator_vapor_pressure[i_last])**((m.gamma -1)/m.gamma) - 273.15
    m.isentropic_temp_calculation = pyo.Constraint(m.j, rule = _isentropic_temp_calculation)
    
    #Maximum possible compression (Equation 53)
    def _maximum_compression_calculation(m, j):
        if j == j_first:
            return m.super_heated_vapor_pressure[j] <= m.CR_max*m.evaporator_vapor_pressure[i_last]
    m.maximum_compression_calculation = pyo.Constraint(m.j, rule = _maximum_compression_calculation)
    
    #Temperature of superheated vapor
    def _temperature_super_heated_vapor_calculation(m, j):
        if j == j_first: 
            return m.super_heated_vapor_temperature[j] == m.evaporator_brine_temperature[i_last] +\
                1/m.eta*(m.isentropic_temperature[j] - m.evaporator_brine_temperature[i_last])
    m.temperature_super_heated_vapor_calculation = pyo.Constraint(m.j, rule = _temperature_super_heated_vapor_calculation)
    
    #(Equation 56, 57)
    def _compressor_pressure_con(m, j):
        if j ==j_first:
            return m.super_heated_vapor_pressure[j]  >= m.evaporator_vapor_pressure[i_last]
    m.compressor_pressure_con = pyo.Constraint(m.j, rule = _compressor_pressure_con)
    
    #Compressor work calculation (Equation 58)
    def _compressor_work_calculation(m, j):
        return 1/100*(m.compressor_work[j]) == 1/100*(m.flow_super_heated_vapor*(m.super_heated_vapor_enthalpy[j] - m.evaporator_vapor_enthalpy[i_last]))
    m.compressor_work_calculation = pyo.Constraint(m.j, rule = _compressor_work_calculation)
    
    def _super_heated_vapor_enthalpy_calculation(m, j):
        return 1/1000*m.super_heated_vapor_enthalpy[j]== (-13470 + 1.84*m.super_heated_vapor_temperature[j])/1000
    m.super_heated_vapor_enthalpy_calculation = pyo.Constraint(m.j, rule = _super_heated_vapor_enthalpy_calculation)
    
    def _total_compressor_work(m):
        return m.total_compressor_work == sum(m.compressor_work[j] for j in m.j)
    m._total_compressor_work = pyo.Constraint(rule = _total_compressor_work)
    
    # m.delta_pos = pyo.Var(domain = pyo.NonNegativeReals)
    # m.delta_neg = pyo.Var(domain = pyo.NonPositiveReals)
    #m.delta.fix(0)
    #Salt outlet condition
    def _salt_outlet_con(m):
        return m.salt[i_first] == m.salt_outlet_spec #+ m.delta_pos + m.delta_neg
    m.salt_outlet_con = pyo.Constraint(rule = _salt_outlet_con)
    
    #Fresh water constraint
    def _water_recovery_con(m):
        return sum(m.flow_vapor_evaporator[i] for i in range(0,I)) == m.water_recovery_fraction*m.flow_feed
    m.water_recovery_con = pyo.Constraint(rule = _water_recovery_con)
    
    #Costing===========================================================================
    #Costing is taken from the couper book same as used in Onishi's paper

    #Converting compressor_capacity from kW to HP 
    def _comp_capacity_con(m,j):
        return m.compressor_capacity[j] == m.compressor_work[j]*1.34
    m.comp_capacity_con = pyo.Constraint(m.j, rule = _comp_capacity_con)
    
    def _annual_fac_con(m):
        return m.annual_fac == m.r*(m.r + 1)**m.years/((1+m.r)**m.years - 1)
    m.annual_fac_con = pyo.Constraint(rule = _annual_fac_con)
    
    def _cepci_ratio_con(m):
        return m.cepci_ratio == m.cepci_2022/m.cepci_2003
    m.cepci_ratio_con = pyo.Constraint(rule = _cepci_ratio_con)
    
    #All costs are in kUS$
    def _compressor_capex_con(m):
        return m.compressor_capex == m.cepci_ratio*(sum(7.9*m.compressor_capacity[j]**0.62 for j in range(N_compr)))
    m.compressor_capex_con = pyo.Constraint(rule = _compressor_capex_con)
    
    def _evaporator_capex_con(m):
        return m.evaporator_capex == m.cepci_ratio*1.218*sum(pyo.exp(3.2362 - 0.0126*pyo.log(m.each_evaporator_area[i]*10.764) + 0.0244*(pyo.log(m.each_evaporator_area[i]*10.764))**2) for i in range(N_evap))
    m.evaporator_capex_con = pyo.Constraint(rule = _evaporator_capex_con)
    
    def _preheater_capex_con(m):
        a = m.preheater_area*10.764
        fd = pyo.exp(-0.9186 + 0.0830*pyo.log(a)) #utube heat exchanger
        fp = 1 #pressure < 4 bar
        fm = 1 #carbon steel
        Cb = pyo.exp(8.821-0.306*pyo.log(a) + 0.0681*(pyo.log(a))**2)
        return m.preheater_capex ==  m.cepci_ratio*1.218/1000*fd*fm*fp*Cb
    m.preheater_capex_con = pyo.Constraint(rule = _preheater_capex_con)
    
    m.capital_cost = pyo.Constraint(expr = m.CAPEX == m.annual_fac*
                                    (m.compressor_capex + m.evaporator_capex + m.preheater_capex))
    
    m.treatment_reward = pyo.Var(initialize = 1)
    m.treatment_reward_con = pyo.Constraint(expr = m.treatment_reward == 15.216310215*sum(m.flow_vapor_evaporator[i] for i in m.i))
    m.operation_cost = pyo.Constraint(expr = m.OPEX == 1.4*m.total_compressor_work)
    m.obj = pyo.Objective(expr = (m.CAPEX*7/365 + m.OPEX*7/365 
                                  - m.treatment_reward))
  
    return m     

def integrated_model_build(solve_network_problem = False, I = 1):
    #Network Model
    file_dir = os.path.dirname(__file__)
    fname = os.path.join(file_dir, "data.pickle")
    with open(fname, 'rb') as f:
        data = pickle.load(f)
    m_network = build_qcp(data)
    
    m_network.TINflowmin['R02_IN', :].deactivate()
    m_network.TINflowmax['R02_IN', :].deactivate()
    m_network.CINmin['R02_IN', :].deactivate()
    
    m_network.v_alphaW['R02_IN', :].fix(0.4)
    
    if solve_network_problem:
        solver = pyo.SolverFactory('gams')
        #ipopt.options['max_iter'] = 10000
        res = solver.solve(m_network, solver = 'conopt', tee= True)
       
        try:
            pyo.assert_optimal_termination(res)
        except:
            print("Constraint lower bound violated")
            for c in m_network.component_data_objects(ctype = pyo.Constraint):
                if c.lb is not None and pyo.value(c.body) <= c.lb - 1e-3:
                    print(c.name, pyo.value(c.body), c.lb)
            print("Constraint upper bound violated")
            for c in m_network.component_data_objects(ctype = pyo.Constraint):
                if c.ub is not None and pyo.value(c.body)>= c.ub + 1e-3:
                    print(c.name, pyo.value(c.body), c.ub)
            print("Variable lower bound violated")
            for v in m_network.component_data_objects(ctype = pyo.Var):
                if v.lb is not None and pyo.value(v)<= v.lb - 1e-3:
                    print(v.name, pyo.value(v), v.lb)
            print("Constraint upper bound violated")   
            import pdb;pdb.set_trace()
            
    for v in m_network.component_data_objects(ctype = pyo.Var):
        if v.ub is not None and pyo.value(v) >= v.ub + 1e-3:
            print(v.name, pyo.value(v), v.ub)
    #import pdb;pdb.set_trace()
    #Treatment Model
    N_evap = I
    m_treatment = mee_svr_model(I = N_evap, flag = 0)
    
    #==============================================================================
               #Getting the treatment models for each time periods
    treatment_models = {}
    for i in m_network.s_T:
        treatment_models[i] = mee_svr_model(I = N_evap, flag = 0)
    
    #Making a concrete model for the integrated model
    m = pyo.ConcreteModel()
   
    
    m.m_network = m_network
    m.m_treatment = pyo.Reference(treatment_models)
    #==============================================================================
                    #Add global variables for capacities
    m.CAPEX = pyo.Var(initialize = 1000, domain = pyo.PositiveReals)
    #Global variable for evaporator area
    m.global_evaporator_capex = pyo.Var(initialize = 1000, domain = pyo.PositiveReals)
    
    #Global preheater area
    m.global_preheater_capex = pyo.Var(initialize = 100, domain = pyo.PositiveReals)
    
    #Global compressor capacity
    m.global_compressor_capex = pyo.Var(initialize = 100, domain = pyo.PositiveReals)
    
    #===============================================================================
        #Processing for differentiating pretreatment units from desalintion units
    
    #R02 is a pretreatment unit. So we fix the recovery fraction to 40%
    #Note the cost for pretreatment is not considered. 
    #Need to think more about this. This is leading to sometimes infeasible solutions
    m.m_network.v_alphaW['R02_IN', :].fix(0.4)
    
    #Deactivate minimum flow and concentration requirements for the pretreatment unit
    m.m_network.TINflowmin['R02_IN', :].deactivate()
    m.m_network.TINflowmax['R02_IN', :].deactivate()
    m.m_network.CINmin['R02_IN', :].deactivate()
    
    #==============================================================================
                    #Linking network model with treatment units
    #Link feed flow to treatment in each time period
    def _feed_flow_link(m, t):
        return m.m_treatment[t].flow_feed == m.m_network.v_F['N03', 'R01_IN', 'Piped', t]
    m.feed_flow_link = pyo.Constraint(m.m_network.s_T, rule = _feed_flow_link)
    
    #Link feed composition to treatment feed composition
    def _feed_conc_link(m, t):
        return m.m_treatment[t].salt_feed == m.m_network.v_C['R01_IN','TDS', t]
    m.feed_conc_link = pyo.Constraint(m.m_network.s_T, rule = _feed_conc_link)
    
    #Link recovery fractions
    def _recovery_fractions_link(m, t):
        return m.m_network.v_alphaW['R01_IN', t] == m.m_treatment[t].water_recovery_fraction
    m.recovery_fractions_link = pyo.Constraint(m.m_network.s_T, rule = _recovery_fractions_link)
    
    #Inequality constraints linking global vars to local treatment vars
    def _global_evap_area_link(m, t):
        return m.global_evaporator_capex >= m.m_treatment[t].evaporator_capex
    m.global_evap_area_link = pyo.Constraint(m.m_network.s_T, rule = _global_evap_area_link)
    
    def _global_ph_area_link(m, t):
        return m.global_preheater_capex >= m.m_treatment[t].preheater_capex
    m.global_preheater_area_link = pyo.Constraint(m.m_network.s_T, rule = _global_ph_area_link)
    
    def _global_comp_capacity_link(m, t):
        return m.global_compressor_capex >= m.m_treatment[t].compressor_capex
    m.global_compressor_capacity_link = pyo.Constraint(m.m_network.s_T, rule = _global_comp_capacity_link)
    
    #Calculating CAPEX based on the global variables
    
    m.capex_con = pyo.Constraint(expr = m.CAPEX == m.m_treatment['T01'].annual_fac*(m.global_evaporator_capex + m.global_compressor_capex + m.global_preheater_capex))
    #Objective for the integrated problem
    #The periods are defined over a year and the newtork costs and treatment opex and capex are annualized
    #Therefore it is necessary to scale them by the number of periods so that we know 
    #that the final costs are still $/year 
    #(This ofcourse is given the flow parameters otherwise it should be reported in terms of some $/flow unit)
    
    #Deactivate objectives
    m.m_network.obj.deactivate()
    [m.m_treatment[t].obj.deactivate() for t in m.m_network.s_T]
    m.obj =  pyo.Objective(expr = 1e-4*(m.m_network.obj
                           + sum(m.m_treatment[t].OPEX*m.m_network.p_dt for t in m.m_network.s_T)/365
                           + m.CAPEX/365*m.m_network.p_dt*len(m.m_network.s_T)))
    return m

if __name__ =="__main__":
    m = integrated_model_build()
    
    ipopt = pyo.SolverFactory('ipopt')
    ipopt.options["max_iter"] = 10000
    ipopt.solve(m, tee = True)