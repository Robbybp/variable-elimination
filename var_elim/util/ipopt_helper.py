import pyomo.environ as pyo
import os

def ipopt_solver_timing_records(m, tee= True):
    with open("ipopt.opt", "w") as f:
        f.write("print_timing_statistics yes")
    ipopt = pyo.SolverFactory('ipopt')
    ipopt.options["output_file"] = "ipopt.out"
    res = ipopt.solve(m, tee = tee)
    pyo.assert_optimal_termination(res)
    with open("ipopt.out", "r") as f: 
        lines = f.readlines()
        
    for line in lines:
        if line.startswith("Total number of variables"):
            num_vars = line.split(':')[-1]
            num_vars.strip()
        elif line.startswith("Total number of equality constraints"):
            num_cons = line.split(":")[-1]
            num_cons.strip()
        elif line.startswith("Total number of inequality constraints"):
            num_ineq =  line.split(":")[-1]
            num_ineq.strip()
        elif line.startswith("Total CPU secs in IPOPT (w/o function evaluations)"):
            ipopt_lin_solve = line.split("=")[-1]
            ipopt_lin_solve.strip()
        elif line.startswith("Total CPU secs in NLP function evaluations"):
            nlp_eval = line.split("=")[-1]
            nlp_eval.strip()
        elif line.startswith("Number of nonzeros in equality constraint Jacobian"):
            nlp_nze = line.split(":")[-1]
            nlp_nze.strip()
        elif line.startswith("Number of nonzeros in inequality constraint Jacobian"):
            nlp_nzie = line.split(":")[-1]
            nlp_nzie.strip()
        elif line.startswith("Number of nonzeros in Lagrangian Hessian"):
            nlp_nzlh = line.split(":")[-1]
            nlp_nzlh.strip()
        elif line.startswith("Number of Iterations"):
            nlp_iter = line.split(":")[-1]
            nlp_iter.strip()
   
    return int(num_vars), int(num_cons), int(num_ineq), float(ipopt_lin_solve), float(nlp_eval), int(nlp_nze), int(nlp_nzie), int(nlp_nzlh), int(nlp_iter) 