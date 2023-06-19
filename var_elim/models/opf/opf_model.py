#  ___________________________________________________________________________
#
#  Variable Elimination: Research code for variable elimination in NLPs
#
#  Copyright (c) 2023. Triad National Security, LLC. All rights reserved.
#
#  This program was produced under U.S. Government contract 89233218CNA000001
#  for Los Alamos National Laboratory (LANL), which is operated by Triad
#  National Security, LLC for the U.S. Department of Energy/National Nuclear
#  Security Administration. All rights in the program are reserved by Triad
#  National Security, LLC, and the U.S. Department of Energy/National Nuclear
#  Security Administration. The Government is granted for itself and others
#  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
#  in this material to reproduce, prepare derivative works, distribute copies
#  to the public, perform publicly and display publicly, and to permit others
#  to do so.
#
#  This software is distributed under the 3-clause BSD license.
#  ___________________________________________________________________________

import os
from egret.parsers.matpower_parser import create_ModelData
from egret.models.acopf import create_psv_acopf_model


def get_model_from_file(fname):
    model_data = create_ModelData(fname)
    model, _ = create_psv_acopf_model(model_data)
    return model


def make_model(fname=None):
    if fname is None:
        dirname = os.path.dirname(__file__)
        # By default, build model corresponding to 5k-bus OPF model
        # from grid optimization competition
        fname = os.path.join(dirname, "pglib_opf_case4917_goc.m")
    return get_model_from_file(fname)


if __name__ == "__main__":
    import pyomo.environ as pyo
    m = make_model()
    solver = pyo.SolverFactory("ipopt")
    solver.solve(m, tee=True)
