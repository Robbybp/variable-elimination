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


DIRNAME = os.path.dirname(__file__)


datafile_dict = {
    "ieee118": os.path.join(DIRNAME, "pglib_opf_case118_ieee.m"),
    "goc4917": os.path.join(DIRNAME, "pglib_opf_case4917_goc.m"),
    "goc10000": os.path.join(DIRNAME, "pglib_opf_case10000_goc.m"),
}


def get_model_from_file(fname):
    model_data = create_ModelData(fname)
    model, _ = create_psv_acopf_model(model_data)
    return model


def make_model(instance_id="goc4917", pf_bounds=True):
    fname = datafile_dict[instance_id]
    model = get_model_from_file(fname)

    # Turn off unnecessary bounds on power flow variables
    if not pf_bounds:
        model.pf.setlb(None)
        model.pf.setub(None)
        model.pt.setlb(None)
        model.pt.setub(None)
        model.qf.setlb(None)
        model.qf.setub(None)
        model.qt.setlb(None)
        model.qt.setub(None)
    return model


if __name__ == "__main__":
    import pyomo.environ as pyo
    m = make_model()
    solver = pyo.SolverFactory("ipopt")
    solver.solve(m, tee=True)
