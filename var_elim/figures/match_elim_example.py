import pyomo.environ as pyo
from pyomo.common.collections import ComponentMap
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix

import scipy.sparse as sps
import matplotlib.pyplot as plt

from var_elim.models.distillation.distill import create_instance as create_distill
from var_elim.heuristics.matching import (
    generate_elimination_via_matching,
    define_elimination_order,
)
from var_elim.algorithms.replace import eliminate_variables 


def solve_original(tee=True):
    m = create_distill()
    var_elim, con_elim = generate_elimination_via_matching(m)
    var_elim, con_elim = define_elimination_order(var_elim, con_elim)
    eliminate_variables(m, var_elim, con_elim)
    solver = pyo.SolverFactory("ipopt")
    solver.solve(m, tee=tee)
    return m


def solve_reduced(tee=True):
    m = create_distill()
    solver = pyo.SolverFactory("ipopt")
    solver.solve(m, tee=tee)
    return m


def print_statistics(m):
    igraph = IncidenceGraphInterface(m)
    linear_igraph = IncidenceGraphInterface(
        m, include_inequality=False, linear_only=True
    )

    nvar = len(igraph.variables)
    ncon = len(igraph.constraints)

    nedge = igraph.n_edges
    nedge_linear = linear_igraph.n_edges

    matching = linear_igraph.maximum_matching()
    n_match = len(matching)

    matched_vars = list(matching.values())
    matched_cons = list(matching.keys())
    vblocks, cblocks = igraph.block_triangularize(matched_vars, matched_cons)
    nblocks = len(vblocks)
    singleton_block_indices = [i for i in range(nblocks) if len(vblocks[i]) == 1]
    nsingleton = len(singleton_block_indices)

    var_elim, con_elim = generate_elimination_via_matching(m)
    nelim = len(var_elim)

    print(f"Problem: {nvar} Var., {ncon} Con.")
    print(f"Incidence: {nedge} edges ({nedge_linear} linear)")
    print(f"Size of maximum linear matching: {n_match}")
    print(f"N. blocks: {nblocks} ({nsingleton} singletons)")
    print(f"N. variables in non-singleton blocks: {n_match-nsingleton}")
    print(f"Size of elimination order: {nelim}")


def print_dm_statistics(var_dmp, con_dmp):
    oc_var = var_dmp.overconstrained
    oc_con = con_dmp.overconstrained + con_dmp.unmatched
    uc_var = var_dmp.unmatched + var_dmp.underconstrained
    uc_con = con_dmp.underconstrained
    wc_var = var_dmp.square
    wc_con = con_dmp.square

    n_oc_var = len(oc_var)
    n_oc_con = len(oc_con)
    n_uc_var = len(uc_var)
    n_uc_con = len(uc_con)
    n_wc_var = len(wc_var)
    n_wc_con = len(wc_con)

    print(f"Under-constrained: {n_uc_con} Con., {n_uc_var} Var.")
    print(f"Over-constrained: {n_oc_con} Con., {n_oc_var} Var.")
    print(f"Well-constrained: {n_wc_con} Con., {n_wc_var} Var.")


def get_matched_entry_submatrix(variables, constraints, matching):
    var2idx = ComponentMap((var, i) for i, var in enumerate(variables))
    con2idx = ComponentMap((con, i) for i, con in enumerate(constraints))
    matched_coords = [(con2idx[con], var2idx[var]) for con, var in matching.items()]
    matched_coord_set = set(matched_coords)
    imat_full = get_structural_incidence_matrix(variables, constraints)
    sm_row = []
    sm_col = []
    sm_data = []
    for row, col, val in zip(imat_full.row, imat_full.col, imat_full.data):
        if (row, col) in matched_coord_set:
            sm_row.append(row)
            sm_col.append(col)
            sm_data.append(val)
    sm = sps.coo_matrix((sm_data, (sm_row, sm_col)), shape=imat_full.shape)
    return sm


def make_plots(m, show=False, save=False, file_prefix="imat", extension=".pdf"):
    igraph = IncidenceGraphInterface(m)
    linear_igraph = IncidenceGraphInterface(
        m, include_inequality=False, linear_only=True
    )
    #matching = linear_igraph.maximum_matching()
    #matched_vars = list(matching.values())
    #matched_cons = list(matching.keys())
    #vblocks, cblocks = igraph.block_triangularize(matched_vars, matched_cons)
    #var_elim, con_elim = generate_elimination_via_matching(m)

    lvar_dmp, lcon_dmp = linear_igraph.dulmage_mendelsohn()
    matched_vars = lvar_dmp.underconstrained + lvar_dmp.square + lvar_dmp.overconstrained
    matched_cons = lcon_dmp.underconstrained + lcon_dmp.square + lcon_dmp.overconstrained
    matching = ComponentMap(zip(matched_cons, matched_vars))
    #matching = linear_igraph.maximum_matching()
    print_dm_statistics(lvar_dmp, lcon_dmp)
    lvdmorder = sum(lvar_dmp, [])
    lcdmorder = sum(lcon_dmp, [])

    vblocks, cblocks = igraph.block_triangularize(matched_vars, matched_cons)
    lvbtorder = sum(vblocks, [])
    lcbtorder = sum(cblocks, [])

    var_elim, con_elim = generate_elimination_via_matching(m)
    var_elim, con_elim = define_elimination_order(var_elim, con_elim)

    imat_full = get_structural_incidence_matrix(igraph.variables, igraph.constraints)
    imat_linear = get_structural_incidence_matrix(lvdmorder, lcdmorder, linear_only=True)
    matching_sm = get_matched_entry_submatrix(lvdmorder, lcdmorder, matching)
    imat_bt = get_structural_incidence_matrix(lvbtorder, lcbtorder)
    imat_elim = get_structural_incidence_matrix(var_elim, con_elim)

    plt.rcParams["font.size"] = 18

    # TODO: Put in Dulmage-Mendelsohn order?
    fig, ax = plt.subplots()
    ax.spy(imat_full, markersize=1)
    ax.xaxis.set_tick_params(length=0)
    ax.yaxis.set_tick_params(length=0)
    if save:
        fig.savefig(file_prefix+"_full"+extension, transparent=True)

    fig, ax = plt.subplots()
    ax.spy(imat_linear, markersize=1)
    ax.spy(matching_sm, markersize=1, color="orange")
    ax.xaxis.set_tick_params(length=0)
    ax.yaxis.set_tick_params(length=0)
    if save:
        fig.savefig(file_prefix+"_linear"+extension, transparent=True)

    fig, ax = plt.subplots()
    ax.spy(imat_bt, markersize=1)
    ax.xaxis.set_tick_params(length=0)
    ax.yaxis.set_tick_params(length=0)
    if save:
        fig.savefig(file_prefix+"_bt"+extension, transparent=True)

    fig, ax = plt.subplots()
    ax.spy(imat_elim, markersize=1)
    ax.spy(sps.identity(len(var_elim)), markersize=1, color="orange")
    ax.xaxis.set_tick_params(length=0)
    ax.yaxis.set_tick_params(length=0)
    ax.set_xticks([0, 1000, 2000, 3000])
    ax.set_yticks([0, 1000, 2000, 3000])
    if save:
        fig.savefig(file_prefix+"_lt"+extension, transparent=True)

    if show:
        plt.show()


def main():
    #solve_original()
    #solve_reduced()
    m = create_distill()
    print_statistics(m)
    make_plots(m, show=False, save=True, file_prefix="distill", extension=".pdf")


if __name__ == "__main__":
    main()
