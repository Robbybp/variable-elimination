#  ___________________________________________________________________________
#
#  variable elimination: research code for variable elimination in nlps
#
#  copyright (c) 2023. triad national security, llc. all rights reserved.
#
#  this program was produced under u.s. government contract 89233218cna000001
#  for los alamos national laboratory (lanl), which is operated by triad
#  national security, llc for the u.s. department of energy/national nuclear
#  security administration. all rights in the program are reserved by triad
#  national security, llc, and the u.s. department of energy/national nuclear
#  security administration. the government is granted for itself and others
#  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
#  in this material to reproduce, prepare derivative works, distribute copies
#  to the public, perform publicly and display publicly, and to permit others
#  to do so.
#
#  this software is distributed under the 3-clause bsd license.
#  ___________________________________________________________________________

import argparse
from collections import namedtuple
import os
import matplotlib
import matplotlib.pyplot as plt
from pyomo.common.collections import ComponentSet
from pyomo.util.subsystems import TemporarySubsystemManager
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix
import var_elim.scripts.config as config
from var_elim.models.testproblems import DistillationTestProblem
from var_elim.dae_utils import DAEInterface


TIME_SETS_LOOKUP = {
    "distill": ("t",),
    "pipeline": ("fs.time",),
    "mb-steady": ("fs.moving_bed.solid_phase.length_domain", "fs.moving_bed.gas_phase.length_domain"),
}


DOF_LOOKUP = {
    "distill": ["u1"],
}


def get_dae_components(m, time_sets, omit=None):
    set_elements = list(time_sets[0])
    # Sets have to have the same elements 
    assert all(list(s) == set_elements for s in time_sets)
    if omit is not None:
        set_elements = [i for i in set_elements if i not in omit]
    interfaces = [DAEInterface(m, s) for s in time_sets]

    discretization_cons = []
    differential_vars = []
    algebraic_cons = []
    algebraic_vars = []
    differential_cons = []
    derivative_vars = []

    for interface in interfaces:
        for i, t in enumerate(set_elements):
            variables, constraints = interface.get_subsystem_at_time(t)
            algvars, algcons = interface.get_algebraic_subsystem_at_time(t)
            derivvars, diffcons = interface.get_differential_subsystem_at_time(t)
            seen = ComponentSet(algvars + algcons + derivvars + diffcons)
            diffvars = [v for v in variables if v not in seen]
            disccons = [c for c in constraints if c not in seen]

            # These should be lists-of-lists
            discretization_cons.append(disccons)
            differential_vars.append(diffvars)
            algebraic_cons.append(algcons)
            algebraic_vars.append(algvars)
            differential_cons.append(diffcons)
            derivative_vars.append(derivvars)

    # I want to return an "indexed list of subsystems"
    DaeTuple = namedtuple("DaeTuple", ["discretization", "algebraic", "differential"])
    subsystems = [
        DaeTuple(
            (discretization_cons[i], differential_vars[i]),
            (algebraic_cons[i], algebraic_vars[i]),
            (differential_cons[i], derivative_vars[i]),
        )
        for i, t in enumerate(set_elements)
    ]
    return subsystems


def _get_rectangle_around_coords(ij1, ij2, linewidth=2, linestyle="-"):
    i1, j1 = ij1
    i2, j2 = ij2
    buffer = 0.5
    ll_corner = (min(i1, i2) - buffer, min(j1, j2) - buffer)
    width = abs(i1 - i2) + 2 * buffer
    height = abs(j1 - j2) + 2 * buffer
    rect = matplotlib.patches.Rectangle(
        ll_corner,
        width,
        height,
        clip_on=False,
        fill=False,
        edgecolor="orange",
        linewidth=linewidth,
        linestyle=linestyle,
    )
    return rect


def plot_dae_order(m, time_sets, dof):
    omit = (time_sets[0].first(), time_sets[0].last())
    with TemporarySubsystemManager(to_fix=dof):
        subsystems = get_dae_components(m, time_sets, omit=omit)

    corder = []
    vorder = []
    alg_coords = []
    for sub in subsystems:
        corder.extend(sub.discretization[0])
        vorder.extend(sub.discretization[1])

        alg_start = len(corder), len(vorder)
        corder.extend(sub.algebraic[0])
        vorder.extend(sub.algebraic[1])
        # I think I expect these coords to be inclusive
        alg_stop = len(corder)-1, len(vorder)-1
        alg_coords.append((alg_start, alg_stop))

        corder.extend(sub.differential[0])
        vorder.extend(sub.differential[1])
    vorder.extend(dof)
    imat = get_structural_incidence_matrix(vorder, corder)

    rectangles = [
        _get_rectangle_around_coords(start, stop)
        for start, stop in alg_coords
    ]

    fig, ax = plt.subplots()
    ax.set_ylabel("Constraints")
    ax.set_xlabel("Variables")
    ax.xaxis.set_label_position("top")
    ax.spy(
        imat,
        markersize=1.5,
    )
    ax.tick_params(length=0)
    for rect in rectangles:
        ax.add_patch(rect)
    return fig, ax


def main(args):
    constructor_lookup = dict(config.CONSTRUCTOR_LOOKUP)
    # Construct a small instance
    constructor_lookup["distill"] = DistillationTestProblem(nfe=5, horizon=5).create_instance
    m = constructor_lookup[args.model]()
    assert args.model in TIME_SETS_LOOKUP
    assert args.model in DOF_LOOKUP, "This script only supports model=distill for now"
    time_set = [m.find_component(s) for s in TIME_SETS_LOOKUP[args.model]]
    dof = [m.find_component(v) for v in DOF_LOOKUP[args.model]]
    dof = [vdata for v in dof for vdata in v.values()]
    fig, ax = plot_dae_order(m, time_set, dof)

    fname = f"{args.model}-daeorder.pdf"
    fpath = os.path.join(args.image_dir, fname)
    if not args.no_save:
        print(f"Saving fig to {fpath}")
        fig.savefig(fpath, transparent=not args.opaque)
    else:
        print(f"no-save set. Would have saved to {fpath}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    argparser = config.get_plot_argparser()
    args = argparser.parse_args()
    plt.rcParams["font.size"] = 14
    plt.rcParams["font.family"] = "sans-serif"
    main(args)
