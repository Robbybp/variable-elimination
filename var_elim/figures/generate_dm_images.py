import pyomo.environ as pyo
from pyomo.common.collections import ComponentMap
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import (
    get_incidence_graph,
    get_structural_incidence_matrix,
)

import networkx.drawing.nx_pylab as nxpl
import networkx.drawing.layout as nx_layout
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import scipy.sparse as sps

from var_elim.figures.dm_model import make_model


plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 18


def project_onto(coo, rows, cols):
    row_set = set(rows)
    col_set = set(cols)
    new_row = []
    new_col = []
    new_data = []
    for r, c, e in zip(coo.row, coo.col, coo.data):
        if r in row_set and c in col_set:
            new_row.append(r)
            new_col.append(c)
            new_data.append(e)
    return sps.coo_matrix((new_data, (new_row, new_col)), shape=coo.shape)


def get_rectangle_around_coords(ij1, ij2):
    i1, j1 = ij1
    i2, j2 = ij2
    buffer = 0.5
    ll_corner = (min(i1, i2)-buffer, min(j1, j2)-buffer)
    width = abs(i1 - i2) + 2*buffer
    height = abs(j1 - j2) + 2*buffer
    rect = Rectangle(
        ll_corner,
        width,
        height,
        clip_on=False,
        fill=False,
        edgecolor="orange",
        linewidth=2,
    )
    return rect


def generate_dm_images(show=True, save=False, transparent=True):
    m = make_model()
    igraph = IncidenceGraphInterface(m)
    var_dmp, con_dmp = igraph.dulmage_mendelsohn()

    variables = (
        var_dmp.unmatched
        + var_dmp.underconstrained
        + var_dmp.square
        + var_dmp.overconstrained
    )
    constraints = (
        con_dmp.underconstrained
        + con_dmp.square
        + con_dmp.overconstrained
        + con_dmp.unmatched
    )

    var_idx_map = ComponentMap((var, i) for i, var in enumerate(variables))
    con_idx_map = ComponentMap((con, i) for i, con in enumerate(constraints))

    matrix = get_structural_incidence_matrix(variables, constraints)

    fig, (ax1, ax2) = plt.subplots(1, 2)

    text = r"""
    \[\begin{array}{c} \\ \mathcal{B}_1 \\ \mathcal{B}_2 \\ \mathcal{B}_3 \\ \end{array}\hspace{-0.2cm} \begin{array}{c} \begin{array}{ccc} \mathcal{A}_1 & \mathcal{A}_2 & \mathcal{A}_3\end{array} \\ \left[\begin{array}{ccc} S_1 & * & * \\ & S_2 & * \\ & & S_3 \\ \end{array}\right] \\ \end{array} \begin{array}{c} \\ \\ = \\ \\ \end{array} \]
    """
    #text = r"""\[\begin{array}{c} a \\ c \\ d \\ \end{array}\]"""
    ax1.set_axis_off()
    ax1.text(0.5, 0.20, text)

    ax2.tick_params(length=0)
    markersize = 6
    ax2.spy(matrix, markersize=markersize)
    tick_locs = [0, 5, 10, 15]
    tick_vals = [r"$\mathrm{0}$", r"$\mathrm{5}$", r"$\mathrm{10}$", r"$\mathrm{15}$"]
    ax2.xaxis.set_ticks(tick_locs, tick_vals)
    ax2.yaxis.set_ticks(tick_locs, tick_vals)

    w, h = fig.get_size_inches()
    fig.set_size_inches(w, h*0.6)
    fig.tight_layout()

    uc_start = (0, 0)
    uc_stop = (
        len(var_dmp.unmatched + var_dmp.underconstrained)-1,
        len(con_dmp.underconstrained)-1,
    )
    rect = get_rectangle_around_coords(uc_start, uc_stop)
    ax2.add_patch(rect)

    x, y = uc_stop
    sq_start = (x+1, y+1)
    sq_stop = (x + len(var_dmp.square), y + len(con_dmp.square))
    rect = get_rectangle_around_coords(sq_start, sq_stop)
    ax2.add_patch(rect)

    x, y = sq_stop
    oc_start = (x+1, y+1)
    oc_stop = (
        x + len(var_dmp.overconstrained),
        y + len(con_dmp.overconstrained) + len(con_dmp.unmatched),
    )
    rect = get_rectangle_around_coords(oc_start, oc_stop)
    ax2.add_patch(rect)

    #subsystems = []
    #colors = []

    #underconstrained_vars = var_dmp.unmatched + var_dmp.underconstrained
    #underconstrained_cons = con_dmp.underconstrained
    #subsystems.append((underconstrained_cons, underconstrained_vars))
    #colors.append("orange")

    #square_vars = var_dmp.square
    #square_cons = con_dmp.square
    #subsystems.append((square_cons, square_vars))
    #colors.append("green")

    #overconstrained_vars = var_dmp.overconstrained
    #overconstrained_cons = con_dmp.overconstrained + con_dmp.unmatched
    #subsystems.append((overconstrained_cons, overconstrained_vars))
    #colors.append("red")

    #for i, (constraints, variables) in enumerate(subsystems):
    #    color = colors[i]
    #    rows = [con_idx_map[con] for con in constraints]
    #    cols = [var_idx_map[var] for var in variables]
    #    proj_matrix = project_onto(matrix, rows, cols)
    #    plt.spy(proj_matrix, color=color, markersize=markersize)

    if save:
        fig.savefig("bkg_dm_order.pdf", transparent=transparent)
    if show:
        plt.show()


if __name__ == "__main__":
    generate_dm_images(save=True, show=False)
