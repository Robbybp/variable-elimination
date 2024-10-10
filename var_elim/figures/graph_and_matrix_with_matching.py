import pyomo.environ as pyo
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import (
    get_structural_incidence_matrix,
    get_bipartite_incidence_graph,
)
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
import networkx as nx

from var_elim.figures.config import (
    get_var_node_color,
    get_con_node_color,
)


TRANSPARENT = True
FONTSIZE = 40
AXESOFF = True
LINEWIDTH = 4.0

VCOLOR = get_var_node_color()
CCOLOR = get_con_node_color()

plt.rcParams["font.size"] = FONTSIZE
plt.rcParams["text.usetex"] = True
#plt.rc("text.latex", preamble=r"\usepackage{amsmath}")


def main():
    m = pyo.ConcreteModel()
    m.x = pyo.Var()
    m.y = pyo.Var()
    m.z = pyo.Var()
    m.eq1 = pyo.Constraint(expr=m.x + m.y + m.z == 1)
    m.eq2 = pyo.Constraint(expr=m.x**2 - 2*m.y == 2)

    igraph = IncidenceGraphInterface(m)

    variables = [m.x, m.y, m.z]
    nvar = len(variables)
    var_labels = [r"$x$", r"$y$", r"$z$"]

    constraints = [m.eq1, m.eq2]
    ncon = len(constraints)
    con_labels = [r"$\mathrm{eq1}$", r"$\mathrm{eq2}$"]

    imat = get_structural_incidence_matrix(variables, constraints)

    graph = get_bipartite_incidence_graph(variables, constraints)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

    if AXESOFF:
        ax1.set_axis_off()
    text = r"""
    \[\begin{array}{cc} \mathrm{eq1:} & x + y + z = 1 \\ \mathrm{eq2:} & x^2 + 2y = 2 \\ \end{array}\rightarrow\]
    """
    ax1.set_xlim((0, 1))
    ax1.set_ylim((0, 1))
    text_xy = (-0.1, 0.18)
    ax1.text(*text_xy, text, fontsize=FONTSIZE*1.1)
    x, y = text_xy
    xrect_xy = (x + 0.42, y + 0.33)
    xrect_width = 0.16
    xrect_height = 0.13
    xrect = Rectangle(
        xrect_xy,
        xrect_width,
        xrect_height,
        fill=False,
        edgecolor="orange",
        linewidth=LINEWIDTH,
    )
    yrect_xy = (x + 0.95, y + 0.12)
    yrect_width = 0.23
    yrect_height = 0.17
    yrect = Rectangle(
        yrect_xy,
        yrect_width,
        yrect_height,
        fill=False,
        edgecolor="orange",
        linewidth=LINEWIDTH,
        clip_on=False,
    )
    ax1.add_patch(xrect)
    ax1.add_patch(yrect)

    text = r"""
    $\mathcal{M}=\{(\mathrm{eq1}, x), (\mathrm{eq2}, y)\}$
    """
    text_xy = (text_xy[0], text_xy[1] + 0.5)
    ax1.text(*text_xy, text, fontsize=FONTSIZE*0.9)


    ybuffer = 0.02
    var_rshift = 4.5
    ax2.set_xlim((-3, 2.5+var_rshift))
    ax2.set_ylim((-1-ybuffer, 1+ybuffer))
    top = list(range(ncon))
    layout = nx.bipartite_layout(graph, top)

    # Orders equations/variables bottom-up for some reason...
    # Flip layout around x axis.
    ymid = 0.0
    for node, (x, y) in list(layout.items()):
        if node < ncon:
            layout[node] = (x, ymid - y)
        else:
            layout[node] = (x+var_rshift, ymid - y)
    node_colors = [
        CCOLOR if node < ncon else VCOLOR for node in list(graph)
    ]
    nx.draw_networkx(
        graph,
        pos=layout,
        ax=ax2,
        node_size=600,
        width=2.0,
        with_labels=False,
        node_color=node_colors,
        style="--",
    )
    nx.draw_networkx_edges(
        graph,
        pos=layout,
        edgelist=[(0, 2), (1, 3)],
        width=5.0,
        ax=ax2,
        edge_color="orange",
    )

    label_layout = {}
    node_labels = {}
    var_offset = 2.5
    con_offset = 3
    for node, (x, y) in layout.items():
        if node >= ncon:
            label_layout[node] = (x + var_offset, y)
            node_labels[node] = var_labels[node - ncon]
        else:
            label_layout[node] = (x - con_offset, y)
            node_labels[node] = con_labels[node]
    nx.draw_networkx_labels(
        graph,
        pos=label_layout,
        labels=node_labels,
        ax=ax2,
        font_size=FONTSIZE,
        clip_on=False,
    )

    var_x, _ = layout[ncon]
    ellipse_w = 3
    ellipse_h = 2.0
    var_ellipse = Ellipse(
        (var_x, ymid),
        ellipse_w,
        ellipse_h,
        fill=False,
        edgecolor=VCOLOR,
        linewidth=LINEWIDTH,
    )
    con_x, _ = layout[0]
    con_ellipse = Ellipse(
        (con_x, ymid),
        ellipse_w,
        ellipse_h,
        fill=False,
        edgecolor=CCOLOR,
        linewidth=LINEWIDTH,
    )
    ax2.add_patch(var_ellipse)
    ax2.add_patch(con_ellipse)

    nodelabel_y = 1.2
    vlabel_offset = 2
    ax2.text(
        var_x-vlabel_offset,
        nodelabel_y,
        r"$\mathrm{Variables}$",
        fontsize=FONTSIZE*0.9,
    )
    clabel_offset = 5.5
    ax2.text(
        con_x-clabel_offset,
        nodelabel_y,
        r"$\mathrm{Constraints}$",
        fontsize=FONTSIZE*0.9,
    )
    ax2.text(10, -0.1, r"$\leftrightarrow$", fontsize=FONTSIZE*1.1)

    if AXESOFF:
        ax2.set_axis_off()

    # TODO: axis tick labels?
    ax3.tick_params(length=0)
    ax3.spy(imat, markersize=30)
    ax3.spy([[1, 0, 0],[0, 1, 0]], markersize=40, color="orange")
    xtick_pos = list(range(nvar))  # NOTE: This depends on implementation of spy
    ytick_pos = list(range(ncon))
    ax3.xaxis.set_ticks(xtick_pos, var_labels)
    ax3.yaxis.set_ticks(ytick_pos, con_labels)
    for axis in (ax3.xaxis, ax3.yaxis):
        for tick in axis.get_major_ticks():
            pad = tick.get_pad()
            tick.set_pad(4*pad)

    w, h = fig.get_size_inches()
    factor = 0.8
    fig.set_size_inches(3.3*w*factor, 1.5*h*factor)
    fig.tight_layout()
    fig.savefig("bkg_imat_matching.pdf", transparent=TRANSPARENT)


if __name__ == "__main__":
    main()
