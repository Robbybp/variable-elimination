#Create bipartite graph for matching heuristic demonstration
import pyomo.environ as pyo
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import (
    get_structural_incidence_matrix,
    get_bipartite_incidence_graph,
)
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import networkx as nx

from var_elim.figures.config import (
    get_var_node_color,
    get_con_node_color,
)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from var_elim.figures.schematics.matching_algorithm_for_figure_generation import generate_elimination_via_matching
from pyomo.common.collections import ComponentSet, ComponentMap
from var_elim.figures.match_elim_example import get_matched_entry_submatrix
import scipy.sparse as sps
TRANSPARENT = True
FONTSIZE = 40
AXESOFF = True
LINEWIDTH = 4.0

VCOLOR = get_var_node_color()
CCOLOR = get_con_node_color()

plt.rcParams["font.size"] = FONTSIZE
plt.rcParams["text.usetex"] = True
#plt.rc("text.latex", preamble=r"\usepackage{amsmath}")

def _make_model_for_nested_replacement():
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3, 4], initialize=1)
    m.y = pyo.Var([1, 2], initialize=1, bounds=(0, 10))

    m.eq1 = pyo.Constraint(expr=m.x[4] == m.x[2] + m.x[1] - m.y[1]**2)
    m.eq2 = pyo.Constraint(expr=m.x[2] == m.y[1]*m.y[2])
    m.eq3 = pyo.Constraint(expr=m.x[3] == m.x[1]*m.x[2] + m.x[4])
    m.eq4 = pyo.Constraint(expr=m.x[1] == m.x[2] - m.y[2])
    # Note that eq5 cannot be used for elimination.
    m.eq5 = pyo.Constraint(expr=m.x[1] + m.x[2] + m.x[3] * m.x[4] == 2)

    # Need this initialization for convergence of the reduced-space problem
    calculate_variable_from_constraint(m.x[2], m.eq2)
    calculate_variable_from_constraint(m.x[1], m.eq4)
    calculate_variable_from_constraint(m.x[4], m.eq1)
    calculate_variable_from_constraint(m.x[3], m.eq3)

    m.obj = pyo.Objective(expr=m.x[1]**2 + 2*m.x[2]**2 + 3*m.x[3]**2 + 4*m.x[4]**2)

    return m

def get_edges(var_list, con_list, con_var_list):
    edge_list = []
    node_id_map = ComponentMap((v, i) for i, v in enumerate(con_var_list))
   
    for v, c in zip(var_list, con_list):
        print(c.name, v.name)
        con_coord = node_id_map[c]
        var_coord = node_id_map[v]
        
        edge_list.append((con_coord, var_coord))
    return edge_list

def get_edges_subgraph(variables, constraints, subgraph, con_var_list):
    edge_list = []
    
    #These lists are for plotting the matrix and are here
    #just because we are computing the edge cars and cons here anyways 
    #It is probably better to
    var_on_edge = []
    con_on_edge = []
    node_id_map_original = ComponentMap((v, i) for i, v in enumerate(con_var_list))
    
    con_var_list_subgraph = constraints + variables
    node_id_map_new = ComponentMap((i, v) for i, v in enumerate(con_var_list_subgraph)) 
    
    for edge in subgraph.edges:
        con_id =  edge[0]
        var_id = edge[1]
        var_name = node_id_map_new[var_id]
        con_name = node_id_map_new[con_id]
        var_on_edge.append(var_name)
        con_on_edge.append(con_name)
        
        var_id_orig = node_id_map_original[var_name]
        con_id_orig = node_id_map_original[con_name]
        edge_list.append((con_id_orig, var_id_orig))
    
    return edge_list, var_on_edge, con_on_edge

def get_submatrix_incidence(variables, constraints, var_on_subgraph_edges, con_on_subgraph_edges):
    #Here variables and constraints are the original variables and constraints
    var2idx = ComponentMap((var, i) for i, var in enumerate(variables))
    con2idx = ComponentMap((con, i) for i, con in enumerate(constraints))
    matched_coords = [(con2idx[con], var2idx[var]) for con, var in zip(con_on_subgraph_edges, var_on_subgraph_edges)]
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
    

def main(highlight_matched_edges = False, induced_subgraph = False, matching_elim_plot = False, file_name = None, edge_colors = 'k'):
    m = _make_model_for_nested_replacement()
    
    igraph = IncidenceGraphInterface(m)

    variables = [m.x[1], m.x[2], m.x[3], m.x[4], m.y[1], m.y[2]]
    nvar = len(variables)
    var_labels = [r"$x_1$", r"$x_2$", r"$x_3$", r"$x_4$", r"$y_1$", r"$y_2$"]

    constraints = [m.eq1, m.eq2, m.eq3, m.eq4, m.eq5]
    ncon = len(constraints)
    con_labels = [r"$\mathrm{eq1}$", r"$\mathrm{eq2}$", r"$\mathrm{eq3}$", r"$\mathrm{eq4}$", r"$\mathrm{eq5}$"]

    imat = get_structural_incidence_matrix(variables, constraints)

    graph = get_bipartite_incidence_graph(variables, constraints)
    con_var_list = constraints + variables
    
    var_list, con_list, var_blocks, con_blocks, var_order, con_order = generate_elimination_via_matching(m)
    
    if highlight_matched_edges:
        matching = ComponentMap(zip(con_list, var_list))
        
        #Matrix to plot the matched entries
        matched_matrix = get_matched_entry_submatrix(variables, constraints, matching)
        
        edge_list = get_edges(var_list, con_list, con_var_list)
        edge_colors = []
        for edge in graph.edges:
            if edge in edge_list:
                edge_colors.append('r')
            else:
                edge_colors.append('tab:gray')
       
    if induced_subgraph: 
        #Induced subgraph
        matched_subgraph = igraph.subgraph(var_list, con_list)
        matching_induced_subgraph = get_bipartite_incidence_graph(matched_subgraph.variables, matched_subgraph.constraints)
       
        edge_list = get_edges(var_list, con_list, con_var_list)
        edge_list_subgraph, var_on_subgraph_edges, con_on_subgraph_edges = get_edges_subgraph(matched_subgraph.variables, matched_subgraph.constraints, matching_induced_subgraph, con_var_list)
        edge_colors = []
        for edge in graph.edges:
            if edge in edge_list:
                edge_colors.append('r')
            elif edge in edge_list_subgraph:
                edge_colors.append('y')
            else:
                edge_colors.append('tab:gray')
        
        induced_incidence_matrix = get_submatrix_incidence(variables, constraints, var_on_subgraph_edges, con_on_subgraph_edges)
   
    if matching_elim_plot:
        matching = ComponentMap(zip(con_order, var_order))
        
        #Induced subgraph
        matched_subgraph = igraph.subgraph(var_order, con_order)
        matching_induced_subgraph = get_bipartite_incidence_graph(matched_subgraph.variables, matched_subgraph.constraints)
       
        edge_list = get_edges(var_order, con_order, con_var_list)
        edge_list_subgraph, var_on_subgraph_edges, con_on_subgraph_edges = get_edges_subgraph(matched_subgraph.variables, matched_subgraph.constraints, matching_induced_subgraph, con_var_list)
        edge_colors = []
        for edge in graph.edges:
            if edge in edge_list:
                edge_colors.append('r')
            elif edge in edge_list_subgraph:
                edge_colors.append('y')
            else:
                edge_colors.append('tab:gray')
        matched_elim_matrix = get_matched_entry_submatrix(variables, constraints, matching)
        induced_elim_matrix = get_submatrix_incidence(variables, constraints, var_on_subgraph_edges, con_on_subgraph_edges)
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1, 2, 1.2]})

    if AXESOFF:
        ax1.set_axis_off()
    text = r"""
    \[\begin{array}{cc}\mathrm{eq1:} & x_1 + x_2 - x_4 - y_1^2 = 0  \\ \mathrm{eq2:} & x_2 - y_1y_2 = 0 \\ \mathrm{eq3:} & x_1x_2 - x_3 + x_4 = 0 \\ \mathrm{eq4:} & x_1 - x_2 + y_2 = 0 \\ \mathrm{eq5:} & x_1 + x_2 + x_3x_4 = 2 \\ \end{array}\rightarrow\]
    """
    ax1.set_xlim((0, 1))
    ax1.set_ylim((0, 1))
    ax1.text(-0.1, 0.18, text, fontsize=FONTSIZE*1.1)

    ybuffer = 0.04
    var_rshift = 8
    ax2.set_xlim((-8, 2.5+var_rshift))
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
        edge_color=edge_colors
    )

    label_layout = {}
    node_labels = {}
    var_offset = 3
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
    clabel_offset = 4.5
    ax2.text(
        con_x-clabel_offset,
        nodelabel_y,
        r"$\mathrm{Constraints}$",
        fontsize=FONTSIZE*0.9,
    )
    ax2.text(13, -0.1, r"$\leftrightarrow$", fontsize=FONTSIZE*1.1)

    if AXESOFF:
        ax2.set_axis_off()

    # TODO: axis tick labels?
    ax3.tick_params(length=0)
    ax3.spy(imat, markersize=40)  
    if induced_subgraph:
        ax3.spy(induced_incidence_matrix, markersize = 40, color = 'y')
    if highlight_matched_edges:
        ax3.spy(matched_matrix, markersize = 40, color = 'r')
    if matching_elim_plot:
        ax3.spy(induced_elim_matrix, markersize = 40, color = 'y')
        ax3.spy(matched_elim_matrix, markersize = 40, color = 'r')
       
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
    fig.set_size_inches(4.5*w*factor, 3*h*factor)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.6, hspace=0)
    fig.savefig(file_name, transparent=TRANSPARENT)
    

if __name__ == "__main__":
    filename = 'Bipartite_incidence_graph.svg'
    main(file_name = filename)
    
    filename = 'Matched_edges.svg'
    main(highlight_matched_edges=True, file_name = filename)
    
    filename = 'Matching_induced_subgraph.svg'
    main(highlight_matched_edges=True, induced_subgraph = True, file_name = filename)
    
    filename = "matching_elim.svg"
    main(matching_elim_plot = True, file_name = filename)
    
    
