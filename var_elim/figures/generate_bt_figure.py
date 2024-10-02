import pyomo.environ as pyo
import scipy.sparse as sps
import random
from pyomo.contrib.incidence_analysis.triangularize import block_triangularize
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import get_structural_incidence_matrix
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


random.seed(100)


plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 18


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


def get_imat():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(range(16))

    # Block of 3
    def con1_rule(m, i):
        if i >= 1:
            return m.x[i] == m.x[i-1]
        else:
            return m.x[0] == m.x[2]
    m.con1 = pyo.Constraint(range(3), rule=con1_rule)

    # Block of 5
    def con3_rule(m, i):
        idep = random.randint(0, 7)
        imod = i % 2
        if i >= 1:
            if imod:
                return m.x[8+i] == m.x[7+i] + m.x[6+i] + m.x[idep]
            else:
                return m.x[12] + m.x[8+i] == m.x[7+i] + m.x[6+i] + m.x[idep]
        else:
            return (
                sum(m.x[j] for j in range(8, 12))
                == m.x[12] + m.x[idep]
            )
    m.con3 = pyo.Constraint(range(5), rule=con3_rule)

    def con4_rule(m, i):
        idep = random.randint(0, 12)
        if i >= 1:
            if not i % 3:
                return m.x[15] + m.x[13+i] == m.x[12+i] + m.x[idep]
            else:
                return m.x[13+i] == m.x[12+i] + m.x[idep]
        else:
            return m.x[13] == m.x[15] + m.x[idep]
    m.con4 = pyo.Constraint(range(3), rule=con4_rule)

    # Block of 5
    def con2_rule(m, i):
        idep = random.randint(0, 2)
        if i >= 1:
            if not i % 3:
                return m.x[7] - m.x[3+i] == m.x[2+i] + m.x[idep]
            else:
                return m.x[3+i] == m.x[2+i] + m.x[idep]
        else:  # i == 0
            return (
                sum(m.x[j] for j in range(3, 7))
                == m.x[7] + m.x[idep]
            )
    m.con2 = pyo.Constraint(range(5), rule=con2_rule)

    igraph = IncidenceGraphInterface(m)
    vblocks, cblocks = igraph.block_triangularize()

    print(vblocks)
    print(cblocks)

    vorder = sum(vblocks, [])
    corder = sum(cblocks, [])
    imat = get_structural_incidence_matrix(vorder, corder)

    return imat


def main():
    imat = get_imat()

    fig, (ax1, ax2) = plt.subplots(1, 2)

    text = r"""$\begin{array}{c}\\ B_1 \\ \vdots \\ B_{n_b} \\ \end{array}\hspace{-0.2cm}\begin{array}{c} \begin{array}{ccc} A_1 & \cdots & A_{n_b} \\ \end{array} \\ \left[\begin{array}{ccc} D_1 & & \\ ** & \ddots & \\ ** & * & D_{n_b} \\ \end{array}\right] \end{array} \begin{array}{c} \\ \\ = \\ \\ \end{array}$
    """
    ax1.text(-0.05, 0.1, text)
    ax1.set_axis_off()

    ax2.spy(imat, markersize=5)
    ax2.tick_params(length=0)
    tick_locs = [0, 5, 10, 15]
    tick_vals = [r"$\mathrm{0}$", r"$\mathrm{5}$", r"$\mathrm{10}$", r"$\mathrm{15}$"]
    ax2.xaxis.set_ticks(tick_locs, tick_vals)
    ax2.yaxis.set_ticks(tick_locs, tick_vals)

    rect = get_rectangle_around_coords((0, 0), (2, 2))
    ax2.add_patch(rect)

    rect = get_rectangle_around_coords((3, 3), (7, 7))
    ax2.add_patch(rect)

    rect = get_rectangle_around_coords((8, 8), (12, 12))
    ax2.add_patch(rect)

    rect = get_rectangle_around_coords((13, 13), (15, 15))
    ax2.add_patch(rect)

    w, h = fig.get_size_inches()
    fig.set_size_inches(w, h*0.5)

    fig.savefig("bkg_bt.pdf", transparent=True)


if __name__ == "__main__":
    main()
