import pyomo.environ as pyo
from pyomo.contrib.incidence_analysis.visualize import spy_dulmage_mendelsohn
import pselib
from var_elim.models.testproblems import (
    DistillationTestProblem,
)
from var_elim.elimination_callbacks import matching_elim_callback
from var_elim.scripts.analyze_structure import get_structural_results
import matplotlib.pyplot as plt

from pyomo.common.timing import HierarchicalTimer

timer = HierarchicalTimer()

#Get steady state mb model
def mb_steady_constructor():
    model = pselib.get_problem("MBCLC-METHANE-STEADY").create_instance()
    pyo.TransformationFactory("core.scale_model").apply_to(model)
    return model


def plot_incidence_matrix():
    model = mb_steady_constructor()
    
    #Plot no-elim incidence
    spy_dulmage_mendelsohn(model,
                           highlight_coarse=False,
                           highlight_fine=False,)
    plt.title("Incidence Matrix - No elimination")
    plt.show()
    
    #Perform matching based elimination
    get_structural_results(model, matching_elim_callback, htimer = timer)
    
    #Plot no-elim incidence
    spy_dulmage_mendelsohn(model,
                           highlight_coarse=False,
                           highlight_fine=False,)
    plt.title("Incidence Matrix - Linear Matching")
    plt.show()
if __name__ == "__main__":
    plot_incidence_matrix()
    
    