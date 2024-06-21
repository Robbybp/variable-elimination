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

from pselib.testproblem import PseTestProblemBase
from var_elim.models.distillation.distill import create_instance as create_distill
from var_elim.models.gas_pipelines.gas_network_model import make_dynamic_model as create_pipeline


#
# Sample test problem implementation
#
# Test problems should inherit from PseTestProblemBase.
# This is not strictly necessary, but good practice, as it will make sure that
# uid and create_instance are implemented
class DistillationTestProblem(PseTestProblemBase):

    # Every test problem should be given a unique ID
    uid = "DISTILLATION"

    # Implementing an __init__ method is not necessary. The purpose of this method
    # is to allow customization of "meta-parameters" (those that we will not use in
    # the parameter sweep)
    def __init__(self, nfe=300, horizon=300):
        self._nfe = nfe
        self._horizon = horizon

        self._parameters = ["vol", "x_Feed"]
        _parameter_range_list = [(1.01, 10.0), (0.01, 0.99)]
        self._parameter_ranges = dict(zip(self._parameters, _parameter_range_list))

    # For our parameter sweeps, we need to implement the parameters property. This
    # should a list of 2 strings or ComponentUIDs corresponding to fixed variables
    # or mutable parameters on the model.
    @property
    def parameters(self):
        return self._parameters

    # We also need to implement the parameter_ranges property. This is a dict
    # mapping parameters (strings or CUIDs) to (lower-bound, upper-bound) tuples.
    # These are used to sample our inputs for the parameter sweep.
    @property
    def parameter_ranges(self):
        return self._parameter_ranges

    # We need a create_instance method to create a model with the correct
    # meta-parameters. This should basically just me a wrapper around a
    # pre-existing model-construction function we already have.
    #
    # This method does not need to set parameters that we will vary in the sweep.
    # Those are handled separately.
    def create_instance(self):
        return create_distill(
            horizon=self._horizon,
            nfe=self._nfe,
        )

class PipelineTestProblem(PseTestProblemBase):
    uid = "GASPIPELINE"

    def __init__(self, nxfe = 4, ntfe = 24, horizon = 24.0):
        self._nxfe = nxfe
        self._ntfe = ntfe
        self._horizon = horizon
    
        self._parameters = ["fs.nodes[0].state[*].temperature", "fs.nodes[0].state[*].pressure"]
        _parameter_range_list = [(273, 313), (40, 120)]
        self._parameter_ranges = dict(zip(self._parameters, _parameter_range_list))

    @property
    def parameters(self):
        return self._parameters

    @property
    def parameter_ranges(self):
        return self._parameter_ranges

    def create_instance(self):
        return create_pipeline(
                nxfe = self._nxfe,
                ntfe = self._ntfe,
                horizon = self._horizon
        )
