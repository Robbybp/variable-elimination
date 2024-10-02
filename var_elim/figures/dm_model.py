import pyomo.environ as pyo
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.incidence_analysis.interface import (
    get_structural_incidence_matrix,
)
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.dae.flatten import flatten_dae_components
from pyomo.util.subsystems import create_subsystem_block
from pyomo.common.collections import ComponentSet

import idaes.core as idaes
from idaes.core.util.model_statistics import degrees_of_freedom

import numpy as np
import matplotlib.pyplot as plt

from var_elim.figures.singular_solid_prop import (
    SingularSolidProperties
)


def make_model():
    m = pyo.ConcreteModel()

    m.component_list = pyo.Set(initialize=["A", "B", "C"])
    n_comp = len(m.component_list)

    # Parameters
    m.area = pyo.Param(initialize=1.0)
    m.void_fraction = pyo.Param(initialize=0.2)
    mw_init_dict = {"A": 0.15969, "B": 0.231533, "C": 0.10196}
    m.mw_comp = pyo.Param(m.component_list, initialize=mw_init_dict)
    m.dens_mass_comp_skeletal = pyo.Param(
        m.component_list,
        initialize={"A": 5250.0, "B": 5000.0, "C": 3987.0},
    )

    # Variables
    # TODO: Should I use units here?
    m.material_accumulation = pyo.Var(m.component_list, initialize=0.0)
    m.energy_accumulation = pyo.Var(initialize=0.0)

    m.flow_mass_comp = pyo.Var(m.component_list, initialize=1.0)
    m.enth_mass = pyo.Var(initialize=1.0)
    m.enth_mol_comp = pyo.Var(m.component_list, initialize=1.0)
    m.flow_mass = pyo.Var(initialize=1.0)
    m.mass_frac_comp = pyo.Var(m.component_list, initialize=1.0/n_comp)
    m.temperature = pyo.Var(initialize=298.15)
    m.dens_mass_particle = pyo.Var(initialize=3252.0)
    m.dens_mass_skeletal = pyo.Var(initialize=3252.0)
    m.volume = pyo.Var(initialize=1.0)
    m.material_holdup = pyo.Var(m.component_list, initialize=1.0)
    m.energy_holdup = pyo.Var(initialize=1.0)
    m.particle_porosity = pyo.Var(initialize=0.27)
    m.velocity = pyo.Var(initialize=1.0)

    m.flow_mass_in = pyo.Var(initialize=1.0)
    m.mass_frac_comp_in = pyo.Var(m.component_list, initialize=1.0/n_comp)
    m.enth_mass_in = pyo.Var(initialize=1.0)

    m.flow_mass_in.fix()
    m.mass_frac_comp_in[:].fix()
    m.enth_mass_in.fix()

    m.material_holdup[:].fix()
    m.energy_holdup.fix()

    m.velocity.fix()
    m.volume.fix()
    m.particle_porosity.fix()

    # Constraints

    @m.Constraint(m.component_list)
    def material_holdup_calculation(m, j):
        return (
            m.material_holdup[j]
            == m.volume * (1 - m.void_fraction)
            * m.dens_mass_particle * m.mass_frac_comp[j]
        )

    @m.Constraint(m.component_list)
    def material_balance(m, j):
        return (
            m.material_accumulation[j]
            == m.flow_mass_in * m.mass_frac_comp_in[j]
            - m.flow_mass * m.mass_frac_comp[j]
        )

    @m.Constraint()
    def energy_holdup_calculation(m):
        return (
            m.energy_holdup
            == m.volume * (1 - m.void_fraction)
            * m.dens_mass_particle * m.enth_mass
        )

    @m.Constraint()
    def enthalpy_balance(m):
        return (
            m.energy_accumulation
            == m.flow_mass_in * m.enth_mass_in
            - m.flow_mass * m.enth_mass
        )

    @m.Constraint(m.component_list)
    def flow_comp_eqn(m, j):
        return (
            m.flow_mass_comp[j]
            == m.flow_mass * m.mass_frac_comp[j]
        )

    conv = 0.001 # Temperature conversion factor
    enthalpy_shomate_coefs = {
        "A": (
            0.1109362, 32.04714*conv**2, 9.192333*conv**3,
            0.901506*conv**4, 5.433677*conv**-1, 843.1471, 825.5032,
        ),
        "B": (
            0.200832, 1.58e-7*conv**2, 6.661e-8*conv**3,
            9.452e-9*conv**4, 3.186e-8*conv**-1, 1174.135, 1120.894,
        ),
        "C": (
            0.102429, 38.7498*conv**2, 15.9109*conv**3,
            2.628181*conv**4, 3.007551*conv**-1, 171793, 1675.69,
        ),
    }

    @m.Constraint(m.component_list)
    def enthalpy_shomate_eqn(m, j):
        A, B, C, D, E, F, G = enthalpy_shomate_coefs[j]
        T = m.temperature
        return (
            m.enth_mol_comp[j]
            == A*T + B*T**2 + C*T**3 + D*T**4 + E/T - F + G
        )

    @m.Constraint()
    def mixture_enthalpy_eqn(m):
        return (
            m.enth_mass
            == sum(
                m.mass_frac_comp[j] * m.enth_mol_comp[j] / m.mw_comp[j]
                for j in m.component_list
            )
        )

    @m.Constraint()
    def sum_component_eqn(m):
        return sum(m.mass_frac_comp[j] for j in m.component_list) == 1.0

    @m.Constraint()
    def density_particle_constraint(m):
        return (
            m.dens_mass_particle
            == (1 - m.particle_porosity) * m.dens_mass_skeletal
        )

    @m.Constraint()
    def density_skeletal_constraint(m):
        return (
            m.dens_mass_skeletal * sum(
                m.mass_frac_comp[j] / m.dens_mass_comp_skeletal[j]
                for j in m.component_list
            )
            == 1.0
        )

    return m


def make_model_from_idaes():
    m = pyo.ConcreteModel()
    fs_config = {"dynamic": True, "time_units": pyo.units.s}
    m.fs = idaes.FlowsheetBlock(default=fs_config)
    time = m.fs.time
    t0 = time.first()
    m.fs.properties = SingularSolidProperties()
    cv_config = {"property_package": m.fs.properties}
    m.fs.cv = idaes.ControlVolume0DBlock(default=cv_config)
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_geometry()
    m.fs.cv.add_material_balances()
    m.fs.cv.add_energy_balances()
    comp_list = m.fs.properties.component_list
    m.fs.cv.flow_mass_comp = pyo.Var(time, comp_list, initialize=1.0)

    @m.fs.cv.Constraint(time, comp_list)
    def flow_comp_eqn(cv, t, j):
        return (
            cv.flow_mass_comp[t, j]
            == cv.properties_out[t0].get_material_flow_terms(None, j)
        )
    disc = pyo.TransformationFactory("dae.finite_difference")
    disc.apply_to(m, wrt=m.fs.time, nfe=1, scheme="BACKWARD")

    # Fix degrees of freedom
    m.fs.cv.volume.fix(1.0)
    for t in time:
        for var in m.fs.cv.properties_in[t].define_state_vars().values():
            var.fix()
    m.fs.cv.properties_out[:].particle_porosity.fix()

    # Fix initial conditions
    m.fs.cv.material_holdup[t0, ...].fix(1.0)
    m.fs.cv.energy_holdup[t0, ...].fix(1.0)

    #
    # Check that we have zero degrees of freedom
    #
    print(degrees_of_freedom(m))

    #
    # Check model for structural singularity
    #
    igraph = IncidenceGraphInterface(m)
    N = len(igraph.constraints)
    M = len(igraph.variables)
    matching = igraph.maximum_matching()
    print(N, M, len(matching))

    #
    # Check subsystems at each point in time for structural singularity
    #
    scalar_vars, dae_vars = flatten_dae_components(m, time, pyo.Var)
    scalar_cons, dae_cons = flatten_dae_components(m, time, pyo.Constraint)

    singular_subsystems = []
    for t in time:
        constraints = [con[t] for con in dae_cons if t in con and con[t].active]
        variables = [var[t] for var in dae_vars if not var[t].fixed]
        n_con = len(constraints)
        n_var = len(variables)
        matching = igraph.maximum_matching(variables, constraints)
        n_matched = len(matching)

        if n_matched != n_var:
            print("Subsystem at %s is structurally singular" % t)
            singular_subsystems.append((constraints, variables))

    constraints, variables = singular_subsystems[0]
    blk = create_subsystem_block(constraints, variables, include_fixed=True)

    ### This is about where I would like to start.

    #
    # Check system for structural nonsingularity
    #
    igraph = IncidenceGraphInterface(blk)
    matching = igraph.maximum_matching()
    M = len(igraph.constraints)
    N = len(igraph.variables)
    print(M, N, len(matching))

    #
    # Display variables, constraints, and parameters in model so I can copy
    # them later
    #
    print("Variables:")
    for var in blk.component_data_objects(pyo.Var):
        print("  %s, %s, %s" % (var.name, var.fixed, var.value))
    print("Constraints:")
    for con in blk.component_data_objects(pyo.Constraint):
        print("  %s" % con.name)
    print("Parameters:")
    for param in blk.component_data_objects(pyo.Objective):
        print("  %s" % param.name)

    return m


if __name__ == "__main__":
    m = make_model()
    igraph = IncidenceGraphInterface(m)
    N = len(igraph.variables)
    M = len(igraph.constraints)
    matching = igraph.maximum_matching()
    print(M, N, len(matching))
