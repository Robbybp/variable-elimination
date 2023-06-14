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
#  ___________________________________________________________________________
#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
This file contains utilities for constructing a flowsheet for Kai's 5 demands,
3 compressors, 1 supply gas pipeline network.
"""
import pyomo.environ as pyo
import idaes.core as idaes
import pyomo.network as network
from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.models_extra.gas_distribution.unit_models.pipeline import GasPipeline
from idaes.models_extra.gas_distribution.unit_models.compressor import (
    IsothermalCompressor as Compressor,
)
from idaes.models_extra.gas_distribution.unit_models.node import PipelineNode

from idaes.core.util.model_statistics import degrees_of_freedom

from nmpc_examples.nmpc.dynamic_data.series_data import (TimeSeriesData)
from nmpc_examples.nmpc.model_helper import DynamicModelHelper
from idaes.apps.nmpc.dynamic_data import (interval_data_from_time_series,
                                          load_inputs_into_model
                                          )
def make_model(
    dynamic=True,
    nxfe=2,
    space_method="dae.finite_difference",
    space_scheme="FORWARD",
    ntfe=40,
    horizon=20.0,
    time_method="dae.finite_difference",
    time_scheme="BACKWARD",
    ):
    
    m = pyo.ConcreteModel()
    default = {"dynamic": dynamic}
    if dynamic:
        default["time_set"] = [0.0, horizon]
        default["time_units"] = pyo.units.hr
    
    m.fs = idaes.FlowsheetBlock(**default)
    m.fs.properties = NaturalGasParameterBlock()


    

    node_configs = [
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 0,
            "n_outlet_pipelines": 1,
            "n_supplies": 1,
            "n_demands": 0,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 2,
            "n_supplies": 0,
            "n_demands": 0,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 2,
            "n_supplies": 0,
            "n_demands": 1,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 0,
            "n_supplies": 0,
            "n_demands": 1,
        },
        
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 0,
            "n_supplies": 0,
            "n_demands": 1,
        },
        
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 1,
            "n_supplies": 0,
            "n_demands": 1,
        },
        
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 0,
            "n_supplies": 0,
            "n_demands": 1,
        }
    ]
    m.fs.node_set = pyo.Set(initialize=list(range(len(node_configs))))
    node_configs = {i: config for i, config in enumerate(node_configs)}
    m.fs.nodes = PipelineNode(m.fs.node_set, initialize=node_configs)
   
    pipeline_config = {
        "property_package": m.fs.properties,
        "finite_elements": nxfe,
        "transformation_method": space_method,
        "transformation_scheme": space_scheme,
        "has_holdup": True,
    }
    m.fs.pipeline_set = pyo.Set(initialize=range(6))
    m.fs.pipeline = GasPipeline(m.fs.pipeline_set, **pipeline_config)
    
    compressor_config = {"property_package": m.fs.properties}
    m.fs.compressor_set = pyo.Set(initialize=range(3))
    m.fs.compressor = Compressor(
        m.fs.compressor_set, **compressor_config
    )
    

    # Connect compressors to pipelines

    m._compressor_to_pipeline_1 = network.Arc(
        ports=(m.fs.compressor[0].outlet_port,
        m.fs.pipeline[0].inlet_port)
    )
    
    m._compressor_to_pipeline_2 = network.Arc(
        ports=(m.fs.compressor[1].outlet_port,
        m.fs.pipeline[1].inlet_port)
    )
    
    m._compressor_to_pipeline_3 = network.Arc(
        ports=(m.fs.compressor[2].outlet_port,
        m.fs.pipeline[4].inlet_port)
    )
    
    m.fs.nodes[0].add_pipeline_to_outlet(m.fs.compressor[0])
    
    m.fs.nodes[1].add_pipeline_to_inlet(m.fs.pipeline[0])
    m.fs.nodes[1].add_pipeline_to_outlet(m.fs.compressor[1])
    m.fs.nodes[1].add_pipeline_to_outlet(m.fs.compressor[2])
    
    m.fs.nodes[2].add_pipeline_to_inlet(m.fs.pipeline[1])
    m.fs.nodes[2].add_pipeline_to_outlet(m.fs.pipeline[2])
    m.fs.nodes[2].add_pipeline_to_outlet(m.fs.pipeline[3])
    
    m.fs.nodes[3].add_pipeline_to_inlet(m.fs.pipeline[2])
    m.fs.nodes[4].add_pipeline_to_inlet(m.fs.pipeline[3])
    
    m.fs.nodes[5].add_pipeline_to_inlet(m.fs.pipeline[4])
    m.fs.nodes[5].add_pipeline_to_outlet(m.fs.pipeline[5])
    
    m.fs.nodes[6].add_pipeline_to_inlet(m.fs.pipeline[5])
    
    
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m)
    
    #Set pipeline length
    cv = m.fs.pipeline[:].control_volume
    m.fs.pipeline[:].diameter.fix(0.92*pyo.units.m)
    cv.length.fix(300*pyo.units.km)
    
    
    # Initial conditions:
   
    x0 = m.fs.pipeline[0].control_volume.length_domain.first()
    xf = m.fs.pipeline[0].control_volume.length_domain.last()
    j = next(iter(m.fs.properties.component_list))
    if dynamic:
        t0 = m.fs.time.first()
        for x in m.fs.pipeline[0].control_volume.length_domain:
            # Here I assume that all three pipelines have the same
            # length domain.
            if x != x0:
                m.fs.pipeline[0].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[1].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[2].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[3].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[4].control_volume.pressure[t0, x].fix()
                m.fs.pipeline[5].control_volume.pressure[t0, x].fix()
            if x != xf:
                m.fs.pipeline[0].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[1].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[2].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[3].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[4].control_volume.flow_mass[t0, x].fix()
                m.fs.pipeline[5].control_volume.flow_mass[t0, x].fix()
            
        cv.momentum_balance[t0, xf].deactivate()
        
        disc = pyo.TransformationFactory(time_method)
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme=time_scheme)

            
    # Fix "dynamic inputs." This needs to be done after a potential
    # discretization transformation.

    m.fs.nodes[0].state[:].mole_frac_comp[j].fix()
    m.fs.nodes[0].state[:].temperature.fix()
  
    
    return m 
    
def make_steady_model(
        nxfe=4,
        to_fix=None,
        input_data=None,
        tee=True,
        ):
    if to_fix is None:
        to_fix = []
    if input_data is None:
        input_data = {}
    m = make_model(nxfe=nxfe, dynamic=False)

    for cuid in to_fix:
        var = m.find_component(cuid)
        var[:].fix()

    for cuid, val in input_data.items():
        var = m.find_component(cuid)
        var[:].set_value(val)

    return m

def get_half_demand_profile(demand_cuid=None):
    """
    """
    if demand_cuid is None:
        demand_cuid = pyo.ComponentUID("target_demand")
    sample_points = [0.0, 6.0, 12.0,18.0, 24.0]

    #
    # Converting "paper units" (1e4 SCM/hr) to "IDAES units" (kg/hr)
    #
    # TODO: Change this later and actually use unit-ed expressions here.
    # 30 (1e4 SCM)/hr -> kmol/hr
    val0 = pyo.value(30 * 1e4 * 0.72 / 18.0)
    # 50 (1e4 SCM)/hr -> kmol/hr
    val1 = pyo.value(45.0 * 1e4 * 0.72 / 18.0)
    val2 = pyo.value(65.0 * 1e4 * 0.72 / 18.0)

    # This is a time series where each value corresponds to the right
    # end point of a piecewise-constant interval (or the first time
    # point).
    demand_data = TimeSeriesData(
        {
            demand_cuid: [val0, val0, val1, val2, val0],
        },
        sample_points,
    )
    return demand_data

def demand_data_helper():
    demand_nodes = ["fs.nodes[2].demands[0].flow_mol[*]",
                         "fs.nodes[3].demands[0].flow_mol[*]",
                         "fs.nodes[4].demands[0].flow_mol[*]",
                         "fs.nodes[5].demands[0].flow_mol[*]",
                         "fs.nodes[6].demands[0].flow_mol[*]"]
    demand = []
    for d in demand_nodes:
        demand_cuid = pyo.ComponentUID(d)
        demand_seq = get_half_demand_profile(demand_cuid = demand_cuid)
        demand.append(demand_seq)
    return demand

def make_dynamic_model(horizon=24.0, ntfe= 24, eta=0.7):
    demand = demand_data_helper()
    ipopt = pyo.SolverFactory("ipopt")
    nxfe = 4
    demand_nodes = [2,3,4,5,6]
    nominal_steady_inputs = {
        "fs.nodes[0].state[*].pressure": 100.0,
        "fs.compressor[0].boost_pressure[*]": 20.0,
        "fs.compressor[1].boost_pressure[*]": 0.0,
        "fs.compressor[2].boost_pressure[*]": 0.0,
    }
    

    #Steady state model with nominal steady state inputs. 
    m_initial = make_steady_model(
        nxfe=nxfe,
        input_data=nominal_steady_inputs,
        to_fix=list(nominal_steady_inputs.keys()),
    )
    
    #Fix compressor efficiency of the model
    
    m_initial.fs.compressor[:].efficiency = 0.7
    
    for d in demand_nodes:
        m_initial.fs.nodes[d].demands[0].flow_mol.fix(12000)
    print("DOF: %s" % degrees_of_freedom(m_initial))
    
    #Solving a square steady state problem with fixed demands and compressor boost pressures
    ipopt.solve(m_initial, tee=True)

    m_initial_time = m_initial.fs.time
    m_initial_helper = DynamicModelHelper(m_initial, m_initial_time)
    
    #Inital data from the steady state model.
    #Will be used in the future in the dynamic model initialization 
    initial_data = m_initial_helper.get_data_at_time()
    init_scalar_data = m_initial_helper.get_scalar_variable_data()
    
   
    m = make_model(
        nxfe=nxfe,
        dynamic=True,
        horizon=horizon,
        ntfe=ntfe,
    )
    m.fs.compressor[:].efficiency = eta

    time = m.fs.time
    
    m_helper = DynamicModelHelper(m, time)
    m_helper.load_data_at_time(initial_data)
    m_helper.load_scalar_data(init_scalar_data)
    
    m.fs.nodes[0].state[:].pressure.fix(100)
    #Fixing demands? Do I want to add demand penalties to the objective instead?
    for i in range(len(demand)):
        demand_data = (demand[i].get_time_points(), demand[i].get_data())
        demand_data = interval_data_from_time_series(demand_data)
        load_inputs_into_model(m, m.fs.time, demand_data)
    
   
    for d in demand_nodes:
        m.fs.nodes[d].demands[0].flow_mol.fix()
        
    t0 = m.fs.time.first()
    m.fs.compressor[:].boost_pressure[:].unfix()
    m.fs.compressor[:].boost_pressure[t0].fix()
    
    for d in demand_nodes:
        m.fs.nodes[d].state[:].pressure[:].setlb(35)
    
    m.fs.nodes[0].supplies[0].state[:].flow_mol.setlb(1000.0)
    m.fs.nodes[0].supplies[0].state[:].flow_mol.setub(1e+10)
    m.fs.pipeline[:].control_volume.pressure[:,:].setlb(0)
    m.fs.pipeline[:].control_volume.flow_mass[:,:].setlb(0)
    m.fs.compressor[:].power.setub(60000)
    m.fs.compressor[:].power[t0].setub(None)
    
    
    
    #Terminal bounds
    tf = m.fs.time.last()
    [m.fs.pipeline[:].control_volume.flow_mass[tf,x].setlb(100000.0) for x in m.fs.pipeline[0].control_volume.length_domain if x!= m.fs.pipeline[0].control_volume.length_domain.last()]
    [m.fs.pipeline[:].control_volume.pressure[tf,x].setlb(20.0) for x in m.fs.pipeline[0].control_volume.length_domain if x!= m.fs.pipeline[0].control_volume.length_domain.first()]
    
    #----------------------------------------------------------------------------------------
    #Objective
    
    m.obj = pyo.Objective(expr = 0.1*(sum(m.fs.compressor[0].power[t] for t in m.fs.time)
                          + sum(m.fs.compressor[1].power[t] for t in m.fs.time)
                          + sum(m.fs.compressor[2].power[t] for t in m.fs.time)))
                             
                              
    m.total_compression_cost = pyo.value(0.1*(sum(m.fs.compressor[0].power[t] for t in m.fs.time)
                          + sum(m.fs.compressor[1].power[t] for t in m.fs.time)
                          + sum(m.fs.compressor[2].power[t] for t in m.fs.time)))
    
    print("DOF of the dynamic problem:", degrees_of_freedom(m))
    
    return m 

def main():
    m = make_dynamic_model()
    print("DOF: %s" % degrees_of_freedom(m))
    
    ipopt = pyo.SolverFactory('ipopt')
    ipopt.solve(m, tee=True)
if __name__ == "__main__":
    main()