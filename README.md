Variable Elimination
--------------------
This repository started as a place for collaborative exploratory work
on variable elimination in nonlinear optimization problems. It now primarily
serves to host the code necessary to reproduce the results in the paper
"Variable aggregation in nonlinear optimization problems":
```bibtex
@misc{naik2026aggregation,
  title={Variable aggregation for nonlinear optimization problems},
  author={Sakshi Naik and Lorenz Biegler and Russell Bent and Robert Parker},
  year={2026},
  eprint={2502.13869},
  archivePrefix={arXiv},
  primaryClass={math.OC},
  url={https://arxiv.org/abs/2502.13869},
}
```

# Installation
This repository is structured as a small Python package to facilitate code
organization and testing. To install, you must first install our non-Python
dependencies
([PyNumero](https://pyomo.readthedocs.io/en/stable/explanation/solvers/pynumero/index.html)
and [IPOPT](https://github.com/coin-or/ipopt))
then install this repository as a Python package.

## Installing non-Python dependencies

### PyNumero

We use PyNumero to interface Pyomo, CyIpopt, and the AMPL Solver Library (ASL).
PyNumero can be installed by the Pyomo Python package:
```
pip install pyomo
pyomo build-extensions
```
This compiles `libpynumero_ASL` and saves it to an OS-dependent Pyomo directory, e.g.,
`$HOME/.pyomo/lib` on Linux.
The end of the output of `pyomo build-extensions` should look like this:
```
INFO: Finished building Pyomo extensions.
INFO: The following extensions were built:
    [ OK ]  ampl_function_demo
    [FAIL]  appsi
    [ OK ]  cspline_external
    [ OK ]  aslfunctions
    [ OK ]  mcpp
    [ OK ]  pynumero
    [SKIP]  ginac
```
As long as `pynumero` has status `OK`, this step was successful for the purpose
of this repository.

### IPOPT
We use the IPOPT solver (version 3.14.10 with linear solver MA27)
to compare the impact of our aggregation techniques.
IPOPT can be installed by following instructions on the github repository,
[](https://github.com/coin-or/ipopt).
See the following instructions to compile an Ipopt-compatible HSL library
(which contains MA27): https://github.com/coin-or-tools/ThirdParty-HSL.

> [!NOTE]
> Depending on your IPOPT install location, you may have to set some environment
> variables for `libipopt` to be discoverable by the CyIpopt Python package.
> Commonly, `PKG_CONFIG_PATH=INSTALL_DIR/lib/pkgconfig`,
> `LD_LIBRARY_PATH=INSTALL_DIR/lib`, and/or `DYLD_LIBRARY_PATH=INSTALL_DIR/lib`
> must be set.

## Installing this package
Download this repository and navigate into it with:
```
git clone https://github.com/Robbybp/variable-elimination.git
cd variable-elimination
```
The following command, run from the root of this repository, will install this
package and the latest versions of all its Python dependencies:
```
pip install -e .[test]
```
(The `[test]` flag installs `pytest`, which is only necessary to run tests.)
Verify that the installation was successful by running `pytest` from the
repository root.
The functionality of this package can then be imported in Python, e.g.:
```python
from var_elim.models.distillation.distill import create_instance
model = create_instance()
```

> [!NOTE]
> The above command installs the *latest* version of dependencies.
> The exact versions of all dependencies used to produce the results
> in our paper are provided in `requirements.txt`. Install these
> dependencies with
> ```
> pip install -r requirements.txt
> ```
> before or after installing this package.

# Reproducing our results

Our paper presents three kinds of results:
1. Structural results describe the size and shape of models before loading any
variable values.
2. Numerical results, in our case, are solution objective values, solve times,
and solve time breakdowns.
3. Convergence results are the success or failure of many instances when performing
a sweep over model parameter values.

Scripts to produce these results are located in the `var_elim/scripts` subdirectory
of this repository.
You can view the command line interface for each script with:
```bash
python myscript.py --help
```
In particular, the `--results-dir=RESULTS_DIR` option specifies the directory
where results are written, and the `--image-dir=IMAGE_DIR` option specifies the
directory where images are saved (for scripts that plot figures).
These default to `results` and `images` respectively, but it might be useful
to set them to custom values to avoid overwriting previous results
(produced with, e.g., a different Pyomo version).

TL;DR: Reproduce the bulk of our results with the following commands:
```bash
python analyze_structure.py
python write_latex_table.py results/structure.csv
python analyze_solvetime.py
python write_latex_table.py results/solvetime.csv
python run_param_sweep.py --model=distill
python run_param_sweep.py --model=mb-steady
python run_param_sweep.py --model=pipeline
python summarize_sweep_results.py --model=distill
python summarize_sweep_results.py --model=mb-steady
python summarize_sweep_results.py --model=pipeline
```
See below for more details, especially on how to speed this up in an HPC environment.

### Producing results in parallel on HPC
The results can be time-consuming to reproduce, so we typically run them in parallel
on multiple-node/core HPC systems. This repository includes scripts to write command
lines that can be run in parallel with utilities like Slurm and GNU Parallel.

## Reproducing structural results
Structural results are produced by the `analyze_structure.py` script:
```bash
python analyze_structure.py --results-dir=RESULTS_DIR
```
Results are written to `RESULTS_DIR/structure.csv`.
Display these results as a Latex table, similar to that displayed in the paper, with:
```bash
python write_latex_table.py RESULTS_DIR/structure.csv
```

To write structure-analysis commands that can be run in parallel:
```bash
python write_command_lines.py structure --results-dir=RESULTS_DIR
```

To run these commands in parallel with multiple subprocesses:
```bash
# This may require installing GNU Parallel for the `parallel` command
parallel -a structure-commands.txt
```

To collect the results into a single file:
```bash
# Only necessary if we have written many small files in parallel runs!
python collect_results.py structure --results-dir=RESULTS_DIR
```
The results are now combined into `RESULTS_DIR/structure.csv` and can be
displayed as above.
Figures may be produced with:
```bash
python plot_structure_bargraphs.py RESULTS_DIR/structure.csv --image-dir=IMAGE_DIR
```

## Reproducing numerical results
Numerical results are produced using the `analyze_solvetime.py` script:
```bash
python analyze_solvetime.py --results-dir=RESULTS_DIR
```
Results are written to `RESULTS_DIR/solvetime.csv`, and can be displayed
with:
```bash
python write_latex_table.py RESULTS_DIR/solvetime.csv
```

Independent commands for parallel runs can be written with:
```bash
python write_command_lines.py solvetime --results-dir=RESULTS_DIR
```

We typically prefer to run scripts that measure solvetimes on independent,
identical compute nodes rather than using multiple processes on the same
node. We do this with a Slurm batch script and the `sbatch` command.
The exact contents of this batch script depends on your HPC environment.

Collect results into a single file with:
```bash
# Only necessary if we have written many small files in parallel runs!
python collect_results.py solvetime --results-dir=RESULTS_DIR
```

The breakdown of solve time can be plotted with:
```bash
python plot_timing_bargraphs.py RESULTS_DIR/solvetime.csv --image-dir=IMAGE_DIR
```

## Reproducing convergence results
Parameter sweeps are run with the `run_param_sweep.py` script:
```bash
python run_param_sweep.py --results-dir=RESULTS_DIR
```
This writes a CSV file of convergence results for each model-method combination,
e.g., `mb-steady-matching-sweep.csv`, into the `RESULTS_DIR/sweep` subdirectory.

A summary of sweep results may be displayed with:
```bash
python summarize_sweep_results.py --results-dir=RESULTS_DIR --model=MODEL
```
where `MODEL` is one of `mb-steady`, `distill`, or `pipeline`.

Parameter sweep success/failure results may be plotted in a grid with:
```bash
python plot_sweep_results.py SWEEP_CSV --image-dir=IMAGE_DIR
```

We typically run the parameter sweep for each model-method combination
on a different compute node (managed by e.g., Slurm), and use multiprocess
parallelism to run individual parameter samples in parallel within each node
(using GNU Parallel).

Commands can be written with:
```bash
python write_sweep_command_lines.py --results-dir=RESULTS_DIR
```
The commands for each model-method combination, written to
`commands/parallel-sweep-commands.txt`, can then be run on multiple compute
nodes with, for example, the Slurm `sbatch` command.

Parameter sweep results for each model-method combination may be combined
into CSV files, e.g., `mb-steady-matching-sweep.csv` with:
```bash
# Alternatively, we could run each command in this file manually
parallel -a commands/collect-sweep-commands.txt
```

Plots may be generated with:
```bash
# Alternatively, we could run each command in this file manually
parallel -a commands/plot-sweep-commands.txt
```

# Using these methods on your own models
The functionality to identify sets of variables and constraints to eliminate, perform
the elimination in-place on Pyomo models, and analyze the resulting models is located
in the `var_elim/algorithms` and `var_elim/heuristics` subdirectories.
Different methods, e.g., `identify_vars_for_elim_ampl` and `get_degree_two_elimination`,
have different call signatures. For ease of comparison among the different methods,
the `var_elim.elimination_callbacks` module has functions with consistent call signatures
for performing each elimination method in-place on a Pyomo model.

For example, to perform the "linear-matching" elimination method on our distillation
column test model, run the following:
```python
from var_elim.models.distillation.distill import create_instance
from var_elim.elimination_callbacks import matching_elim_callback
model = create_instance()
matching_elim_callback(model)
```

`pyomo.contrib.incidence_analysis` provides some useful tools to inspect the
structure of resulting models:
```python
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
igraph = IncidenceGraphInterface(model, linear_only=False)
nvar = length(igraph.variables)
ncon = length(igraph.constraints)
nnz = igraph.n_edges
```

# License
This code is copyrighted by Los Alamos National Laboratory and distributed under a
BSD 3-clause license. The license and copyright are provided in `LICENSE.md` and
`COPYRIGHT.txt`.
