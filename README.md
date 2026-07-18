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
We use the IPOPT solver (version 3.14.17 with linear solver MA27)
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
>
> The results for our paper were generated with Python 3.11.5.

# Reproducing our results

Our paper presents three kinds of results:
1. Structural results describe the size and shape of models before loading any
variable values.
2. Numerical results, in our case, are solution objective values, solve times,
and solve time breakdowns.
3. Convergence results are the success or failure of many instances when performing
a sweep over model parameter values.

Scripts to run analyses, summarize results, and generate plots and tables are
located in the `var_elim/scripts` subdirectory of this repository.
You can view the command line interface for each script with:
```bash
python var_elim/scripts/myscript.py --help
```

## Using the `reproduce.py` script

For convenience, we have provided a Python script, `reproduce.py`, that allows
each set of results to be produced with a single command.
This script should be run from the root of this repository.
To make sure this script is working properly, run the "smoke tests"
with:
```
python reproduce.py structure --smoke
python reproduce.py solvetime --smoke
python reproduce.py convergence --smoke
```
By default, this script saves results in the `runs/DATE` subdirectory, where
`DATE` is today's date in YYYYMMDD format. After running the smoke tests,
we should have the following files:
```
runs
└── DATE
    ├── images
    │   ├── distill-matching-sweep-convergence.pdf
    │   ├── fraction-elim.pdf
    │   ├── fraction-solvetime.pdf
    │   └── nnz-per-con.pdf
    └── results
        ├── solvetime-distill-matching.csv
        ├── solvetime-distill-matching.txt
        ├── structure-distill-matching.txt
        ├── structure-distill.csv
        ├── structure-distill.txt
        └── sweep
            ├── distill-matching-sweep-summary.csv
            └── distill-matching-sweep.csv
```
The full results can be produced with:
```
python reproduce.py structure
python reproduce.py solvetime
python reproduce.py convergence
```
The convergence results use 11 samples per parameter by default.
To run an abridged parameter sweep, use e.g., `--nsamples=2`.
After running these commands, we should have the following files:
```
TODO
```

## Producing results in parallel on HPC
The results can be time-consuming to reproduce, so we typically run them in parallel
on multiple-node/core HPC systems. We have provided a makefile with the commands
used to run these results in parallel.

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
