# variable-elimination
Repository for collaborative exploratory work on variable elimination in NLPs

# Copyright
The copyright header in `header.txt` must be added to every source (Python)
file in this repository.

# Installation
This repository is structured as a small Python package to facilitate code
organization and testing. It can be installed with:
```bash
git clone https://github.com/Robbybp/variable-elimination.git
cd variable-elimination
pip install -r requirements.txt
python setup.py develop
```
Then functionality can be imported in Python:
```python
from var_elim.distill import create_instance
model = create_instance()
```

# Dependencies

## Python dependencies

This repository was developed and tested using Python 3.11.5.
See `requirements.txt` for a list of Python dependencies. For stability, we pin to
specific version of our dependencies, e.g. Pyomo and IDAES.
Non-PyPI dependencies are:
- [`nmpc_examples`](https://github.com/robbybp/nmpc_examples)
- [`pselib`](https://github.com/robbybp/pselib)

When the main branch of this repository requires a specific branch of some dependency
(e.g. the Pyomo main branch, rather than the latest release), an issue should be opened.

## Non-Python dependencies
Results are generated with Ipopt 3.14.17, using the Pyomo-CyIpopt interface. This
interface depends on the following non-Python dependencies:

1. `libpynumero_ASL`. This is a library that allows us to make calls to the AMPL
solver library (ASL) from Python. It can be installed on most systems with the
following commands:
```bash
pip install pyomo
pyomo build-extensions
```

2. [Ipopt 3.14](https://github.com/coin-or/ipopt), with linear solver MA27.
See the following instructions to compile an Ipopt-compatible HSL library
(which contains MA27): https://github.com/coin-or-tools/ThirdParty-HSL.

# Reproducing our results

This repository was used to generate results for the paper:
```bibtex
@misc{naik2025aggregation,
      title={Variable aggregation for nonlinear optimization problems},
      author={Sakshi Naik and Lorenz Biegler and Russell Bent and Robert Parker},
      year={2025},
      eprint={2502.13869},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2502.13869},
}
```
This paper presents three kinds of results:
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

### Producing results in parallel on HPC
The results can be time-consuming to reproduce, so we typically run them in parallel
on multiple-node/core HPC systems. This repository includes scripts to write command
lines that can be run in parallel with utilities like Slurm and GNU Parallel.

## Reproducing structural results
Structural results are produced by the `analyze_structure.py` script:
```bash
python analyze_structure.py
```
Results are written to `RESULTS_DIR/structure.csv`.
Display these results as a Latex table, similar to that displayed in the paper, with:
```bash
python write_latex_table.py `RESULTS_DIR/structure.csv`
```

To write structure-analysis commands that can be run in parallel:
```bash
python write_command_lines.py structure --results-dir=RESULTS_DIR
```

To run these commands in parallel on multiple cores:
```bash
parallel -a structure-commands.txt
```

To collect the results into a single file:
```bash
python collect_results.py structure --results-dir=RESULTS_DIR
```

## Reproducing numerical results

## Reproducing convergence results
