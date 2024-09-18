# variable-elimination
Repository for collaborative exploratory work on variable elimination in NLPs

# Copyright
The copyright header in `header.txt` must be added to every source (Python)
file in this repository.

# Installation
This repository is structured as a small Python package to facilitate code
organization and testing. It can be installed with:
```bash
python setup.py develop
```
Then functionality can be imported in Python:
```python
from var_elim.distill import create_instance
model = create_instance()
```

# Dependencies
See `requirements.txt` for a list of dependencies. For stability, we pin to specific
version of our dependencies, e.g. Pyomo and IDAES.
Non-PyPI dependencies are:
- [`nmpc_examples`](https://github.com/robbybp/nmpc_examples)
- [`pselib`](https://github.com/robbybp/pselib)

We anticipate keeping PSELib as a dependency, but would like to remove the NMPC-examples
dependency.
When the main branch of this repository requires a specific branch of some dependency
(e.g. the Pyomo main branch, rather than the latest release), an issue should be opened.
Results are generated with Ipopt 3.14.11.
