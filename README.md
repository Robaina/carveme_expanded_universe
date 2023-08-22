[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![PyPI version](https://badge.fury.io/py/carveme.svg)](https://badge.fury.io/py/carveme) [![Documentation Status](http://readthedocs.org/projects/carveme/badge/?version=latest)](http://carveme.readthedocs.io/en/latest/?badge=latest)

![CarveMe](logo_300px.png)

# THIS FORK:

Mainly aims at two things:

1) Modifying the universal model curation pipeline so it includes a larger number of reactions. Essentially by i) not removing reactions that are located in eukaryotic compartments, since many prokaryotes contain some of these reactions (i.e., cobalamin synthesis) and ii) modifying the curation pipeline so futile energy cycles are not checked in the universal model, but rather in each individual model. This is to avoid removing reactions that are part of futile cycles only in the universal model but not in the individual models (because they do not co-occur).

2) Adding new functions to the curation pipeline.

## Genome-scale metabolic model reconstruction with CarveMe

CarveMe is a python-based tool for genome-scale metabolic model reconstruction.

### Documentation

For more details please check: http://carveme.readthedocs.io/
 
### Installation

Can be easily installed using **pip**:
```
pip install carveme
```

Additionally, you must install:
   - diamond (conda install -c bioconda diamond)
   - IBM CPLEX or Gurobi solver (full versions with respective academic licenses)

### Credits and License

For citation purposes please refer to our paper:

D. Machado et al, "Fast automated reconstruction of genome-scale metabolic models for microbial species and communities", Nucleic Acids Research, gky537, 2018. doi: https://doi.org/10.1093/nar/gky537

Developed at:
- European Molecular Biology Laboratory (2017-2019).
- Norwegian University of Science and Technology (Since 2020)

CarveMe is now an official service of ELIXIR Norway.

![ELIXIR](elixir.png)
