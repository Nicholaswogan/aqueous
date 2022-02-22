# Aqueous

`aqueous` computes Gibbs energies of formation of aqueous chemical species, and computes chemical equilibrium of aqueous solutions.

## Installation

You need a Fortran compiler (`gfortran>=9.30`, [install instructions here](https://fortran-lang.org/learn/os_setup/install_gfortran)) and C compiler (e.g. install with `conda install -c conda-forge clang`)

Create a `conda` environment with all dependencies

```sh
conda create -n aqueous_env -c conda-forge python numpy numba scikit-build cython nlopt
```

Clone or download this Gitub repository: 

```sh
git clone --depth=1 https://github.com/Nicholaswogan/aqueous.git
```

Navigate to the root directory with a terminal, activate your new `conda` environment, then install:

```sh
conda activate aqueous_env
python -m pip install --no-deps --no-build-isolation .
```

## Usage

**Compute Gibbs energy of formation:**

```python
from aqueous import gibbs_energy

T = 298 # K
P = 1 # bar
G = gibbs_energy("NH3,aq", T, P)
# Gibbs energy in J/mol
```

All Gibbs energies are computed using the [SUPCRT database](https://www.sciencedirect.com/science/article/pii/009830049290029Q).

**Compute equilibrium composition:**

```python
import numpy as np
from aqueous import AqueousSolution

species = ['H+','OH-'] # species in aqueous solution
s = AqueousSolution(species) # make aqueous solution object containing species

m = np.array([1.0e-12, 1.0e-12]) # initial molality of each species (mol species/kg of water)
T = 298 # K
P = 1 # bar
s.equilibrate(m, T, P)
# `m` now contains equilibrium molality of each species
```

The `equilibrate` method finds equilibrium using Gibbs energy minimization. Minimization is done with [NLopt](https://nlopt.readthedocs.io/en/latest/). There are several optimizer settings that can be tweaked:

```python
s.algorithm = 'LD_MMA' # Default is NLopt's LD_MMA. This seems to work best.
s.conserv_tol = 1.0e-9 # relative tolerance for conservation of atoms/charge. Default is 1.0e-9, which works well.
s.ftol = 1.0e-17 # Stopping tolerance for NLopt. Default is 1.0e-17, which works well.
s.xtol = 1.0e-6 # Another stopping tolerance for NLopt. Default is 1.0e-6, which works well.
s.maxtime = 1.0 # Maximum time (seconds) permitted for optimization. Default is 1 second.
```


