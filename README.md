# sidmhalo

**sidmhalo** is a Python package for modeling the structure of dark matter halos using the squashed Jeans model described in [this paper](link). **sidmhalo** can model SIDM and CDM halos in one or two dimensions for the spherical and nonspherical cases. The Jeans model code computes the density, rotation curve, shape, potentials, and so on, for a dark matter halo for a given baryon potential, mass, concentration, adiabatic contraction perscription, and outer halo shape profile.

---

## Features

- Nonspherical Jeans equation solver for SIDM halos
- Supports spherical, axisymmetric, and general baryon potentials
- Easy integration with simulation or observational data
- Includes tools for mass profiles, density fits, and rotation curves

---

## Installation

Clone the repository and install in editable (development) mode:

```bash
git clone https://github.com/AdamSmithOrlik/nonspherical-sidm-jeans-model.git
cd nonspherical-sidm-jeans-model
pip install -e .

```

## Examples

See the examples forlder for jupyter notebooks that show how to use the sidmhalo package.
