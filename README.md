# jeans

**jeans** is a Python package for modeling the structure of dark matter halos using the squashed Jeans model described in [this paper](link). It supports both self-interacting dark matter (SIDM) and cold dark matter (CDM) halos, in spherical and nonspherical geometries. The package provides tools to compute density profiles, rotation curves, shapes, and gravitational potentials for dark matter halos, given a baryon potential, mass, concentration, adiabatic contraction prescription, and outer halo shape profile.

This package is intended for researchers and students working in astrophysics and cosmology, especially those interested in the structure and dynamics of dark matter halos.

## Features

- Nonspherical Jeans equation solver for SIDM halos
- Supports spherical, axisymmetric, and general baryon potentials
- Easy integration with simulation or observational data
- Tools for mass profiles, density fits, and rotation curves
- Example notebooks and scripts included

## Requirements

- Python 3.9+
- numpy
- scipy
- h5py (for EAGLE-50 data)
- matplotlib (for plotting)

All dependencies are listed in `requirements.txt` and will be installed automatically.

## Installation

Clone the repository and install:

```bash
git clone https://github.com/AdamSmithOrlik/nonspherical-sidm-jeans-model.git
cd nonspherical-sidm-jeans-model
pip install .
```

## Usage Example

```python
import jeans

# Example: create a nonspherical SIDM squashed halo with baryons and adiabatic contraction
# rm: matching radius
# M200: mass enclosed within 200 kpc
# c: concentration within 200 kpc
# q0: shape of outer halo at matching radius
# Phi_b: baryon potential as a function of (r, th)
# AC_prescription: prescription for adiabatic contraction

profile = jeans.squashed(rm, M200, c, q0=q0, Phi_b=Phi_b, **{'AC_prescription':'Cautun'})
r = np.logspace(-1, 3, 100)  # kpc
rho_sph = profile.rho_sph_avg(r)

# Plot the density profile
import matplotlib.pyplot as plt
plt.loglog(r, rho_sph)
plt.xlabel('Radius [kpc]')
plt.ylabel('Density [Msun/kpc^3]')
plt.show()
```

See the `examples/` folder for Jupyter notebooks and scripts demonstrating more advanced usage.

To recreate examples involving EAGLE-50 data, download the data from https://zenodo.org/records/16331984 and save it in the `examples/data/` folder.

---

## License

This project is licensed under the terms of the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this package in your research, please cite the associated paper (see [this paper](link)).

---

## Contact

For questions or support, please open an issue on GitHub or contact the maintainer at [asorlik@yorku.ca].

---

## Affiliation

This project is developed and maintained at the Department of Physics and Astronomy, York University, Toronto, Canada.
