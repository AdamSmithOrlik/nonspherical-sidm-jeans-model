from analytic_potentials import Phi_MN

# fmt: off
run_dictionary = {
    "model": "spherical",  # Options: 'spherical', 'cdm', 'squashed' or 'isothermal'. Setting to cdm will override r1 st r1=0. Setting to spherical will override q0 so that q0=1.
    "r1": 10.0,  # kpc. Matching radius for SIDM and CDM halos. If r1=0 CDM halo is returned.
    "M200": 1e12,  # Msun
    "c": 10.0,  # Dimensionless. Concentration parameter. c = r200/rs.
    "q0": 1.0,  # Dimensionless. Initial outer halo shape. If q0=1, spherical halo is assumed.
    "alpha": None,  # Dimensionless. Halo flattening parameter for Einasto profile. If None, NFW profile is used.
    "Phi_b": None,  # (km/s)^2. Baryon potential. Must be a function with signature Phi_b(r, theta), even if spherical.
    "AC_prescription": None,  # Adiabatic contraction prescription. Options: 'Cautun' or 'Gnedin'.
    "Gnedin_params": (1.6, 0.8,),  # Only used if AC_prescription='Gnedin'. (A, w) parameters.
    "save_profile": True,  # If True, saves the profile to a .npz file.
    "save_dir": "data/",  # Relative path to save the profile .npz file.
    "verbose": False,  # If True, prints progress and warnings.
    "L_list": [0],  # Angular momentum modes to include in the isothermal model.
    "M_list": [0],  # M > 0 not yet implemented.
    "q_mode": "smooth",  # Only used if model='squashed'. Options: 'uniform' or 'smooth'.
}

# Baryon potential function
Md = 6e10  # Msun
a = 3.0  # kpc
b = 0.28  # kpc

Phi_b = Phi_MN(Md, a, b)  # Miyamoto-Nagai potential function

# ~~~~~~~~~~~~~~update the run dictionary~~~~~~~~~~~~~~ #
# run_dictionary["Phi_b"] = Phi_b # uncomment to include baryon function

# CDM halo shape function
q0 = 0.8  # Dimensionless. Outer halo shape.
r200 = 200  # kpc. Virial radius.
gamma = 0.2  # Dimensionless. Shape parameter.

# Non-constant CDM halo shape function
def q_cdm(q0, r200=200, gamma=0.3):
    return lambda r: q0 * (r / r200) ** gamma

# ~~~~~~~~~~~~~~update the run dictionary~~~~~~~~~~~~~~ #
# run_dictionary["q0"] = q_cdm(q0, r200, gamma) # uncomment to include non-constant halo shape function
