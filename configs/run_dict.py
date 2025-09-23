run_dictionary = {
    "r1": 10.0,  # kpc. Matching radius for SIDM and CDM halos. If r1=0 CDM halo is returned.
    "M200": 1e12,  # Msun
    "c": 10.0,  # Dimensionless. Concentration parameter. c = r200/rs.
    "q0": 1.0,  # Dimensionless. Initial outer halo shape. If q0=1, spherical halo is assumed.
    "alpha": None,  # Dimensionless. Halo flattening parameter for Einasto profile. If None, NFW profile is used.
    "Phi_b": None,  # (km/s)^2. Baryon potential. Must be a function with signature Phi_b(r, theta), even if spherical.
    "AC_prescription": None,  # Adiabatic contraction prescription. Options: 'Cautun' or 'Gnedin'.
    "Gnedin_params": (
        1.6,
        0.8,
    ),  # Only used if AC_prescription='Gnedin'. (A, w) parameters.
    "save_profile": True,  # If True, saves the profile to a .npz file.
    "save_dir": "profiles/",  # Directory to save the profile .npz file.
    "nonspherical_model": "squashed",  # Options: 'squashed' or 'isothermal'.
    "verbose": True,  # If True, prints progress and warnings.
    "L_list": [0],  # Angular momentum modes to include in the isothermal model.
    "M_list": [0],  # M > 0 not yet implemented.
    "q_mode": "smooth",  # Only used if nonspherical_model='squashed'. Options: 'uniform' or 'smooth'.
}
