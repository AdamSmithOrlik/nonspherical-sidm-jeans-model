from .classes import profile, CDM_profile, isothermal_profile
from . import tools
from . import spherical as sphmodel
from . import nonspherical

# Generate profile object from inputs


# Spherical isothermal Jeans model
def spherical(r1, *outer_halo_params, Phi_b=None, **kwargs):

    # CDM profile
    if r1 == 0:
        outer_halo = CDM_profile(*outer_halo_params, q0=1, Phi_b=Phi_b, **kwargs)
        return profile(outer=outer_halo)

    # Spherical Jeans profile
    else:

        # CDM outer profile (spherical)
        outer_halo = CDM_profile(*outer_halo_params, q0=1, Phi_b=Phi_b, **kwargs)

        # Compute spherically averaged potential from enclosed mass profile (computed within outer_halo)
        Phi_b_sph = tools.compute_Phi_b_spherical(outer_halo.M_b, 1e-6 * outer_halo.r200, outer_halo.r200)

        # Spherical Jeans model (with spherically averaged potential)
        # Matched onto spherical outer halo
        inner_halo, success = sphmodel.relaxation(r1, outer_halo, Phi_b=Phi_b_sph, **kwargs)

        # Matching was successful
        if success:
            halo = profile(inner=inner_halo, outer=outer_halo)
            return halo

        # Unsuccessful matching
        else:
            return None


# Squashed Jeans model
def squashed(r1, *outer_halo_params, q0=1, Phi_b=None, q_mode="smooth", **kwargs):

    # CDM profile
    if r1 == 0:
        return cdm(*outer_halo_params, q0=q0, Phi_b=Phi_b, **kwargs)

    # Squashed Jeans profile
    # Start with spherical profiles
    else:

        # CDM outer profile (spherical)
        outer_halo = CDM_profile(*outer_halo_params, q0=1, Phi_b=Phi_b, **kwargs)

        # Compute spherically averaged potential from enclosed mass profile (computed within outer_halo)
        Phi_b_sph = tools.compute_Phi_b_spherical(outer_halo.M_b, 1e-6 * outer_halo.r200, outer_halo.r200)

        # Spherical Jeans model (with spherically averaged potential)
        # Matched onto spherical outer halo
        inner_halo, success = sphmodel.relaxation(r1, outer_halo, Phi_b=Phi_b_sph, **kwargs)

        # Matching was successful
        if success:

            if q_mode == "uniform":

                # Apply uniform squashing by q0
                halo = profile(inner=inner_halo, outer=outer_halo, q=q0)
                return halo

            elif q_mode == "smooth":

                # Calculate q(r_sph) from spherical Jeans model profile
                sph_halo = profile(inner=inner_halo, outer=outer_halo)
                q_eff = tools.compute_q_eff(sph_halo, q0, **kwargs)

                halo = profile(inner=inner_halo, outer=outer_halo, q=q_eff)
                return halo

            else:
                raise Exception("Unknown q_mode=%s. Must be 'uniform' or 'smooth' only.")

        # Unsuccessful matching
        else:
            return None


# Nonspherical isothermal Jeans model
def isothermal(r1, *outer_halo_params, q0=1, Phi_b=None, **kwargs):

    # CDM profile
    if r1 == 0:
        return cdm(*outer_halo_params, q0=q0, Phi_b=Phi_b, **kwargs)

    # Isothermal nonspherical Jeans profile
    else:

        # CDM outer halo
        outer_halo = CDM_profile(*outer_halo_params, q0=q0, Phi_b=Phi_b, **kwargs)

        # Nonspherical Jeans model matched onto nonspherical outer halo
        inner_halo, success = nonspherical.relaxation(r1, outer_halo, Phi_b=outer_halo.Phi_b, **kwargs)

        # Matching was successful
        if success:
            return profile(inner=inner_halo, outer=outer_halo)

        # Unsuccessful matching
        else:
            return None


# CDM profile
def cdm(*outer_halo_params, q0=1, Phi_b=None, **kwargs):

    # Spherically symmetric CDM halo
    outer_halo = CDM_profile(*outer_halo_params, Phi_b=Phi_b, **kwargs)

    # Return profile object squashed by q0
    return profile(outer=outer_halo, q=q0)
