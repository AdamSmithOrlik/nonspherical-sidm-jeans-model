import numpy as np

from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import fsolve

from sidmhalo.definitions import GN

################
# NFW profiles #
################


def rho_NFW(*params, mass_concentration=False):

    if mass_concentration == True:
        # Assume inputs are M200,c200
        M200, c200 = params
        rho_s, rs, r200 = mass_concentration_to_NFW_parameters(M200, c200)
    else:
        # Assume inputs are rho_s, rs
        rho_s, rs = params

    def rho_function(r):
        x = r / rs
        return rho_s / x / (1 + x) ** 2

    return rho_function


def M_NFW(*params, mass_concentration=False):

    if mass_concentration == True:
        # Assume inputs are M200,c200
        M200, c200 = params
        rho_s, rs, r200 = mass_concentration_to_NFW_parameters(M200, c200)
    else:
        # Assume inputs are rho_s, rs
        rho_s, rs = params

    def M_function(r):
        x = r / rs
        if x > 1e-5:
            return 4 * np.pi * rho_s * rs**3 * (np.log(1 + x) - x / (1 + x))
        else:
            return 2 * np.pi * rho_s * rs**3 * x**2

    return np.vectorize(M_function)


def NFW_profiles(*params, **kwargs):

    return rho_NFW(*params, **kwargs), M_NFW(*params, **kwargs)


# f_baryon=0.156352 (Cautun et al value)

#########################################
# Adiabatically contracted NFW profiles #
#########################################


def AC_profiles(
    M200, c200, M_baryon, AC_prescription="Cautun", Gnedin_params=(1.6, 0.8)
):

    rho_s, rs, r200 = mass_concentration_to_NFW_parameters(M200, c200)

    # Use cosmological baryon fraction, not estimated value
    # This value matches one used by Cautun et al.
    f_b = 0.156352

    # This is the baryon fraction estimated directy using input baryon profile
    # f_b = M_baryon(r200)/(M200 + M_baryon(r200))

    # range of r values considered
    rmin = 1e-10 * r200
    rmax = 1e4 * r200
    num_points = 1000

    # DM profile without AC
    M_CDM = M_NFW(rho_s, rs)

    # AC prescription following Cautun et al [1911.04557]
    if AC_prescription == "Cautun":

        # Create table of r values
        r_list = np.logspace(np.log10(rmin), np.log10(rmax), num=num_points)
        M_values = np.array(
            [
                M_CDM(r)
                * (
                    0.45
                    + 0.38 * ((1 - f_b) / f_b * M_baryon(r) / M_CDM(r) + 1.16) ** 0.53
                )
                for r in r_list
            ]
        )

        M_values = np.append(0, M_values)
        r_list = np.append(0, r_list)

        M_function = InterpolatedUnivariateSpline(r_list, M_values)

    # AC prescription following Gnedin et al [1108.5736]
    elif AC_prescription == "Gnedin":

        A0, w = Gnedin_params

        # Define orbit-averaged radius bar(r)/r0 = A0*(r/r0)**w
        def bar(r):
            r0 = 0.03 * r200
            return A0 * r0 * (r / r0) ** w

        # Mass within ri -> mass within rf
        # Use iterative approach to obtain original radius ri
        def find_ri(rf, rtol=1e-6, weight=0.5, max_iter=1000):

            # Initialize
            ri = rf
            ri_old, ri_new = 0, 0
            iter_num = 0

            while (abs(1 - ri_old / ri) > rtol) & (iter_num < max_iter):

                ri_old = ri
                ri_new = rf * (1 - f_b) * (1 + M_baryon(bar(rf)) / M_CDM(bar(ri)))
                ri = weight * ri_new + (1 - weight) * ri

                # Check that ri never becomes negative
                if ri < 0:
                    weight = 0.5 * weight
                    ri = ri_old
                else:
                    weight = 0.5

            if iter_num == max_iter:
                raise Exception("AC error: finding r_i never converged within 'rtol'")

            return ri, M_CDM(ri)

        # Create table of (ri, M_CDM(ri)) values
        r_list = np.logspace(np.log10(rmin), np.log10(rmax), num=num_points)
        table = np.array([find_ri(rf) for rf in r_list])

        M_values = np.append(0, table[:, 1])
        r_list = np.append(0, r_list)

        M_function = InterpolatedUnivariateSpline(r_list, M_values)

    else:
        raise Exception("AC_prescription=" + AC_prescription + " not found.")

    # Compute density rho(r) numerically from M(r):
    # Extrapolate beyond rmin and rmax using power law
    dlogM_dlogr = np.gradient(np.log(M_values[1:]), np.log(r_list[1:]))
    log_rho = np.log(M_values[1:] / (4 * np.pi * r_list[1:] ** 3) * dlogM_dlogr)
    log_rho_function = InterpolatedUnivariateSpline(
        np.log(r_list[1:]), log_rho, k=1, ext=0
    )

    def rho_function(r):
        return np.exp(log_rho_function(np.log(r)))

    return rho_function, M_function


###################
# Other functions #
###################


def mass_concentration_to_NFW_parameters(
    Mvir, c, h=0.7, del_c=200, Omega_m=0.3, Omega_Lambda=0.7, z=0
):

    # Constants
    H0 = h * 100 * 1e-3  # km/kpc/s
    rho_crit = (3 * H0**2) / (8 * np.pi * GN) * (Omega_m * (1 + z) ** 3 + Omega_Lambda)

    # Overdensity
    over_density = (del_c / 3) * ((c**3) / (np.log(1 + c) - (c / (1 + c))))

    # Evaluate parameters
    rvir = np.cbrt(3 * Mvir / (del_c * 4 * np.pi * rho_crit))
    rs = rvir / c
    rhos = rho_crit * over_density

    return rhos, rs, rvir


def NFW_parameters_to_mass_concentration(
    rhos, rs, h=0.7, del_c=200, Omega_m=0.3, Omega_Lambda=0.7, z=0
):

    # Constants
    H0 = h * 100 * 1e-3  # km/kpc/s
    rho_crit = (3 * H0**2) / (8 * np.pi * GN) * (Omega_m * (1 + z) ** 3 + Omega_Lambda)

    # Solve for c
    g = lambda c: c**3 / (np.log(1 + c) - c / (1 + c))
    eqn = lambda c: 3 * rhos / (rho_crit * del_c) - g(c)
    c = fsolve(eqn, 10)[0]

    # Solve for rvir
    rvir = rs * c

    # Sove for Mvir
    Mvir = del_c * rho_crit * 4 * np.pi / 3 * rvir**3

    return Mvir, c, rvir
