import numpy as np
from inspect import signature
from scipy.integrate import solve_ivp
from findiff import FinDiff

from sidmhalo.definitions import GN, Z, integrate

# Derivitive helper function
def central_derivative_findiff(f, x, dx):
    D = FinDiff(0, dx, 1)
    f_vals = np.array([f(x - dx), f(x), f(x + dx)])
    return D(f_vals)[1]  

# Baryon contributions to rotation curve
def Vsq_baryon(Phi_b,r):
        
    # Baryon potential
    num_Phi_b_variables = len(signature(Phi_b).parameters)
    dPhi_b_dr = np.zeros_like(r)

    r_arr = np.array(r)
    pos = r_arr > 0

    # if num_Phi_b_variables == 1:
    #     dPhi_b_dr[pos] = np.array([derivative(Phi_b, ri, dx=1e-2*ri) for ri in r_arr[pos]])

    # elif num_Phi_b_variables == 2:
    #     theta = np.pi/2
    #     dPhi_b_dr[pos] = np.array([derivative(lambda r: Phi_b(r,theta), ri, dx=1e-2*ri) for ri in r_arr[pos]])

    if num_Phi_b_variables == 1:
        dPhi_b_dr[pos] = np.array([
            central_derivative_findiff(Phi_b, ri, dx=1e-2 * ri)
            for ri in r_arr[pos]
        ])

    elif num_Phi_b_variables == 2:
        theta = np.pi / 2
        dPhi_b_dr[pos] = np.array([
            central_derivative_findiff(lambda r: Phi_b(r, theta), ri, dx=1e-2 * ri)
            for ri in r_arr[pos]
        ])
    else:
        raise Exception('Case with %d Phi_b arguments not supported.' % num_Phi_b_variables)

    Vsq_out = r_arr * dPhi_b_dr

    if np.ndim(r) == 0:
        return float(Vsq_out)
    else:
        return Vsq_out

# Halo contribution to rotation curve for given L,M mode

def Vsq_LM(rho_LM,r,L,M=0):

    if M != 0:
        raise Exception("M != 0 not implemented.")

    # Handle case where r is just a number
    # Code below assumes r is array or list
    if np.ndim(r) == 0:
        return Vsq_LM(rho_LM,[r],L,M=M)[0]
    
    elif np.ndim(r) > 1:
        raise Exception("r must be number or 1D array/list.")
    else:
        pass

    # Make sure r is ordered
    order = np.argsort(r)
    r_arr = np.array(r,dtype='float')[order]
    r_eval = r_arr[r_arr > 0]

    # Prefactor
    prefactor = 4*np.pi*GN / (2*L+1)

    # Initialize output
    Vsq_out = np.zeros_like(r_arr)

    if max(r) > 0:

        # First term: 
        # G(r) = int_0^r dx x^(L+2) rho_LM(x)

        # End points for integration
        rmin, rmax = 0, max(r_eval)

        # integrand = lambda r,y: r**(L+2) * rho_LM(r,L)
        def integrand(r,y): 
            if r > 0:
                return r**(L+2) * rho_LM(r)
            else:
                return 0

        # Calculate integrals
        solution = solve_ivp(integrand,[rmin,rmax],[0],rtol=1e-6,atol=1e-6,t_eval=r_eval)
        G_vals = solution.y[0]

        Vsq_out[r_arr > 0] += prefactor * (L+1) / r_eval**(L+1) * G_vals

        # Second term: only needed if L > 0
        # H(r) = int_r^inf x^(1-L) rho_LM(x) = - F(r) + H0
        # where F(r) = int_rmin^r x^(1-L) rho_LM(x)
        # and H0 = int_rmin^inf x^(1-L) rho_LM(x)
        if L > 0:

            integrand = lambda r,y: r**(1-L) * rho_LM(r)
            rmin, rmax = r_eval[0], r_eval[-1]

            H_vals = np.zeros_like(r_eval)
            integral_list = []
            for i in range(len(r_eval)-1):
                integral = integrate(lambda r: r**(1-L) * rho_LM(r),r_eval[i],r_eval[i+1],rtol=1e-6,atol=1e-6)
                integral_list.append(integral)

            for i in range(len(H_vals)):
                H_vals[i] = np.sum(integral_list[i:])

            Vsq_out[r_arr > 0] +=  - prefactor * L * r_eval**L * H_vals 

            # Finally add extra piece int_rmax^infty dx x^(1-L) rho_LM(x)
            max_iter = 100
            for i in range(max_iter):
                rmin, rmax = rmax, 2*rmax
                solution = solve_ivp(integrand,[rmin,rmax],[0],rtol=1e-6,atol=1e-6)
                extra = solution.y[0][-1]

                Vsq_out[r_arr > 0] += - prefactor * L * r_eval**L * extra

                if np.allclose(Vsq_out[r_arr > 0],Vsq_out[r_arr > 0]-extra):
                    break

    # Unorder
    unorder = np.argsort(order)

    return Vsq_out[unorder]