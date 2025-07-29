import numpy as np
from inspect import signature
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

from sidmhalo.definitions import no_baryons, integrate, GN

import time as t

# Compute baryon enclosed mass profile from baryon potential
def compute_Mb(Phi_b,rmin,rmax,num=100):
    
    # Number of arguments in Phi_b
    num_variables = len(signature(Phi_b).parameters)

    # Create list of r values
    r_list = np.geomspace(rmin,rmax,num=num)

    # Calculate list of spherically averaged Phi_b
    if num_variables == 1:
        Phi_b_sph_avg_list = np.array([Phi_b(r) for r in r_list])

    elif num_variables == 2:
        def Phi_b_sph_avg_func(r):
            integrand = lambda th: Phi_b(r,th)*np.sin(th)
            return 0.5*integrate(integrand, 0, np.pi)
        Phi_b_sph_avg_list = np.array([Phi_b_sph_avg_func(r) for r in r_list])

    else:
        raise Exception('Case with num_sph_coords=%d not supported.' % num_sph_coords)
        
    # Check if Phi_b is nonzero
    if np.all(Phi_b_sph_avg_list < 0):
        
        # Make interpolation function for y = log(|Phi_b_sph_avg|) as function of x = log(r)
        x = np.log(r_list)
        y = np.log(-Phi_b_sph_avg_list)
        y_interp = InterpolatedUnivariateSpline(x,y)
        dy_dx_interp = y_interp.derivative()
        
        # Phi_b and dPhi_b/dr
        Phi_b_interp = lambda r: - np.exp(y_interp(np.log(r)))

        # dPhi_b/dr
        dPhi_b_dr_interp = lambda r: Phi_b_interp(r) / r * dy_dx_interp(np.log(r))
        
    # Method if previous doesn't work e.g. Phi_b = 0, return 0
    else:
        
        # Make interpolation of Phi_b directly (use linear interp)
        dPhi_b_dr_interp = lambda r: 0

    # Calculate M_b
    # Define M_b recursively to handle case where r is a number or a list/array
    def M_b(r):
        
        # r is a list or array
        if np.ndim(r) > 0:
            return np.array([M_b(ri) for ri in r])
        
        # r is a number
        elif r < rmin:
            Mmin = rmin**2*dPhi_b_dr_interp(rmin)/GN
            return float(Mmin * (r/rmin)**3)

        elif r > r_list[-1]:
            rmax = r_list[-1]
            Mmax = rmax**2*dPhi_b_dr_interp(rmax)/GN
            return float(Mmax)

        else:
            return float(r**2*dPhi_b_dr_interp(r)/GN)

    # Return M_b function
    return M_b

# Calculate spherically averaged potential

def compute_Phi_b_spherical(Mb,rmin,rmax):

    Mmin = Mb(rmin)
    Mmax = Mb(rmax)

    # Let y = Phi(r) - Phi(0)
    RHS = lambda r, y: GN*Mb(r)/r**2

    # Boundary condition y(rmin) = 0.5*GN*M(rmin)/rmin
    # Assuming constant density sphere for r < rmin
    ymin = 0.5*GN*Mmin/rmin
    solution = solve_ivp(RHS, [rmin,rmax], [ymin], dense_output=True, atol=1e-8, rtol=1e-8)
    y_int = solution.sol
    ymax = y_int(rmax)[0]

    # Assume Phi(r) like a point charge for r >= rmax
    # Then Phi(rmax) = y(rmax) + Phi(0) = -GN*Mmax/rmax
    # which implies Phi(0) = -GN*Mmax/rmax - ymax
    Phi_0 = -GN*Mb(rmax)/rmax - ymax


    def Phi_b_out(r):
        
        # Case where r is an array or list
        if np.ndim(r) > 1:
            return np.array([Phi_b(ri) for ri in r])
        
        # Case where r is a single number
        if (r <= rmin) and (r >= 0):
            return float(ymin*(r/rmin)**2 + Phi_0)
            
        elif (r > rmin) and (r < rmax):
            return float(y_int(r)[0] + Phi_0)
            
        elif r >= rmax:
            return float(-GN*Mb(rmax)/r)
        else:
            print(r, 'value for r not valid')
            return 0
        
    return Phi_b_out

def compute_q_baryon(sph_halo,grid=None,R0=0,z0=0,**extraneous):
    
    # Baryon potential
    Phi_b = sph_halo.outer.Phi_b
    
    if len(signature(Phi_b).parameters) == 1:
        return lambda r: np.zeros_like(r)
    
    else:
    
        # Calculate nonsphericity of baryon potential
        if grid is None:
            r200 = sph_halo.outer.r200
            grid = np.geomspace(1e-6*r200,r200)

        R_list = np.array(grid)
        r_list = np.sqrt(R_list**2 + z0**2)
        th1_list = np.arctan2(R_list,z0)
        th2_list = np.arctan2(R_list,-z0)
        Phi_radial_slice = np.abs(np.array([0.5*(Phi_b(r,th1) + Phi_b(r,th2)) for r,th1,th2 in zip(r_list,th1_list,th2_list)]))

        z_list = np.array(grid)
        r_list = np.sqrt(R0**2 + z_list**2)
        th1_list = np.arctan2(R0,z_list)
        th2_list = np.arctan2(R0,-z_list)
        Phi_axial_slice = np.abs(np.array([0.5*(Phi_b(r,th1) + Phi_b(r,th2)) for r,th1,th2 in zip(r_list,th1_list,th2_list)]))


        log_Phi_radial_interp = InterpolatedUnivariateSpline(np.log(R_list),np.log(Phi_radial_slice),ext=3)
        log_Phi_axial_interp = InterpolatedUnivariateSpline(np.log(z_list),np.log(Phi_axial_slice),ext=3)

        rmin = 2 * np.amax([R0,z0,grid[0]])
        rmax = grid[-1]

        r_list = np.geomspace(rmin,rmax,num=100)
        q_list = []

        for r in r_list:

            def f(q):
                logz = 0.5 * np.log(r**2 * q**(4/3) - R0**2 * q**2)
                logR = 0.5 * np.log(r**2 * q**(-2/3) - z0**2 * q**-2)
                return log_Phi_radial_interp(logR) - log_Phi_axial_interp(logz)

            qmin, qmax = 0.5, 2

            while True:

                try:
                    q_list.append(brentq(f,qmin,qmax))
                    break
                except:
                    qmin = 0.5*qmin
                    qmax = 2*qmax

                if (qmax > 10**10):
                    q_list.append(1)
                    break

        return InterpolatedUnivariateSpline(r_list,q_list,ext=3)

def compute_q_iso(sph_halo,q_model=None,**kwargs):
    
    # sph_halo is profile object from spherical Jeans model
    r1 = sph_halo.r1
    r200 = sph_halo.outer.r200
    
    # Define array of r points for calculating q_iso
    r = np.geomspace(1e-6*r200,r200)
    
    # Calculate DM fraction vs radius
    # Equivalent to fdm = Mdm(r) / (Mdm(r) + Mb(r))
    fdm = 1 / (1 + sph_halo.outer.M_b(r)/sph_halo.M_encl(r))
    
    # Nonsphericity q for baryon potential
    q_baryon = compute_q_baryon(sph_halo,**kwargs)
    qb = q_baryon(r)
    # print(qb)
    
    if not(q_model):
        
        def q_model(qb,f):
            F = f + f * (1 - f) * (-1 + 0.5 * qb)
            return qb * (1 - F) + F
    else:
        pass
    
    q_iso = InterpolatedUnivariateSpline(r,q_model(qb,fdm),ext=3)
    
    return q_iso

def compute_q_eff(sph_halo,q0,Nm=1,**kwargs):
    
    # sph_halo is profile object from spherical Jeans model
    r1 = sph_halo.r1
    r200 = sph_halo.outer.r200
    
    # Define array of r points for calculating q_iso
    r = np.geomspace(1e-6*r200,r200)

    # Define collisionless q0
    if callable(q0):
        q0_func = lambda r: q0(r)
    else:
        q0_func = lambda r: q0
    
    # Calculate number of scatters vs r
    if r1 > 0:
        
        N = lambda r: Nm * sph_halo.rho_sph_avg(r) / sph_halo.rho_sph_avg(r1)
        
        # Isothermal q profile
        q_iso = compute_q_iso(sph_halo,**kwargs)
        
        # Result from toy model calculation
        # k0 = 1 / (np.sqrt(2) * 5)
        k0 = np.sqrt(2) / 5
        
        def q_eff(r):
            log_qiso = np.log(q_iso(r))
            log_q0 = np.log(q0_func(r))
            log_q = log_qiso + (log_q0 - log_qiso) * np.exp( - k0 * N(r) )
            return np.exp(log_q)
        
    else:
        q_eff = np.vectorize(q0_func)
        
    return q_eff

def compute_r_sph_grid(outer_halo,q_func,numr=30,numth=20):
    
    r200 = outer_halo.r200
    r_list = np.geomspace(1e-6*r200,r200,num=numr)
    th_list = np.linspace(0,np.pi/2,num=numth)
    
    sph_grid = np.zeros((numr,numth))
    
    # Spheroidal radius function    
    def r_sph(r,th,max_iter=200,k=1,**kwargs):

        steps_to_damp = max_iter / 10
    
        R = r * np.sin(th)
        z = r * np.cos(th)
    
        rsph_old = r
        # start = t.time()
        for i in range(max_iter):
            
            rsph_new_calculated = np.sqrt(R**2 * q_func(rsph_old)**(2/3) + z**2 * q_func(rsph_old)**(-4/3))

            rsph_new = (1 - k) * rsph_old + k * rsph_new_calculated

            if np.allclose(rsph_new,rsph_old,**kwargs):
                break

            rsph_old = rsph_new

            # if not converged within tolerance reduce the damping factor k by 0.1 every max_iter/10 iterations
            if i % steps_to_damp == 0: 
                k -= 0.1 
                
        else:
            raise Exception('rsph did not converge within tolerance with max_iter=%d iterations for (r,th)=(%f,%f)' % (max_iter,r,th))
        
        # end = t.time()
        # print('Time taken for r=%f, th=%f is %f' % (r,th,end-start))
        return rsph_new
    
    # Calculate 
    for i in range(numr):
        for j in range(numth):
            sph_grid[i,j] = r_sph(r_list[i],th_list[j])/r_list[i]
    
    # Make interpolating function
    logr_list = np.log(r_list)
    interp = RectBivariateSpline(logr_list,th_list,sph_grid)
    
    def r_sph_interp(r,th):
        # Symmetrize in theta about np.pi/2
        th_abs = np.arccos(np.abs(np.cos(th)))
        return interp(np.log(r),th_abs,grid=False) * r
    
    return r_sph_interp

