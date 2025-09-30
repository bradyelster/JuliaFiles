# -*- coding: utf-8 -*-
"""
Created on Wed Jul 16 12:00:21 2025

@author: gaboa
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.interpolate import PPoly
from netCDF4 import Dataset
import time
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.markers import JoinStyle
from matplotlib.markers import MarkerStyle
from copy import deepcopy
from diffeqpy import de 
from scipy.integrate import quad
from scipy.integrate import simpson


def parityfunc(x,flag):
    if flag==1:
        return x 
    else:
        return 1
    
def siestaAngleInterpolation(quantity, half_grid):
    # Want to interpolate for negative s which is the equivalent of
    # shifting the poloidal angle by pi. Also want to be able to interpolate
    # across the axis since we do not have this information
    half_grid_neg = - half_grid[::-1]
    x = np.append(half_grid_neg, half_grid)
    interpolations = np.empty((n_mode,m_mode), dtype='object_')
    for n in range(0,n_mode, 1):
        for m in range(0,m_mode,1):
            temp_quantity = quantity[1:,n,m]
            
            if m % 2 == 0 :
                fac = 1 
            else:
                fac = -1 
            y = np.append(fac*temp_quantity[::-1],temp_quantity)
            
            interpolations[n,m] = CubicSpline(x, y, bc_type='natural')
            
    return interpolations

def siestaRadialInterpolation(quantity,u_inter,v_inter):
    interpolations = np.empty((n_mode,m_mode),dtype='object_')
    grid = u_inter[0,0].x
    coeffs = np.zeros((5,len(grid)-1))
    for n in range(0,n_mode, 1):
        for m in range(0,m_mode,1):
            # s^3 terms for radial, and s^2 terms for angular
            coeffs[0,:] = (m*u_inter[n,m].c[0,:]+nfp*tor_modes[n]*v_inter[n,m].c[0,:])*1/4
            # s^3 terms for radial, and s^2 terms for angular
            coeffs[1,:] = (m*u_inter[n,m].c[1,:]+nfp*tor_modes[n]*v_inter[n,m].c[1,:])*1/3
            # s^2 terms for radial, and s^1 terms for angular
            coeffs[2,:] = (m*u_inter[n,m].c[2,:]+nfp*tor_modes[n]*v_inter[n,m].c[2,:])*1/2
            # s^1 terms for radial, and s^0 terms for angular
            coeffs[3,:] =  m*u_inter[n,m].c[3,:]+nfp*tor_modes[n]*v_inter[n,m].c[3,:]
            # s terms for radial, which is just A(s) on the grid
            temp_quantity = quantity[1:,n,m]
            if m % 2 == 0 :
                fac = 1 
            else:
                fac = -1 
            temp_quantity = np.append(fac*temp_quantity[::-1],temp_quantity)
            # The coefficients that are not multiplying s are the left points of 
            # each spline, hence the last point on the right will not be used.
            coeffs[4,:] = temp_quantity[:-1]
            
            interpolations[n,m] = deepcopy(PPoly(coeffs,grid))
            
    return interpolations

def siestaRZInterpolation(quantity, full_grid):
    full_grid = full_grid[1::]
    full_grid_neg = - full_grid[::-1]
    x = np.append(full_grid_neg, full_grid)
    interpolations = np.empty((n_mode,m_mode), dtype='object')
    for n in range(0,n_mode, 1):
        for m in range(0,m_mode,1):
            temp_quantity = quantity[1:,n,m]
            
            if m % 2 == 0 :
                fac = 1 
            else:
                fac = -1 
            y = np.append(fac*temp_quantity[::-1],temp_quantity)
            
            interpolations[n,m] = CubicSpline(x, y)
    
    return interpolations

def siestaTransform(quantity, ss, uu, vv, dflag, parity, ihalf ,normflag):
    """
    Transforms the Fourier quantity to real space.
    Parameters
    ----------
    quantity : Array of interpolation objects with shape (n_mode,m_mode)
        The array of all the interpolated Fourrier coefficients.
    ss : float between (-1,1)
        Radial coordinate to evaluate the quantity at.
    uu : float
        Poloidal angle to evaluate the quantity at.
    vv : float
        Toroidal angle to evaluate the quantity at.
    dflag : INT
        Designates kind of derivative being taken. 
        0 takes no derivative
        1 takes a radial derivative
        2 takes a poloidal derivative
        3 takes a toroidal derivative
    parity : INT
        Designates sine or cosine parity. 
        0 indicates cos parity 
        1 indicates sine parity
    ihalf : INT
        Designates if the quantity is on the half grid or on the full grid.
        0 indicates it is on the full grid
        1 indicates it is on the half grid
    normflag : INT
        Designates if quantity is normalized and includes the jacobian. 
        0 indicates it is not normalized
        1 indicates it is normalized

    Returns
    -------
    None.

    """
    
    # First we take care of the normalization facor
    if normflag == True:
        ntor_0_mode_ind = np.where(tor_modes == 0)[0][0]
        ofac = np.full((n_mode,m_mode),np.sqrt(2))
        ofac[ntor_0_mode_ind,0] = 1
        bfac = bfactor
    else:
        ofac = np.ones((n_mode,m_mode))
        bfac = 1
    
    if dflag == 1:
        sdfac = 1
        u_parity_flag = 0
        v_parity_flag = 0  
    
    elif dflag == 2:
        parity += 1
        sdfac = 0
        u_parity_flag = 1
        v_parity_flag = 0    
        
    elif dflag == 3:
        parity += 1
        sdfac = 0
        u_parity_flag = 0
        v_parity_flag = 1
    else:
        sdfac = 0
        u_parity_flag = 0
        v_parity_flag = 0           
    
    temp = 0
    if parity % 2 == 1:
        for n in range(0,n_mode):
            for m in range(0, m_mode):
                temp += ofac[n,m]*quantity[n,m](ss,sdfac)      \
                        *np.sin(m*uu+nfp*tor_modes[n]*vv)      \
                        *parityfunc(-m,u_parity_flag)          \
                        *parityfunc(-nfp*tor_modes[n],v_parity_flag)
        value = temp / bfac
    else:
        for n in range(0,n_mode):
            for m in range(0, m_mode):
                temp += ofac[n,m]*quantity[n,m](ss,sdfac)      \
                        *np.cos(m*uu+nfp*tor_modes[n]*vv)      \
                        *parityfunc(m,u_parity_flag)          \
                        *parityfunc(nfp*tor_modes[n],v_parity_flag)
        value = temp / bfac
    
    return value

def siestaBint(ss,uu,vv):
    f_cos = 0; f_sin = 1
    f_none = 0; 
    f_half = 1; f_norm = 1

    bs = siestaTransform(s_interpolation, ss, uu, vv, f_none, f_sin, f_half, f_norm) # SUM MN-series at target point to get contravariant Bs
    bu = siestaTransform(u_interpolation, ss, uu, vv, f_none, f_cos, f_half, f_norm) # SUM MN-series at target point to get contravariant Bu
    bv = siestaTransform(v_interpolation, ss, uu, vv, f_none, f_cos, f_half, f_norm) # SUM MN-series at target point to get contravariant Bv
    return [bs,bu,bv]

def divBint(ss,uu,vv):
    f_cos = 0; f_sin = 1
    f_none = 0; f_ds = 1; f_du = 2; f_dv = 3
    f_half = 1
    f_norm = 1
    
    dJBsds = siestaTransform(s_interpolation, ss, uu, vv, f_ds, f_sin, f_none ,f_norm)
    dJBudu = siestaTransform(u_interpolation, ss, uu, vv, f_du, f_cos, f_none ,f_norm)
    dJBvdv = siestaTransform(v_interpolation, ss, uu, vv, f_dv, f_cos, f_none ,f_norm)
    
    return dJBsds+dJBudu+dJBvdv

def siestaBprimeint(ss,uu,vv,dflag):
    f_cos = 0; f_sin = 1
    f_none = 0; f_ds = 1; f_du = 2; f_dv = 4
    f_half = 1
    f_norm = 1

      
    bs_prime = siestaTransform(s_interpolation, ss, uu, vv, dflag, f_sin, f_half, f_norm) # SUM MN-series at target point to get contravariant Bs
    bu_prime = siestaTransform(u_interpolation, ss, uu, vv, dflag, f_cos, f_half, f_norm) # SUM MN-series at target point to get contravariant Bu
    bv_prime = siestaTransform(v_interpolation, ss, uu, vv, dflag, f_cos, f_half, f_norm) # SUM MN-series at target point to get contravariant Bv
    
    return [bs_prime,bu_prime,bv_prime]

def integration_fun_int(phi, stheta):
    parr = np.array(phi)
    mag_stp = siestaBint(stheta[0],stheta[1],parr)
    
    # return the rhs of the differential equation
    return np.array((mag_stp[0]/mag_stp[2], mag_stp[1]/mag_stp[2]))

def integration_fun_newton(phi, y):
    f_ds = 1 
    f_du = 2
    phi_arr = np.array(phi)
    
    M_temp = np.zeros((2,2))
    ss = y[0]
    uu = y[1]
    
    M_temp[0][0] = y[2]
    M_temp[0][1] = y[3]
    M_temp[1][0] = y[4]
    M_temp[1][1] = y[5]
    
    # Create C matrix
    C = np.zeros((2,2))
    
    # Extract magnetic field and its derivatives at mapped point
    B    = siestaBint(ss,uu,phi_arr)
    dBds = siestaBprimeint(ss,uu,phi_arr, f_ds)
    dBdu = siestaBprimeint(ss,uu,phi_arr, f_du)
    
    # Calculate the C matrix
    C[0][0] = (B[2]*dBds[0]-B[0]*dBds[2] )/np.square(B[2])
    C[0][1] = (B[2]*dBdu[0]-B[0]*dBdu[2] )/np.square(B[2])
    C[1][0] = (B[2]*dBds[1]-B[1]*dBds[2] )/np.square(B[2])
    C[1][1] = (B[2]*dBdu[1]-B[1]*dBdu[2] )/np.square(B[2])
    
    cdotm = np.matmul(C,M_temp)
    
    # return the rhs of the differential equation
    return np.array([B[0]/B[2], 
                     B[1]/B[2], 
                     cdotm[0][0],
                     cdotm[0][1],
                     cdotm[1][0],
                     cdotm[1][1]
                     ])

def convert_int(ss_in,uu_in,vv_in):
    f_cos = 0 
    f_sin = 1
    f_none = 0
    f_half = 1

    rr_in = siestaTransform(r_interpolation, ss_in, uu_in, vv_in, f_none, f_cos, f_half,f_none)  
    zz_in = siestaTransform(z_interpolation, ss_in, uu_in, vv_in, f_none, f_sin, f_half,f_none)
    
    return [rr_in, zz_in]

def field_line_follow(ss_in,uu_in,vv_in=0.0, niter=100, nfp=1, nstep=100, save_all=False, convert=False, method='LSODA', rtol=1e-6):
    """
    Description:
        Follows the magnetic field lines of points in the equilibrium to generate
        data for a poincare plot. At the same time approximately calculates the
        rotational transform, iota, using a least squares fit.
    ----------
    Parameters
    ----------
    ss_in : List
        List of intitial s values that will be itterated over.
    uu_in : List
        List of coresponding u values for each s coordinate, typically array of 0's.
    vv_in : float, optional
        The intial toroidal angle used in the field line following. The default is 0.0.
    niter : int, optional
        Number of toroidal revolutions across a field period. The default is 100.
    nfp : int, optional
        Number of field periods. The default is 1.
    nstep : int, optional
        Number of steps taken in each field period. The default is 100.
    save_all : bool, optional
        Option to save all points in tracing. Otherwise only save the points
        at each field period. The default is False
    convert : bool, optional
        Option to convert the coordinates from (s,u,v) to (R,Z,phi). Defaults
        to False.
    method : string, optional
        Method solve_ivp will use. See documentation of solve_ivp for specific optrions.
        Defaults to LSODA
    rtol   : float, optional
        Relative tolerance to use in solve_ivp. Defaults to 1e-06.
        
    Returns
    -------
    lines, a list of points that to be used to generate a poincare plot.
    
    iota, a list of iota values for each surface with the associated s value

    """
    # Set up run parameters
    dphi = 2 * np.pi / nfp / nstep # Step size
    phi  = vv_in + dphi * nstep * np.arange(niter) # phi angles at each field period
    nlines = len(ss_in) # Number of surfaces
    lines  = [0]*nlines # List to contain points in each surface
    iota = [0]*nlines
    for i in range(nlines):  # loop over each field-line
        least_squares = [0, 0, 0, uu_in[0],1]
        points = [0]*(niter+1)
        points[0] = [ss_in[i], uu_in[i],0]
        for j in range(niter):  # loop over each toroidal iteration
           print_progress(i * niter + j + 1, nlines * niter)
           su = [points[j][0],points[j][1]]
           phi_start = phi[j]
           for k in range(nstep):  # loop inside one iteration
               sol = solve_ivp(integration_fun_int, (phi_start, phi_start + dphi), su, rtol=rtol, method=method)
               su = sol.y[:, -1]
               phi_start += dphi
           points[j+1] = [su[0],su[1],2*np.pi/nfp*(j+1)]
           # Here we begin the summations that will go into the least squares 
           # calculation of the rotational transform. The equation we are
           # finding the slope of is theta(phi) = theta_0 + iota* phi. Although
           # in our coordinate system our theta is u and phi is v. Using
           # the formula for slope from a least squares regression we have that
           # x is phi, and y is theta
           least_squares[0] += phi_start**2       # Sum of x^2
           least_squares[1] += phi_start          # Sum of x
           least_squares[2] += phi_start*su[-1]    # Sum of x*y
           least_squares[3] += su[-1]              # Sum of y
           least_squares[4] += 1                  # Number of points
        # Below we calculate the value of iota, the formula being:
        # (N * Sum(x*y) - Sum(x)*Sum(y)) / ( N * Sum(x^2) - (Sum(x))^2 )
        iota_temp = (least_squares[4] * least_squares[2] - least_squares[1] * least_squares[3] ) / (least_squares[4] * least_squares[0] - least_squares[1]*least_squares[1])
        iota[i] =  [[ss_in[i],uu_in[i],0], iota_temp]
        lines[i] = np.array(points)
        
        # TODO implement conversion and save all functionalities
    
    return lines, iota

def fixed_point_solver(su,n_mode,m_mode,nfp,tol=1e-8, max_iter=100, nstep=50, damp=0.1):
    """

    Parameters
    ----------
    su : TYPE
        DESCRIPTION.
    n_mode : TYPE
        DESCRIPTION.
    m_mode : TYPE
        DESCRIPTION.
    nfp : TYPE
        DESCRIPTION.
    tol : TYPE, optional
        DESCRIPTION. The default is 1e-8.
    max_iter : TYPE, optional
        DESCRIPTION. The default is 100.
    nstep:
        
    damp: 

    Returns
    -------
    None.

    """
    ss_in = su[0]
    uu_in = su[1] 
    vv_in = 0.0
    # phi_start = vv_in
    dphi = 2*np.pi/ nfp / nstep
    #M = np.eye(2)
    identity_matrix= np.eye(2)
    
    
    for i in range(max_iter):
        phi_start = vv_in
        M = np.eye(2)
        y0 = [ss_in,uu_in,M[0][0], M[0][1],  M[1][0], M[1][1]]
        # Map the point m field periods away 
        for j in range(m_mode): # Loop over all the field periods
            for k in range(nstep):  # Loop inside each field period
                sol = solve_ivp(integration_fun_newton, (phi_start, phi_start + dphi), y0, rtol=1e-8, atol=1e-8, method='LSODA')
                y0 = sol.y[:, -1]
                phi_start += dphi
        
        end_time = time.time()
        # Compute elapsed time
        elapsed_time = end_time - start_time  
        hours = int(elapsed_time // 3600)
        minutes = int((elapsed_time % 3600) // 60)
        seconds = elapsed_time % 60

        # Print result
        print(f"Time elapsed: {hours} hours, {minutes} minutes, {seconds:.2f} seconds")
        
        # Extract mapped point
        ss_map = y0[0]
        uu_map = y0[1] % (2*np.pi)
        
        if uu_map > (np.pi):
            uu_map -= 2*np.pi
        
        print("Input")
        print(ss_in,uu_in)
        
        print("Map")    
        print(ss_map, uu_map)
        
        #v_map = phi_start
        
        # Below is the linearized mapping "Jacobian" matrix. AKA tangent map.
        M[0][0] = y0[2]
        M[0][1] = y0[3]
        M[1][0] = y0[4]
        M[1][1] = y0[5]
        
        # We can calculate the residue with this 
        residue = 1/2 - 1/4 * np.trace(M)
        
        # Get numerator of Newton term
        delta_s = ss_map - ss_in
        delta_u = uu_map - uu_in 
        
        # g(x) = f(x) - x
        g = np.array([delta_s, delta_u])
        
        
        # Get derivative of g(x) = f(x) - x; dgdx = dfdx - 1
        g_prime_temp = np.subtract(M,identity_matrix)
        
        
        # Get determinant
        g_det = np.linalg.det(g_prime_temp)
        
        # Get inverse of derivative of g(x) 
        if abs(g_det) <= 1*10^(-12):
            print("Small determinant")
            g_det = 1*10^(-12)
            g_prime_inv = np.zeros((2,2))
            g_prime_inv[0][0] = M[1][1] / g_det
            g_prime_inv[0][1] = - M[0][1] / g_det
            g_prime_inv[1][0] = - M[1][0] / g_det
            g_prime_inv[1][1] = M[0][0] / g_det
            
            # g(x)/g'(x)
            g_frac = -np.matmul(g_prime_inv,g)
        
        else:
            # Solve x J = g, for x
            g_frac = np.linalg.solve(g_prime_temp,-g)
        
        su_start = np.array([ss_in,uu_in])
        
        # New s and u coordinates
        su_new = np.add(su_start,damp*g_frac)
        #su_new[1] = su_new[1] # % (2*np.pi)
        delta_s = abs(su_new[0] - su_start[0])
        delta_u = abs(su_new[1] - su_start[1])
        print('New Point Change')
        print(delta_s, delta_u )
        print()
        
        
        if delta_s < (tol*damp) and delta_u < (tol*damp*10) :
            print("Desired tolerance was reached.")
            M_det = np.linalg.det(M)
            return su_new, residue, M_det
        
        ss_in = su_new[0]
        uu_in = su_new[1]
        
        #print('New Point')
        #print( ss_in, uu_in)
    
    print("Reached maximum number of iterations. Run with larger max_iter")
    M_det = np.linalg.det(M)
    return su_new, residue, M_det

def int_fun_newton_julia(du,u,p,t):
    f_ds = 1 
    f_du = 2
    phi_arr = t
    M_temp = np.zeros((2,2))
    
    ss = u[0]
    uu = u[1]
    
    M_temp[0][0] = u[2]
    M_temp[0][1] = u[3]
    M_temp[1][0] = u[4]
    M_temp[1][1] = u[5]
    

    #siestaBint, siestaBprimeint = p
    
    # Create C matrix
    C = np.zeros((2,2))
    
    # Extract magnetic field and its derivatives at mapped point
    B    = siestaBint(ss,uu,phi_arr)
    dBds = siestaBprimeint(ss,uu,phi_arr, f_ds)
    dBdu = siestaBprimeint(ss,uu,phi_arr, f_du)
    
    # Calculate the C matrix
    C[0][0] = (B[2]*dBds[0]-B[0]*dBds[2] )/np.square(B[2])
    C[0][1] = (B[2]*dBdu[0]-B[0]*dBdu[2] )/np.square(B[2])
    C[1][0] = (B[2]*dBds[1]-B[1]*dBds[2] )/np.square(B[2])
    C[1][1] = (B[2]*dBdu[1]-B[1]*dBdu[2] )/np.square(B[2])
    
    cdotm = np.matmul(C,M_temp)
    
    # return the rhs of the differential equation
    du[0]= B[0]/B[2] 
    du[1]= B[1]/B[2] 
    du[2]= cdotm[0][0]
    du[3]= cdotm[0][1]
    du[4]= cdotm[1][0]
    du[5]= cdotm[1][1]
    
def fixed_point_solver_julia(su,n_mode,m_mode,nfp,tol=1e-8, max_iter=100, nstep=50, damp=0.1):
    """

    Parameters
    ----------
    su : TYPE
        DESCRIPTION.
    n_mode : TYPE
        DESCRIPTION.
    m_mode : TYPE
        DESCRIPTION.
    nfp : TYPE
        DESCRIPTION.
    tol : TYPE, optional
        DESCRIPTION. The default is 1e-8.
    max_iter : TYPE, optional
        DESCRIPTION. The default is 100.
    nstep:
        
    damp: 

    Returns
    -------
    None.

    """
    ss_in = su[0]
    uu_in = su[1] 
    vv_in = 0.0
    # phi_start = vv_in
    dphi = 2*np.pi/ nfp / nstep
    #M = np.eye(2)
    identity_matrix= np.eye(2)
    #p = (siestaBint, siestaBprimeint)
    tsidas_alg = de.AutoVern7(de.Rodas5())
    
    for i in range(max_iter):
        phi_start = vv_in
        M = np.eye(2)
        y0 = [ss_in,uu_in,M[0][0], M[0][1],  M[1][0], M[1][1]]
        # Map the point m field periods away 
        for j in range(m_mode): # Loop over all the field periods
            for k in range(nstep):  # Loop inside each field period
                prob = de.ODEProblem(int_fun_newton_julia,y0,(phi_start, phi_start + dphi))
                sol = de.solve(prob, tsidas_alg,abstol=1e-8,reltol=1e-8)
                y0 = sol.u[-1]
                phi_start += dphi
        
        end_time = time.time()
        # Compute elapsed time
        elapsed_time = end_time - start_time  
        hours = int(elapsed_time // 3600)
        minutes = int((elapsed_time % 3600) // 60)
        seconds = elapsed_time % 60

        # Print result
        print(f"Time elapsed: {hours} hours, {minutes} minutes, {seconds:.2f} seconds")
        
        # Extract mapped point
        ss_map = y0[0]
        uu_map = y0[1] % (2*np.pi)
        
        if uu_map > (np.pi):
            uu_map -= 2*np.pi
        
        print("Input")
        print(ss_in,uu_in)
        
        print("Map")    
        print(ss_map, uu_map)
        
        #v_map = phi_start
        
        # Below is the linearized mapping "Jacobian" matrix. AKA tangent map.
        M[0][0] = y0[2]
        M[0][1] = y0[3]
        M[1][0] = y0[4]
        M[1][1] = y0[5]
        
        # We can calculate the residue with this 
        residue = 1/2 - 1/4 * np.trace(M)
        
        # Get numerator of Newton term
        delta_s = ss_map - ss_in
        delta_u = uu_map - uu_in 
        
        # g(x) = f(x) - x
        g = np.array([delta_s, delta_u])
        
        
        # Get derivative of g(x) = f(x) - x; dgdx = dfdx - 1
        g_prime_temp = np.subtract(M,identity_matrix)
        
        
        # Get determinant
        g_det = np.linalg.det(g_prime_temp)
        
        # Get inverse of derivative of g(x) 
        if abs(g_det) <= 1*10^(-12):
            print("Small determinant")
            g_det = 1*10^(-12)
            g_prime_inv = np.zeros((2,2))
            g_prime_inv[0][0] = M[1][1] / g_det
            g_prime_inv[0][1] = - M[0][1] / g_det
            g_prime_inv[1][0] = - M[1][0] / g_det
            g_prime_inv[1][1] = M[0][0] / g_det
            
            # g(x)/g'(x)
            g_frac = -np.matmul(g_prime_inv,g)
        
        else:
            # Solve x J = g, for x
            g_frac = np.linalg.solve(g_prime_temp,-g)
        
        su_start = np.array([ss_in,uu_in])
        
        # New s and u coordinates
        su_new = np.add(su_start,damp*g_frac)
        #su_new[1] = su_new[1] # % (2*np.pi)
        delta_s = abs(su_new[0] - su_start[0])
        delta_u = abs(su_new[1] - su_start[1])
        print('New Point Change')
        print(delta_s, delta_u )
        print()
        
        
        if delta_s < (tol*damp) and delta_u < (tol*damp*10) :
            print("Desired tolerance was reached.")
            M_det = np.linalg.det(M)
            return su_new, residue, M_det
        
        ss_in = su_new[0]
        uu_in = su_new[1]
        
        #print('New Point')
        #print( ss_in, uu_in)
    
    print("Reached maximum number of iterations. Run with larger max_iter")
    M_det = np.linalg.det(M)
    return su_new, residue, M_det

def int_fun_julia(du,u,p,t):
    mag_stp = siestaBint(u[0],u[1],t)
    du[0] = mag_stp[0]/mag_stp[2]
    du[1] = mag_stp[1]/mag_stp[2]

def field_line_follow_julia(ss_in,uu_in,vv_in=0.0, niter=100, nfp=1, nstep=100, save_all=False, convert=False, method='LSODA', rtol=1e-6):
    """
    Description:
        Follows the magnetic field lines of points in the equilibrium to generate
        data for a poincare plot. At the same time approximately calculates the
        rotational transform, iota, using a least squares fit.
    ----------
    Parameters
    ----------
    ss_in : List
        List of intitial s values that will be itterated over.
    uu_in : List
        List of coresponding u values for each s coordinate, typically array of 0's.
    vv_in : float, optional
        The intial toroidal angle used in the field line following. The default is 0.0.
    niter : int, optional
        Number of toroidal revolutions across a field period. The default is 100.
    nfp : int, optional
        Number of field periods. The default is 1.
    nstep : int, optional
        Number of steps taken in each field period. The default is 100.
    save_all : bool, optional
        Option to save all points in tracing. Otherwise only save the points
        at each field period. The default is False
    convert : bool, optional
        Option to convert the coordinates from (s,u,v) to (R,Z,phi). Defaults
        to False.
    method : string, optional
        Method solve_ivp will use. See documentation of solve_ivp for specific optrions.
        Defaults to LSODA
    rtol   : float, optional
        Relative tolerance to use in solve_ivp. Defaults to 1e-06.
        
    Returns
    -------
    lines, a list of points that to be used to generate a poincare plot.
    
    iota, a list of iota values for each surface with the associated s value

    """
    # Set up run parameters
    pi2 = 2*np.pi
    phi_end = pi2*niter
    dphi = pi2 / nfp / nstep # Step size
    phis = np.arange(0,phi_end, dphi)
    phi_save = np.arange(0,phi_end, pi2)
    #phi  = vv_in + dphi * nstep * np.arange(niter) # phi angles at each field period
    nlines = len(ss_in) # Number of surfaces
    lines  = [0]*nlines # List to contain points in each surface
    iota = [0]*nlines
    
    def prob_func(prob, i, rep):
        return de.remake(prob, u0=[ss_in[i-1],uu_in[i-1]])

    tsidas_alg = de.AutoVern7(de.Rodas5())
    u0 = [ss_in[0],uu_in[0]]
    phi_span = (0,phi_end)
    prob = de.ODEProblem(int_fun_julia, u0, phi_span)
    ensembleprob = de.EnsembleProblem(prob, prob_func=prob_func, safetycopy=False)
    sol = de.solve(ensembleprob,tsidas_alg,de.EnsembleThreads(),trajectories=nlines, saveat=phi_save, tstops=phis)
    
    for i in range(nlines):
        line = [0]*niter
        for j in range(niter):
            line[j] = [sol.u[i].u[j][0],sol.u[i].u[j][1], sol.u[i].t[j]]
        lines[i] = line 
        temp_array = np.asarray(line)
      
        iota[i] = [temp_array[0,0],simpson(temp_array[:,1],x=temp_array[:,2])/simpson(temp_array[:,2],x=temp_array[:,2])]
      
    return np.asarray(lines), np.asarray(iota)

if __name__ == '__main__':
    '''
    Read in the data from the SIESTA restart file
    '''
    # Insert file path
    #filein = "C:\\Users\\gaboa\\Documents\\School_Stuff\\Research\\Ware\\SIESTA\\siesta_restart_rotatingellipse_island.nc"
    #filein = "C:\\Users\\gaboa\\Documents\\School_Stuff\\Research\\Ware\\SIESTA\\siesta_restart_cth_140_M5_N7_island_pert_1Em3.nc"
    filein = "C:\\Users\\gaboa\\Documents\\School_Stuff\\Research\\Ware\\SIESTA\\siesta_restart_cth_120_island_M5_N7_1Em3_resistive.nc"
    
        
    siesta = Dataset(filein,"r")
    m_mode = siesta.dimensions["m-mode"].size 
    n_mode = siesta.dimensions["n-mode"].size 
    radius = siesta.dimensions["radius"].size 
    # print('m_mode,n_mode,radius=',m_mode,n_mode,radius)

    nfp = 1*siesta.variables["nfp"][:]
    ns = 1*siesta.variables["nrad"][:]
    mpol = 1*siesta.variables["mpol"][:]
    ntor = 1*siesta.variables["ntor"][:]
    tor_modes = 1*siesta.variables["tor_modes"][:].data

    rmnc = siesta.variables["rmnc_m_n_r_"][:].data            # (radius, n-mode, m-mode)
    zmns = siesta.variables["zmns_m_n_r_"][:].data            # (radius, n-mode, m-mode)

    bsupsmnsh = siesta.variables["bsupsmnsh_m_n_r_"][:].data    # (radius, n-mode, m-mode)
    bsupumnch = siesta.variables["bsupumnch_m_n_r_"][:].data    # (radius, n-mode, m-mode)
    bsupvmnch = siesta.variables["bsupvmnch_m_n_r_"][:].data    # (radius, n-mode, m-mode)
    
    # Below are the values of the magnetic field with the jacobian. These
    # come straight from the internals of siesta. The other numbers are
    # processed to remove the jacobian
    jbsupsmnsh = siesta.variables["JBsupssh_m_n_r_"][:].data    # (radius, n-mode, m-mode)
    jbsupumnch = siesta.variables["JBsupuch_m_n_r_"][:].data    # (radius, n-mode, m-mode) 
    jbsupvmnch = siesta.variables["JBsupvch_m_n_r_"][:].data    # (radius, n-mode, m-mode)
    
    # Normalization factor in jb terms
    bfactor = siesta.variables["b_factor"][:].data
    
    siesta_full = np.linspace(0,1, ns)
    siesta_half = (siesta_full[:-1] + siesta_full[1:])/2
    
    u_interpolation = siestaAngleInterpolation(jbsupumnch, siesta_half)
    v_interpolation = siestaAngleInterpolation(jbsupvmnch, siesta_half)
    s_interpolation = siestaRadialInterpolation(jbsupsmnsh,u_interpolation,v_interpolation)
    
    r_interpolation = siestaRZInterpolation(rmnc, siesta_full)
    z_interpolation = siestaRZInterpolation(zmns, siesta_full)
    
    #%%
    '''
    Finding the Fixed Point
    '''
    # Number of field periods
    #nfp = 5
    # Fixed Point M mode number
    m_mode_f = 2
    n_mode_f = 1
    
    # First set initial guess for the axis, assume toroidal angle of 0
    # Rottating Ellipse
    # ss_g, uu_g = 0.6678, 0
    # ss_g, uu_g = 0.67, 0
    
    # CTH Case 120
    ss_g  = 9.0723236703960200E-01  #0.929 
    uu_g  = 1.9946645720803600E+00  #2.011583498026671801e+00 
    
    print(siestaBint(ss_g,uu_g,0))
    print(divBint(ss_g,uu_g,0))
    # CTH Case 140
    # ss_g =  0.8551050130771877 #0.86
    # uu_g =  1.8668276838258855 #1.87
     
    start_time = time.time()
    # 2D
    fixed_point, residue, M_det = fixed_point_solver_julia([ss_g,uu_g],n_mode=1 ,m_mode=m_mode_f, nfp=1 , tol=1e-8, max_iter=200, damp=1)
    # 1D
    #fixed_point = fixed_point_solver_1D(ss_g,n_mode_f,m_mode_f,1)
    
    end_time = time.time()
    # Compute elapsed time
    elapsed_time = end_time - start_time  
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60

    # Print result
    print(f"Time elapsed: {hours} hours, {minutes} minutes, {seconds:.2f} seconds")
    #%%
    """
    Fixed Point in real space and accuracy
    """
    # 2D
    real_fixed_point = convert_int(fixed_point[0],fixed_point[1],0)
    
    # 1D
    #real_fixed_point = convert(fixed_point,0,0)
    
    # We can check if a revolution gives the same or close enough to the same
    nsteps = 100
    phi_start = 0
    dphi = 2*np.pi/nfp/nsteps
    
    # 2D
    st = fixed_point
    # 1D
    #st = [fixed_point, 0.0]
    
    for j in range(m_mode_f):
        for k in range(nsteps):  # loop inside one iteration
            sol = solve_ivp(integration_fun_int, (phi_start, phi_start + dphi), st, method='LSODA')
            st = sol.y[:, -1]
            phi_start += dphi
        #rev = solve_ivp(integration_fun, t_span=[0, 2*np.pi/5], y0=[axis[0], axis[1]], t_eval=[0, 2*np.pi/5],  rtol=rtol, atol=atol, method=method)
    fixed_point_m_rev = [st[0],st[1],phi_start]
    
    fixed_point_m_rev_real = convert_int(fixed_point_m_rev[0],fixed_point_m_rev[1],fixed_point_m_rev[2])
    diff_r = real_fixed_point[0]-fixed_point_m_rev_real[0]
    diff_z = real_fixed_point[1]-fixed_point_m_rev_real[1]
    distance = np.sqrt(diff_r**2+diff_z**2)
    
    # Print out the distance between the points after one revolution
    print("After m revolutions, the points are, ", distance, " meters away")
    print(residue)
    print(fixed_point[0], fixed_point[1])
    print(M_det)
    
    #%%
    """
    Field line following
    """
    #nsurf = 5
    #s = np.linspace(0.88, 0.92, nsurf)
    s = siesta_half[70:]
    theta = [0]*len(s)
    phi0 = 0.0
    
    # Define number of steps to save and nuber of field periods
    nstep = 1 # Number of steps in one field period
    niter = 750   # Number of toroidal periods in tracing
    nfp = 1 # nfp is extracted earlier but in SIESTA cases with symmetry breaking islands it will be set to 1
    method = 'LSODA' #'RK45' 
    
    # some settings
    print("Begin field-line tracing: ")
    start_time = time.time()
    #surfaces_s, iota_s = field_line_follow(s,theta, phi0, niter, nfp, nstep, rtol=1e-6)
    surfaces, iota = field_line_follow_julia(s,theta, phi0, niter, nfp, nstep, rtol=1e-6)
    end_time = time.time()
    # Compute elapsed time
    elapsed_time = end_time - start_time  
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60
    
    # Print result
    print(f"Time elapsed: {hours} hours, {minutes} minutes, {seconds:.2f} seconds")
    
    #%% Plotting the above field line following
    '''
    Plotting of Field Lines
    '''
    # Plot all elements of all_points with different colors
    fig = plt.figure()
    #ax = fig.add_subplot(projection='polar')
    ax = fig.add_subplot()
    
    # Generate colormap
    #colors = plt.cm.viridis(np.linspace(0, 1, len(surfaces)))  # Using Viridis colormap
    # twilight_shifted, tab20b, and tab10 best so far
    #colors = plt.cm.brg(np.linspace(0, 1, len(surfaces)))  # Using Viridis colormap
    colors = plt.cm.Dark2(np.arange(len(surfaces)) % 8) 

    R = []
    Z = []
    s_values = []
    u_values = []
    plot_data =[]
    #minor_r = []
    #theta   = []
    
    for i, points in enumerate(surfaces):
        # Extract s and u from data
        svals = [sublist[0] for sublist in points]  # Extracting s values
        uvals = [sublist[1] for sublist in points]  # Convert to radians
        
        # Store s and u data
        s_values.append(svals)
        u_values.append(uvals)
        
        # Convert each (s, u, v) to (r, z)
        converted_points = [convert_int(sublist[0], sublist[1], sublist[2]) for sublist in points]
        
        # Extract r and z values for plotting
        r_vals = [pt[0] for pt in converted_points]  # Extract r values
        z_vals = [pt[1] for pt in converted_points]  # Extract z values
        
        # Store R and Z values
        R.append(r_vals)
        Z.append(z_vals)
        
        # Calculate r and theta, need magnetic axis for this part
        #for ind in range(len(r_vals)):
            #r_diff = r_vals[ind] - real_axis[0]
            #z_diff = z_vals[ind] - real_axis[1]
            #disc = (r_diff)**2 + (z_diff)**2
            #temp_theta = np.arctan2(z_diff,r_diff)
            #rad  = np.sqrt(disc)
            #minor_r.append(rad)
            #theta.append(temp_theta)
            #plot_data.append([r_vals[ind], z_vals[ind], svals[ind], uvals[ind], 2*np.pi/nfp*ind])
        color_num= i
        # if (i%8)==3:
        #     color_num = i + 3
        # else:
        #     color_num = i
        m = MarkerStyle('.', joinstyle=JoinStyle('round'))
        ax.scatter(r_vals, z_vals, color=colors[color_num % 8], label=f"s={s[i]:.2f}",s=0.5, marker=m)
        #ax.plot(r_vals, z_vals, color=colors[color_num % 8])

    # Add legend and show plot
    #ax.legend(loc="upper right")
    
    # Figure for the iota plot
    fig_i = plt.figure()
    ax_i = fig_i.add_subplot()
    
    R_i = [0]*len(iota)
    iota_data = [0]*len(iota)
    iota_plot_data = [0]*len(iota)
    
    for i, data in enumerate(iota):
        s_cord = data[0]
        iota_temp = data[1]
        r_i, z_i = convert_int(s_cord,0,0)
        R_i[i] = r_i
        iota_data[i] = iota_temp
        iota_plot_data[i] = [r_i, z_i, iota_temp]
        #ax_i.scatter(r_i,iota_temp)
    ax_i.scatter(R_i,iota_data)
    plt.show()
    
    
    #%%
    """
    Plotting divb
    """
    # Define x and y ranges
    x_vals = np.linspace(0.01, 0.95, 95)
    y_vals = np.linspace(0, 2 * np.pi, 100)
    X, Y = np.meshgrid(x_vals, y_vals)

    # Evaluate Z over the grid with z=0
    Z = np.zeros_like(X)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Z[i, j] = divBint(X[i, j], Y[i, j], 0)

    # Make contour plot with lines
    contours = plt.contour(X, Y, Z, levels=20, cmap='viridis')
    plt.clabel(contours, inline=True, fontsize=8)  # Add labels to contour lines

    # Optional: fill between contour lines
    plt.contourf(X, Y, Z, levels=20, cmap='viridis', alpha=0.6)

    plt.colorbar(label='divB(s, u, 0)')
    plt.xlabel('s')
    plt.ylabel('u')
    plt.title('Contour Plot of divB(s, u, 0) with Labels')
    plt.show()
    
    
    # Define x and y ranges
    x_vals = np.linspace(0.01, 0.95, 95)
    y_vals = np.linspace(0, 2 * np.pi, 100)
    X, Y = np.meshgrid(x_vals, y_vals)
    
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    # Evaluate Z over the grid with z=0
    Z = np.zeros_like(X)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Z[i, j] = divBint(X[i, j], Y[i, j], 0)
            
    # Plot the surface
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.xlabel('s')
    plt.ylabel('u')
    plt.title('3D Plot of divB(s, u, 0)')
    plt.show()
    






















