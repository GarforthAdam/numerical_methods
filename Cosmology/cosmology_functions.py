import numpy as np
import math

#Define all the functions that'll be used
def e(z, O_M, O_L, O_k):
    #pass in vars
    #O-M is matter dens
    #O_L is dark matter dens (density due to cosmo const)
    #O_k is curvature "density"
    #Assuming O_R is zero for this const

    return (O_M*(1+z)**3 + O_L + O_k*(1+z)**2)**0.5

# Define boxes for riemann sum
def zgrid(zmax, ngrid = 1000):
    return np.expm1(np.linspace(0, np.log(1+zmax),ngrid)) #Logarithmically scales grid

#Compute the integral
def sum(z, O_M, O_L, O_k, ngrid=1000):
    zg = zgrid(z,ngrid=ngrid)
    ez = e(zg,O_M, O_L, O_k)

    return np.trapezoid(1/ez, zg)



#Define the curvature function depending on kappa
def S_k(kap, O_M, O_L, z, O_k, DH):
    s = sum(z,O_M,O_L, O_k, ngrid=1000)    
    if(kap > 0):
        #In scenarios where S_kappa is present divide by DH to get DL/DH
        R_0 = (DH) / (abs(1-O_M-O_L)**0.5)
        return ((R_0)*math.sin((s*DH)/R_0)) / DH
    if(kap < 0):
        R_0 = (DH) / (abs(1-O_M-O_L)**0.5)
        return ((R_0)*math.sinh((s*DH)/R_0)) / DH
    else: 
        #In this case just have DH factored out (hence why it isn't in integral function
        return s



'''
Additional functions for chi-squared fitting:
'''
def f(z, O_M, O_L, O_k):
    dc = sum(z,O_M,O_L,O_k, ngrid=1000)
    if (O_k > 0):
        dL = (1 + z) * (1/math.sqrt(abs(O_k))) * math.sinh(math.sqrt(abs(O_k)) * dc)
    elif (O_k < 0):
        dL = (1 + z) * (1/math.sqrt(abs(O_k))) * math.sin(math.sqrt(abs(O_k)) * dc)
    #For other scenario when O_k is 0
    else:
        dL = (1 + z) * dc
    return 5*math.log10(dL)


def chat(zCMB, m_b_corr, m_b_corr_err_DIAG, N, O_M, O_L, O_k):
 #need some sort of n for sum
 #Using for loops to compute the summation
 sum_unc = 0
 for i in range(N):
    sum_unc = sum_unc + (1/(m_b_corr_err_DIAG[i]**2))
 #Inverse it as equation wants
 sum_unc = sum_unc ** (-1)
 #Then the big equation that's just m_i - f func.
 #Just good practice to set the sum to be 0 before using it's loop
 sum_c = 0

 for i in range(N):
     #Defining values to make code a lil bit cleaner
     z = zCMB[i]
        
     f_val = f(z, O_M, O_L, O_k)


     unc = m_b_corr_err_DIAG[i]

     sum_c = sum_c + ((m_b_corr[i] - f_val) / unc**2)
 return sum_unc * sum_c

#Fix optimization issues by computing across slices of the array

def comoving_distance_interp(z_targets, O_M, O_L, O_k, ngrid=2000):
    zmax = np.max(z_targets)
    zg = np.linspace(0, zmax, ngrid)
    integrand = 1.0 / e(zg, O_M, O_L, O_k)
    #Use cumulative sum to compute the integrand, saving memory across computations
    cum = np.concatenate(([0.0],
            np.cumsum((integrand[:-1] + integrand[1:]) / 2 * np.diff(zg)))) 
    return np.interp(z_targets, zg, cum)  # D_C at every SN z, in one shot

#similar f function but now integration is done separately, O_k accounts for python math errors
def f_vec(z, O_M, O_L, O_k, dc):
    if O_k > 1e-8:
        k = np.sqrt(O_k)
        dL = (1+z) * np.sinh(k*dc) / k
    elif O_k < -1e-8:
        k = np.sqrt(-O_k)
        dL = (1+z) * np.sin(k*dc) / k
    else:
        dL = (1+z) * dc
    return 5*np.log10(dL)

#Chi-sq is put into a separate function instead of calculated in the pipeline
def chisq_vec(zCMB, m_b_corr, m_b_corr_err, O_M, O_L, O_k):
    dc = comoving_distance_interp(zCMB, O_M, O_L, O_k)
    f_vals = f_vec(zCMB, O_M, O_L, O_k, dc)
    w = 1.0 / m_b_corr_err**2
    C = np.sum(w * (m_b_corr - f_vals)) / np.sum(w)
    resid = m_b_corr - f_vals - C
    return np.sum((resid / m_b_corr_err)**2)