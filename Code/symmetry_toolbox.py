#=================================================================================
#=================================================================================
# Script:"symmetry_toolbox"
# Date: 2022-07-12
# Implemented by: Johannes Borgqvist
# Description:
# This script contains all the functions for generating the symmetries of the
# separable phase plane models. 
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
#!python
# Numpy for numerical packages and vectors
from numpy import *
# Scipy to do all the computations and integrations
from scipy import integrate
from scipy.optimize import fsolve
from scipy import integrate # For solving ODEs.
# Import the lambert W function from scipy
from scipy.special import lambertw
# Import quad so we can evaluate our integral as well
from scipy.integrate import quad
# Import matplotlib
import matplotlib
import matplotlib.pyplot as plt
#=================================================================================
#=================================================================================
# The Functions
#=================================================================================
#=================================================================================
#==============================================================
# Phase plane symmetries
#==============================================================
# Function 1: dX_dt
# This is a help function defining the ODE system we want to solve.
def dX_dt_LV(X, t=0,a=1):
    """ Return the growth rate of fox and rabbit populations. """
    return array([ X[0]*(1-X[1]) ,
                   a*X[1]*(X[0]-1)])
# Function 2: u_transf
def u_transf(u, epsilon,alpha):
    func = lambda u_hat :  alpha*(u - log(u)) + epsilon - alpha*(u_hat - log(u_hat))
    if abs(u-1)<0.000000000001:
        u_hat_solution = u
    else:
        u_hat_solution = fsolve(func, u)[0]
    return u_hat_solution
# Function 3: v_transf
def v_transf(v, epsilon, alpha):
    func = lambda v_hat :  log(v**(1/alpha)) - (v/alpha) + epsilon - (log(v_hat**(1/alpha)) - (v_hat/alpha))
    if abs(v-1)<0.00000000001:
        v_hat_solution = v
    else:
        v_hat_solution = fsolve(func, v)[0]
    return v_hat_solution
# Function 4: S_transf
def S_transf(S, epsilon, a, r):
    # Incorrectly implemented?
    #func = lambda S_hat :  S - (a/r)*log(S)+epsilon-(S_hat - (a/r)*log(S_hat))
    func = lambda S_hat :  (a/r)*log(S)-S+epsilon-((a/r)*log(S_hat) - S_hat)    
    S_hat_solution = fsolve(func, S)[0]
    return S_hat_solution
# Function 5: IC
def transform_IC_LV(u0, v0, alpha, H):
    func = lambda v:  v+(alpha*u0)-log((u0**alpha)*v) - H
    v0 = fsolve(func, v0)[0]
    return v0
# Function 6: dX_dt_linear
# This is a help function defining the ODE system we want to solve.
def dX_dt_linear(X, t=0,*parameters):
    # Extract the parameters
    a, b, c, d = parameters
    # Return the dynamics of the linear system
    return array([ a*X[0]+b*X[1],
                   c*X[0]+d*X[1]])
# Function 7: Gamma_r_ODE for the linear model
def Gamma_r_ODE(X, epsilon=0,*parameters):
    # Extract the parameters
    a, b, c, d = parameters
    # Return the dynamics of the linear system
    return array([ X[0],X[1]])
                   
# Function 8: Gamma_theta_ODE for the linear model
def Gamma_theta_ODE(X, epsilon=0,*parameters):
    # Extract the parameters
    a, b, c, d = parameters
    # Return the dynamics of the linear system
    return array([ -((c*X[0]**2+(d-a)*X[0]*X[1]-b*X[1]**2)/(a*X[0]**2+(b+c)*X[0]*X[1]+d*X[1]**2))*X[1],
                   ((c*X[0]**2+(d-a)*X[0]*X[1]-b*X[1]**2)/(a*X[0]**2+(b+c)*X[0]*X[1]+d*X[1]**2))*X[0]])
# Function 9: dX_dt_os for the oscillatory model
def dX_dt_os(X, t=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    # Define the radius
    r = sqrt(X[0]**2+X[1]**2)    
    # Return the dynamics of the linear system
    return array([ X[0]*(1-r)-omega*X[1],
                   X[1]*(1-r)+omega*X[0]])
# Function 10: dX_dt_os_polar for the oscillatory model
def dX_dt_os_polar(X, t=0,*parameters):
    # Extract the parameters
    omega = parameters[0]   
    # Return the dynamics of the linear system
    return array([ omega,
                   X[1]*(1-X[1])])
# Function 11: Gamma_r_os for the oscillatory model
def Gamma_r_os(X, epsilon=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    # Define the radius
    r = sqrt(X[0]**2+X[1]**2)
    # Return the dynamics of the biological oscillator
    return array([(1-r)*X[0],(1-r)*X[1]])                   
# Function 12: Gamma_theta_os for the oscillatory model
def Gamma_theta_os(X, epsilon=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    # Return the dynamics of the biological oscillator
    return array([ -(X[1]/omega),(X[0]/omega)])
#==============================================================
# Push forward symmetries
#==============================================================
# Function 12: For our tricky integrand
def integrand_u(s,u_val,v_val,alpha,branch_number):
    # Calculate the internal energy
    H = alpha*u_val+v_val - log((u_val**alpha)*v_val)
    # Define the invariant
    I_s = -((exp(s)/s)**(alpha))*(1/exp(H))
    # Calculate our main factor in the tangent which depends on the Lambert W function 
    factor = 1+lambertw(I_s,branch_number).real
    # Now we can define the denominator
    denom = factor*((s-1)**2)
    # Now, we can return the integral
    return -1/denom
# Function 13: ODE for the u-directional symmetry of the LV model
def dX_deps_LV_u(X, t=0,*parameters):
    # Extract the parameters
    alpha = parameters[0]
    branch_number = parameters[1]
    u_min = parameters[2]
    # Solve the integral for the time tangent
    xi_u = quad(integrand_u, u_min, X[1], args=(X[1],X[2],alpha,branch_number),epsrel = 1e-012)[0]
    # Return the dynamics of the linear system
    return array([xi_u,
                  (1/alpha)*(X[1]/(X[1]-1)),
                 0])
# Function 14: For our tricky integrand
def integrand_v(s,u_val,v_val,alpha,branch_number):
    # Calculate the internal energy
    H = alpha*u_val+v_val - log((u_val**alpha)*v_val)
    # Define the invariant
    I_s = -(1/alpha)*(((exp(s-H))/(alpha*s))**(1/alpha))
    # Calculate our main factor in the tangent which depends on the Lambert W function 
    factor = 1+lambertw(I_s,branch_number).real
    # Now we can define the denominator
    denom = factor*((1-s)**2)
    # Now, we can return the integral
    return -1/denom
# Function 15: ODE for the v-directional symmetry of the LV model
def dX_deps_LV_v(X, t=0,*parameters):
    # Extract the parameters
    alpha = parameters[0]
    branch_number = parameters[1]
    v_min = parameters[2]
    # Solve the integral for the time tangent
    xi_v = quad(integrand_v, v_min, X[2], args=(X[1],X[2],alpha,branch_number),epsrel = 1e-012)[0]
    # Return the dynamics of the linear system
    return array([xi_v,
                  0,
                 ((X[2])/(1-X[2]))])
# Function 16: Gamma_r_os_time for the oscillatory model
def Gamma_r_os_time(X, epsilon=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    H = parameters[1]
    # Calculate the radius
    r = sqrt(X[2]**2+X[1]**2)    
    # Return the dynamics of the biological oscillator
    return array([H, (1-r)*X[1],(1-r)*X[2]])
# Function 16.5: Gamma_r_os_time_polar for the oscillatory model
def Gamma_r_os_time_polar(X, epsilon=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    H = parameters[1]
    # Return the dynamics of the biological oscillator
    return array([H, 0,X[2]*(1-X[2])])                   
# Function 17: Gamma_theta_os for the oscillatory model
def Gamma_theta_os_time(X, epsilon=0,*parameters):
    # Extract the parameters
    omega = parameters[0]
    H = parameters[1]        
    # Return the dynamics of the biological oscillator
    return array([H, -omega*X[2],omega*X[1]])
# Function 18: dX_dt_SIR
# This is a help function defining the ODE system we want to solve.
def dX_dt_SIR(X, t=0,*parameters):
    # Extract the parameters
    a = parameters[0]
    r = parameters[1]
    #Return the ODE s
    return array([ -r*X[0]*X[1] ,
                   r*X[0]*X[1]-a*X[1]])
# Function 19: The integrand for the S-directional symmetry of the
# SIR model
def integrand_S(s,H_SIR,r,p):
    # Now we can define the denominator
    denom = r*((p-s)**2)*(s-p*log(s)-H_SIR)
    # Now, we can return the integral
    return 1/denom
# Function 20: ODE for the S-directional symmetry of the SIR model
def dX_deps_SIR_S(X, t=0,*parameters):
    # Extract the parameters
    H_SIR = parameters[0]
    r = parameters[1]
    p = parameters[2]
    S0 = parameters[3]
    # Solve the integral for the time tangent
    xi_S = quad(integrand_S, S0, X[1], args=(H_SIR,r,p))[0]
    # Return the dynamics of the linear system
    return array([xi_S,
                  X[1]/(p-X[1]),
                  0])
# Function 21: The integrand for the I-directional symmetry of the
# SIR model
def integrand_I(s,H_SIR,r,p,branch_number):
    # Define the thingy we're taking the lambert function of
    I_I = -(exp(s/p)/p)
    #Define our lovely factor
    factor_I = 1+lambertw(I_I,branch_number).real
    # Now we can define the denominator
    denom = r*p*(s**2)*factor_I
    # Now, we can return the integral
    return 1/denom
# Function 22: ODE for the I-directional symmetry of the SIR model
def dX_deps_SIR_I(X, t=0,*parameters):
    # Extract the parameters
    H_SIR = parameters[0]
    r = parameters[1]
    p = parameters[2]
    I0 = parameters[3]
    branch_number = parameters[4]
    # Solve the integral for the time tangent
    xi_I = quad(integrand_I, I0, X[2], args=(H_SIR,r,p,branch_number),epsrel = 1e-012)[0]
    # Return the dynamics of the linear system
    return array([xi_I,
                  0,
                  1])
#==============================================================
# Function 7: plot_LaTeX_2D
# This functions enable us to reproduce our plots using pgfplots in LaTeX
def plot_LaTeX_2D(t,y,file_str,plot_str,legend_str):
    # Open a file with the append option
    # so that we can write to the same
    # file multiple times
    f = open(file_str, "a")
    # Create a temporary string which
    # is the one that does the plotting.
    # Here we incorporate the input plot_str
    # which contains the color, and the markers
    # of the plot at hand
    if len(legend_str)==0:
        temp_str = "\\addplot[\nforget plot,\n" + plot_str+ "\n]\n"
    else:
        temp_str = "\\addplot[\n" + plot_str+ "\n]\n"
    # Add the coordinates
    temp_str += "coordinates {%\n"
    # Loop over the input files and add
    # them to the file
    for i in range(len(t)):
        temp_str += "(" + str(t[i]) + "," + str(y[i]) + ")\n"
    # The plotting is done, let's close the shop    
    temp_str += "};\n"
    # Add a legend if one is provided
    if len(legend_str) > 0:
        temp_str += "\\addlegendentry{" + legend_str + "}\n"
    # Finally, we write the huge string
    # we have created
    f.write("%s"%(temp_str))
    # Close the file
    f.close()
