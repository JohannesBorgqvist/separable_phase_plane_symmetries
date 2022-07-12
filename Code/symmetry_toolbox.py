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
from numpy import *
import pylab as p
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import fsolve
#=================================================================================
#=================================================================================
# The Functions
#=================================================================================
#=================================================================================
# Function 1: dX_dt
# This is a help function defining the ODE system we want to solve.
def dX_dt_LV(X, t=0,a=1):
    """ Return the growth rate of fox and rabbit populations. """
    return array([ X[0]*(1-X[1]) ,
                   a*X[1]*(X[0]-1)])
# Function 2: u_transf
def u_transf(u, epsilon):
    func = lambda u_hat :  u - log(u) + epsilon - (u_hat - log(u_hat))
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
    #func = lambda S_hat :  S - (1/r)*log(a-r*S)+epsilon-S_hat+(1/r)*log(a-r*S_hat)
    func = lambda S_hat :  S - (a/r)*log(S)+epsilon-(S_hat - (a/r)*log(S_hat))
    S_hat_solution = fsolve(func, S)[0]
    return S_hat_solution
# Function 5: IC
def transform_IC_LV(u0, v0, alpha, H):
    func = lambda v:  v+(alpha*u0)-log((u0**alpha)*v) - H
    v0 = fsolve(func, v0)[0]
    return v0
