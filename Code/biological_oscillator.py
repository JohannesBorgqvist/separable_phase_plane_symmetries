#=================================================================================
#=================================================================================
# Script:"biological_oscillator"
# Date: 2022-12-20
# Implemented by: Johannes Borgqvist
# Description:
# The script plots the phase plane symmetries as well as the corresponding symmetries
# of the biological oscillator in the time domain.
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
from symmetry_toolbox import * # Script containing all functions
# Matplotlib for plotting
import matplotlib.pyplot as plt # For plotting
# Set all parameters to tex
plt.rcParams['text.usetex'] = True
#=================================================================================
#=================================================================================
# Phase plane symmetries of the biological oscillator
#=================================================================================
#=================================================================================
# ---------------------------------------------------------------------------
# Overall properties
# ---------------------------------------------------------------------------
# Transformation parameters
epsilon_r=0.33
epsilon_theta = pi/6
# The indices we wish to transform
lin_indices = arange(0,60,6)
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_r-0.005,epsilon_r/50)
# The frequency rate
omega = 1
# Define the time vector
t = linspace(0, 10, 500)              # Time
X0 = array([2, 2])                  # ICs
X1_os, infodict = integrate.odeint(dX_dt_os, X0, t, args = ((omega),),full_output=True)
# Extract the original solution with the defined parameters and initial conditions
u_3, v_3 = X1_os.T
# Allocate our empty lists
Gamma_os_r_u_1 = []
Gamma_os_r_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_os, array([u_3[lin_index],v_3[lin_index]]), epsilon_vec, args = ((omega),),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_os_r_u_temp, Gamma_os_r_v_temp = X_r.T
    # Save our transformations
    Gamma_os_r_u_1.append(Gamma_os_r_u_temp)
    Gamma_os_r_v_1.append(Gamma_os_r_v_temp)
# Solve to get the transformed angle solution
X1_r, infodict = integrate.odeint(dX_dt_os, array([Gamma_os_r_u_1[0][-1], Gamma_os_r_v_1[0][-1]]), t, args = ((omega,)),full_output=True)    
# Extract the solution
u_os_r_2, v_os_r_2 = X1_r.T
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_theta-0.005,epsilon_theta/50)
# Allocate our empty lists
Gamma_os_theta_u_1 = []
Gamma_os_theta_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_os, array([u_3[lin_index],v_3[lin_index]]), epsilon_vec, args = ((omega),),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_os_theta_u_temp, Gamma_os_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_os_theta_u_1.append(Gamma_os_theta_u_temp)
    Gamma_os_theta_v_1.append(Gamma_os_theta_v_temp)
# Solve to get the transformed angle solution
X1_theta, infodict = integrate.odeint(dX_dt_os, array([Gamma_os_theta_u_1[0][-1], Gamma_os_theta_v_1[0][-1]]), t, args = ((omega,)),full_output=True)    
# Extract the solution
u_os_theta_2, v_os_theta_2 = X1_theta.T
#=================================================================================
#=================================================================================
# Plot the symmetries in the phase plane
#=================================================================================
#=================================================================================
#Define the first figure
f1, ax_1 = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
# Plot 1: Radial symmetry on the biological oscillator
ax_1[0].plot(u_3,v_3,color=(0/256,68/256,27/256),label="$(u,v)$",linewidth=4.0)
ax_1[0].plot(u_os_r_2,v_os_r_2,color=(0/256,109/256,44/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
for index in range(len(Gamma_os_r_u_1)):
    if index == 0:
        ax_1[0].plot(Gamma_os_r_u_1[index], Gamma_os_r_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{r}$",linewidth=0.5)
    else:
        ax_1[0].plot(Gamma_os_r_u_1[index], Gamma_os_r_v_1[index],'--',color=(0,0,0),linewidth=0.5)
ax_1[0].grid()
ax_1[0].legend(loc='best',prop={"size":20})
ax_1[0].set_xlabel(xlabel="$u{(t)}$",fontsize=25)
ax_1[0].set_ylabel(ylabel="$v{(t)}$",fontsize=25)
ax_1[0].tick_params(axis='both', which='major', labelsize=20)
ax_1[0].tick_params(axis='both', which='minor', labelsize=20)
# Plot 2: Angular symmetry on the biological oscillator
ax_1[1].plot(u_3,v_3,color=(103/256,0/256,31/256),label="$(u,v)$",linewidth=4.0)
ax_1[1].plot(u_os_theta_2,v_os_theta_2,color=(206/256,18/256,86/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
for index in range(len(Gamma_os_theta_u_1)):
    if index == 0:
        ax_1[1].plot(Gamma_os_theta_u_1[index], Gamma_os_theta_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{\\theta}$",linewidth=0.5)
    else:
        ax_1[1].plot(Gamma_os_theta_u_1[index], Gamma_os_theta_v_1[index],'--',color=(0,0,0),linewidth=0.5)
ax_1[1].grid()
ax_1[1].legend(loc='best',prop={"size":20})
ax_1[1].set_xlabel(xlabel="$u{(t)}$",fontsize=25)
ax_1[1].set_ylabel(ylabel="$v{(t)}$",fontsize=25)
ax_1[1].tick_params(axis='both', which='major', labelsize=20)
ax_1[1].tick_params(axis='both', which='minor', labelsize=20)
# We have a title of this figure as well
f1.suptitle('Phase plane symmetries of the biological oscillator',fontsize=30,weight='bold')
f1.savefig('../Figures/phase_plane_symmetries_biological_oscillator.png')
#=================================================================================
#=================================================================================
# Time domain symmetries of the biological oscillator
#=================================================================================
#=================================================================================
# ---------------------------------------------------------------------------
# Overall properties
# ---------------------------------------------------------------------------
# The frequency rate
omega = 1
# Define the time vector
t = linspace(0, 6, 500)              # Time
X0 = array([2, 2])                  # ICs
# Calculate the internal energy based on these initial conditions
# Define the radius
r = sqrt(X0[1]**2+X0[0]**2)
# Define the angle
theta = arctan(X0[1]/X0[0])
# Define the internal energy
H = ((r)/(abs(1-r)))*exp(-((theta)/(omega)))
# Solve the system of ODEs
X1_os, infodict = integrate.odeint(dX_dt_os_polar, array([theta,r]), t, args = (omega,),full_output=True)
# Extract the original solution with the defined parameters and initial conditions
#u_os, v_os = X1_os.T
theta_os, r_os = X1_os.T
u_os = r_os*cos(theta_os)
v_os = r_os*sin(theta_os)
# Transformation parameters
epsilon_r=1
epsilon_theta = pi/2
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_r,epsilon_r/100)
#=================================================================================
# The indices we wish to transform
lin_indices = concatenate([arange(0,50,14),arange(50,450,50)])
# Allocate our empty lists
Gamma_os_r_t_1 = []
Gamma_os_r_u_1 = []
Gamma_os_r_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_os_time, array([t[lin_index],u_os[lin_index],v_os[lin_index]]), epsilon_vec, args = ((omega),H),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'    
    # Extract transformations
    Gamma_os_r_t_temp,Gamma_os_r_u_temp, Gamma_os_r_v_temp = X_r.T
    # Save our transformations
    Gamma_os_r_t_1.append(Gamma_os_r_t_temp)
    Gamma_os_r_u_1.append(Gamma_os_r_u_temp)
    Gamma_os_r_v_1.append(Gamma_os_r_v_temp)
# The transformed time as well
t_transformed_r = linspace(Gamma_os_r_t_1[0][-1],t[-1],len(t))
# Solve to get the transformed angle solution
X1_r, infodict = integrate.odeint(dX_dt_os, array([Gamma_os_r_u_1[0][-1], Gamma_os_r_v_1[0][-1]]), t_transformed_r, args = ((omega),),full_output=True)    
# Extract the solution
u_os_r_2, v_os_r_2 = X1_r.T
#--------------------------------------------------------------------------------------------------
# The indices we wish to transform
lin_indices = concatenate([arange(0,50,14),arange(50,450,50)])
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_theta,epsilon_theta/50)
# Allocate our empty lists
Gamma_os_theta_t_1 = []
Gamma_os_theta_u_1 = []
Gamma_os_theta_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_os_time, array([t[lin_index],u_os[lin_index],v_os[lin_index]]), epsilon_vec, args = ((omega),H),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_os_theta_t_temp,Gamma_os_theta_u_temp, Gamma_os_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_os_theta_t_1.append(Gamma_os_theta_t_temp)
    Gamma_os_theta_u_1.append(Gamma_os_theta_u_temp)
    Gamma_os_theta_v_1.append(Gamma_os_theta_v_temp)
# The transformed time as well
t_transformed_theta = linspace(Gamma_os_theta_t_1[0][-1],t[-1],len(t))
# Solve to get the transformed angle solution
X1_theta, infodict = integrate.odeint(dX_dt_os, array([Gamma_os_theta_u_1[0][-1], Gamma_os_theta_v_1[0][-1]]), t_transformed_theta, args = ((omega,H)),full_output=True)    
# Extract the solution
u_os_theta_2, v_os_theta_2 = X1_theta.T
#=================================================================================
#=================================================================================
# Plot the symmetries in the time domain
#=================================================================================
#=================================================================================
#Define the first figure
f2, ax_2 = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
# Plot 1: Radial symmetry on the biological oscillator
ax_2[0].plot(t,u_os,color=(0/256,68/256,27/256),label="$u(t)$",linewidth=4.0)
ax_2[0].plot(t_transformed_r,u_os_r_2,color=(35/256,139/256,69/256),label="$\\hat{u}(t)$",linewidth=4.0)
ax_2[0].plot(t,v_os,color=(103/256,0/256,31/256),label="$v(t)$",linewidth=4.0)
ax_2[0].plot(t_transformed_r,v_os_r_2,color=(206/256,18/256,86/256),label="$\\hat{v}(t)$",linewidth=4.0)
# Plot the symmetry
for index in range(len(Gamma_os_r_u_1)):
    if index == 0:
        ax_2[0].plot(Gamma_os_r_t_1[index], Gamma_os_r_u_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{r}$",linewidth=2.5)
        ax_2[0].plot(Gamma_os_r_t_1[index], Gamma_os_r_v_1[index],'--',color=(0,0,0),linewidth=2.5)        
    else:
        ax_2[0].plot(Gamma_os_r_t_1[index], Gamma_os_r_u_1[index],'--',color=(0,0,0),linewidth=2.5)
        ax_2[0].plot(Gamma_os_r_t_1[index], Gamma_os_r_v_1[index],'--',color=(0,0,0),linewidth=2.5)   
ax_2[0].grid()
ax_2[0].legend(loc='best',prop={"size":20})
ax_2[0].set_xlabel(xlabel="Time, $t$",fontsize=25)
ax_2[0].set_ylabel(ylabel="States",fontsize=25)
ax_2[0].tick_params(axis='both', which='major', labelsize=20)
ax_2[0].tick_params(axis='both', which='minor', labelsize=20)
# Plot 2: Angular symmetry on the biological oscillator
ax_2[1].plot(t,u_os,color=(0/256,68/256,27/256),label="$u(t)$",linewidth=4.0)
ax_2[1].plot(t_transformed_theta,u_os_theta_2,color=(35/256,139/256,69/256),label="$\\hat{u}(t)$",linewidth=4.0)
ax_2[1].plot(t,v_os,color=(103/256,0/256,31/256),label="$v(t)$",linewidth=4.0)
ax_2[1].plot(t_transformed_theta,v_os_theta_2,color=(206/256,18/256,86/256),label="$\\hat{v}(t)$",linewidth=4.0)
# Plot the symmetry
for index in range(len(Gamma_os_theta_u_1)):
    if index == 0:
        ax_2[1].plot(Gamma_os_theta_t_1[index], Gamma_os_theta_u_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{\\theta}$",linewidth=2.5)
        ax_2[1].plot(Gamma_os_theta_t_1[index], Gamma_os_theta_v_1[index],'--',color=(0,0,0),linewidth=2.5)        
    else:
        ax_2[1].plot(Gamma_os_theta_t_1[index], Gamma_os_theta_u_1[index],'--',color=(0,0,0),linewidth=2.5)
        ax_2[1].plot(Gamma_os_theta_t_1[index], Gamma_os_theta_v_1[index],'--',color=(0,0,0),linewidth=2.5)   
ax_2[1].grid()
ax_2[1].legend(loc='best',prop={"size":20})
ax_2[1].set_xlabel(xlabel="Time, $t$",fontsize=25)
ax_2[1].set_ylabel(ylabel="States",fontsize=25)
ax_2[1].tick_params(axis='both', which='major', labelsize=20)
ax_2[1].tick_params(axis='both', which='minor', labelsize=20)
# We have a title of this figure as well
f2.suptitle('Symmetries of the biological oscillator',fontsize=30,weight='bold')
f2.savefig('../Figures/time_domain_symmetries_biological_oscillator.png')
# Show the plot in the end
plt.show()


