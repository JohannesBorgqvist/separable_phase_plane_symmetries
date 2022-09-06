#=================================================================================
#=================================================================================
# Script:"plot_phase_plane_symmetries"
# Date: 2022-07-12
# Implemented by: Johannes Borgqvist
# Description:
# The script plots all the phase plane symmetries of the LV-model, the SIR-model,
# the linear model and the biological oscillator.
#=================================================================================
#=================================================================================
#=================================================================================
#=================================================================================
# Import Libraries
#=================================================================================
#=================================================================================
# Home made all the built in function
from symmetry_toolbox import *
# Matplotlib for plotting
import matplotlib.pyplot as plt # For plotting,
#=================================================================================
#=================================================================================
# Plot the solutions
#=================================================================================
#=================================================================================
# Define the time vector and the initial conditions
t = linspace(0, 10, 500)              # time
X0 = array([1, 0.10])                     # initials conditions: 10 rabbits and 5 foxes
# Define our parameter alpha
alpha = 1
# Solve the ODE at hand
X1, infodict = integrate.odeint(dX_dt_LV, X0, t, args = (alpha,),full_output=True)
infodict['message'] # >>> 'Integration successful.'
# Split the solution into its component parts
u, v = X1.T
# Plot our lovely solutions
fig_1 = plt.figure(constrained_layout=True, figsize=(20, 8))
plt.plot(t, u, '-', label="Prey, $u(\\tau)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
plt.plot(t, v  , '-', label='Predator, $v(\\tau)$',color=(77/256,0/256,75/256),linewidth=3.0)
plt.grid()
plt.legend(loc='best',prop={"size":20})
plt.xlabel(xlabel='Time, $\\tau$',fontsize=25)
plt.ylabel(ylabel='Population size',fontsize=25)
# Change the size of the ticks
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
plt.title('Solutions of the Lotka-Volterra model',fontsize=30,weight='bold')
#plt.savefig('../Figures/LV_solutions.png')
#=================================================================================
#=================================================================================
# Plot the action of the u-directional symmetry
#=================================================================================
#=================================================================================
# Set the value of alpha
alpha = 1
# Epsilon value
epsilon = 0.5
# The transformation parameter
epsilon_vec = linspace(0,epsilon,200)              # epsilon
# We know that Lambertz w satisfies the following (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lambertw.html):
# The principal branch (k=0) is real if z>-1/e,
# The branch (k>1) is real if -1/e < z < 0.
# Let's make an experiment.
# Magical index
magical_index = 5
# Calculate the minimal u value
u_min = 1
# Take a point
X0 = array([t[magical_index], u[magical_index], v[magical_index]])
# Try to solve the ODE at hand
Gamma_epsilon, infodict = integrate.odeint(dX_deps_LV_u, X0, epsilon_vec, args = (alpha,0,u_min),full_output=True)
# Split the solution into its component parts
Gamma_u_t, Gamma_u_u, Gamma_u_v = Gamma_epsilon.T
# Define a new time vector 
t_2 = linspace(Gamma_u_t[-1], 10, 500)
# We need to integrate backwards in time as well to start at 0
t_2_start = linspace(Gamma_u_t[-1], 0, 10)
# Define new initial conditions for the transformed solutions
X02 = array([Gamma_u_t[-1], Gamma_u_u[-1], Gamma_u_v[-1]])  
X0_2 = array([Gamma_u_u[-1], Gamma_u_v[-1]])  
# Solve the ODE at hand with the new initial conditions
# Integrate backwards to get to 0
X2_start, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2_start, args = (1,),full_output=True)
# Find the rest of the solution
X2, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2, args = (1,),full_output=True)
#infodict['message'] # >>> 'Integration successful.'
# Split the solution into its component parts
u_2, v_2 = X2.T
# Split the start bit as well
u_2_start, v_2_start = X2_start.T
# Concatenate to get the full solution curves and the full time
u_2 = concatenate((flip(u_2_start,0), u_2), axis=0)
v_2 = concatenate((flip(v_2_start,0), v_2), axis=0)
t_2 = concatenate((flip(t_2_start,0), t_2), axis=0)
# Plot the symmetry again for some other point on the solution curves
magical_indices = [50, 60, 70, 80, 85, 90, 95, 100, 105, 110]
# The corresponding branches
branch_indices = [0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1]
# Allocate memory for our lovely symmetry
Gamma_u_t_vec = []
Gamma_u_u_vec = []
# Loop over the indices and plot the symmetry transformation
for branch_index,new_magical_index in enumerate(magical_indices):
    # Take a point
    X0 = array([t[new_magical_index], u[new_magical_index], v[new_magical_index]])  
    # Try to solve the ODE at hand
    Gamma_epsilon_temp, infodict = integrate.odeint(dX_deps_LV_u, X0, epsilon_vec, args = (alpha,branch_indices[branch_index],u_min),full_output=True)    
    # Split the solution into its component parts
    Gamma_u_t_temp, Gamma_u_u_temp, Gamma_u_v_temp = Gamma_epsilon_temp.T    
    # Append our solutions
    Gamma_u_t_vec.append(Gamma_u_t_temp)    
    Gamma_u_u_vec.append(Gamma_u_u_temp)    
#-------------------------------------------------------------------------------------------------------------------------------
# Se if we can plot our symmetry?
fig_2 = plt.figure(constrained_layout=True, figsize=(20, 8))
# The original solution
plt.plot(t, u, '-', label="Prey, $u(\\tau)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
# The transformed solution
plt.plot(t_2, u_2, '-', label="Prey, $\\hat{u}(\\hat{\\tau})$" ,color=(35/256,139/256,69/256),linewidth=3.0)
# The symmetry
for index in range(len(Gamma_u_t_vec)):
    if index == 0:
        plt.plot(Gamma_u_t_vec[index],Gamma_u_u_vec[index], '--', label="Symmetry, $\\Gamma_{\\epsilon}^{\\mathrm{LV},u}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
    else:
        plt.plot(Gamma_u_t_vec[index],Gamma_u_u_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
plt.grid()
plt.legend(loc='best',prop={"size":20})
plt.xlabel(xlabel='Time, $\\tau$',fontsize=25)
plt.ylabel(ylabel='Population size',fontsize=25)
# Change the size of the ticks
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
plt.title('Symmetry $\\Gamma_{\epsilon}^{\\mathrm{LV},u}$',fontsize=30,weight='bold')
#plt.show()
#=================================================================================
#=================================================================================
# Plot the action of the v-directional symmetry
#=================================================================================
#=================================================================================
# Epsilon value
epsilon = 0.5
# The transformation parameter
epsilon_vec = linspace(0,epsilon,200)              # epsilon
# Magical index
magical_index = 130
# Print the value of v
print(v[magical_index])
# Take a point
X0 = array([t[magical_index], u[magical_index], v[magical_index]])
# Define the minimal u value
v_min = 1
# Try to solve the ODE at hand
Gamma_epsilon, infodict = integrate.odeint(dX_deps_LV_v, X0, epsilon_vec, args = (alpha,0,v_min),full_output=True)
# Split the solution into its component parts
Gamma_v_t, Gamma_v_u, Gamma_v_v = Gamma_epsilon.T
# Define a new time vector 
t_2 = linspace(Gamma_v_t[-1], 10, 100)
# We need to integrate backwards in time as well to start at 0
t_2_start = linspace(Gamma_v_t[-1], 0, 400)
# Define new initial conditions for the transformed solutions
X02 = array([Gamma_v_t[-1], Gamma_v_u[-1], Gamma_v_v[-1]])  
X0_2 = array([Gamma_v_u[-1], Gamma_v_v[-1]])  
# Solve the ODE at hand with the new initial conditions
# Integrate backwards to get to 0
X2_start, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2_start, args = (1,),full_output=True)
# Find the rest of the solution
X2, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2, args = (1,),full_output=True)
#infodict['message'] # >>> 'Integration successful.'
# Split the solution into its component parts
u_2, v_2 = X2.T
# Split the start bit as well
u_2_start, v_2_start = X2_start.T
# Concatenate to get the full solution curves and the full time
u_2 = concatenate((flip(u_2_start,0), u_2), axis=0)
v_2 = concatenate((flip(v_2_start,0), v_2), axis=0)
t_2 = concatenate((flip(t_2_start,0), t_2), axis=0)
# Plot the symmetry again for some other point on the solution curves
magical_indices = [magical_index-15,magical_index-10,magical_index-5,magical_index,magical_index+5,magical_index+10,magical_index+15]
# The corresponding branches
branch_indices = [-1,-1,-1,0,0,0,0,0]
# Allocate memory for our lovely symmetry
Gamma_v_t_vec = []
Gamma_v_v_vec = []
# Loop over the indices and plot the symmetry transformation
for branch_index,new_magical_index in enumerate(magical_indices):
    # Take a point
    X0 = array([t[new_magical_index], u[new_magical_index], v[new_magical_index]])  
    # Try to solve the ODE at hand
    Gamma_epsilon_temp, infodict = integrate.odeint(dX_deps_LV_v, X0, epsilon_vec, args = (alpha,branch_indices[branch_index],v_min),full_output=True)    
    # Split the solution into its component parts
    Gamma_v_t_temp, Gamma_v_u_temp, Gamma_v_v_temp = Gamma_epsilon_temp.T    
    # Append our solutions
    Gamma_v_t_vec.append(Gamma_v_t_temp)    
    Gamma_v_v_vec.append(Gamma_v_v_temp)  
#-------------------------------------------------------------------------------------------------------------------------------
# Se if we can plot our symmetry?
fig_3 = plt.figure(constrained_layout=True, figsize=(20, 8))
# The original solution
plt.plot(t, v, '-', label="Predator, $v(\\tau)$" ,color=(77/256,0/256,75/256),linewidth=3.0)
# The transformed solution
plt.plot(t_2, v_2, '-', label="Predator, $\\hat{v}(\\hat{\\tau})$" ,color=(136/256,65/256,157/256),linewidth=3.0)
# Plot the symmetry
for index in range(len(Gamma_v_t_vec)):
    if index == 0:
        plt.plot(Gamma_v_t_vec[index],Gamma_v_v_vec[index], '--', label="Symmetry, $\\Gamma_{\\epsilon}^{\\mathrm{LV},v}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
    else:
        plt.plot(Gamma_v_t_vec[index],Gamma_v_v_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
plt.grid()
plt.legend(loc='best',prop={"size":20})
plt.xlabel(xlabel='Time, $\\tau$',fontsize=25)
plt.ylabel(ylabel='Population size',fontsize=25)
# Change the size of the ticks
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
plt.title('Symmetry $\\Gamma_{\epsilon}^{\\mathrm{LV},v}$',fontsize=30,weight='bold')
plt.show()
