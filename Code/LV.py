# =================================================================================
# =================================================================================
# Script:"LV"
# Date: 2022-12-20
# Implemented by: Johannes Borgqvist
# Description:
# The script plots the phase plane symmetries as well as the corresponding symmetries of the Lotka-Volterra model in the time domain.
# =================================================================================
# =================================================================================
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
from symmetry_toolbox import *  # Script containing all functions
# Matplotlib for plotting
import matplotlib.pyplot as plt  # For plotting
# Set all parameters to tex
plt.rcParams['text.usetex'] = True
# =================================================================================
# =================================================================================
# Solving the ODEs for the LV model
# =================================================================================
# =================================================================================
# Dimensionless LV parameter
alpha = 1
t = linspace(0, 10, 500)              # time
# initials conditions: 10 rabbits and 5 foxes
X0 = array([1.00, 0.10])
X1, infodict = integrate.odeint(dX_dt_LV, X0, t, args=(alpha,), full_output=True)
infodict['message']                     # >>> 'Integration successful.'
#!python
u, v = X1.T
# Testing to find ICs
# Calculate internal energy
H = v[0] + alpha*u[0] - log((u[0]**alpha)*v[0])
# Epsilon value
epsilon = 0.5
# The transformation parameter
epsilon_vec = linspace(0,epsilon,200)

# Define the u-coordinate in the new initial condition
u0 = X0[0]
v0 = X0[0]
# Find the other v-coordinate in the transformed initial condition
v0_new = transform_IC_LV(u0, v0, alpha, H-alpha*epsilon)
v0_new_2 = transform_IC_LV(u0, v0, alpha, H+alpha*epsilon)
# Try to solve the LV-model again with these new ICs
X2, infodict = integrate.odeint(dX_dt_LV, array(
    [u0, v0_new]), t, args=(alpha,), full_output=True)
u_transformed, v_transformed = X2.T
X3, infodict = integrate.odeint(dX_dt_LV, array(
    [u0, v0_new_2]), t, args=(alpha,), full_output=True)
u_transformed_2, v_transformed_2 = X3.T
# Transformed solutions
t_trans = asarray([t_temp + epsilon for t_temp in list(t)])
# Transformations time translation
t_sym = []
epsilon_vec = arange(0, epsilon, epsilon/50)
t_indices = list(arange(60, 78, 3))
for t_index in list(t_indices):
    trans_vec = [t[t_index]+epsilon_temp for epsilon_temp in list(epsilon_vec)]
    t_sym.append(trans_vec)
# Transformations u
u_sym = []
u_indices = [40, 60, 75, 90, 100, 110, 117, 130, 134, 145, 150, 160, 290, 307]
#u_indices = [134]
for u_index in list(u_indices):
    #trans_vec = [u_transf(u[u_index], epsilon_temp, alpha) for epsilon_temp in list(epsilon_vec)]
    X, infodict = integrate.odeint(u_transf_ODE, array([u[u_index], v[u_index]]), epsilon_vec, args = (alpha,),full_output=True)
    trans_vec, v_trans = X.T    
    u_sym.append(list(trans_vec))
# Transformations v
v_sym = []
v_indices = [30, 50, 64, 106, 112, 117,125, 140, 155, 170, 300, 350, 400, -1]
for v_index in list(v_indices):
    trans_vec = [v_transf(v[v_index], epsilon_temp, alpha)
                 for epsilon_temp in list(epsilon_vec)]
    v_sym.append(trans_vec)
# =================================================================================
# =================================================================================
# Plotting the solutions
# =================================================================================
# =================================================================================
# Define the first figure
f1, ax_1 = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
# u-directional symmetry
ax_1[0].plot(u, v, '-', label='Original population, $(u,v)$',
             color=(0/256, 68/256, 27/256), linewidth=3.0)
ax_1[0].plot(u_transformed_2, v_transformed_2, '-',
             label='Transformed population, $(\\hat{u},v)$', color=(35/256,139/256,69/256), linewidth=3.0)
ax_1[0].plot(asarray(u_sym[0]), asarray([v[u_indices[0]]*(index+1)/(index+1) for index in range(len(epsilon_vec))]), '--',
             label="$\\Gamma_{2,\\epsilon}^{\\mathrm{LV},u}$", color=(0, 0, 0), linewidth=2.0)
for index, u_index in enumerate(list(u_indices)):
    ax_1[0].plot(asarray(u_sym[index]), asarray([v[u_index]*((index+1)/(index+1))
                                                 for index in range(len(epsilon_vec))]), '--', color=(0, 0, 0), linewidth=2.0)
ax_1[0].grid()
ax_1[0].legend(loc='best', prop={"size": 20})
# v-directional symmetry
ax_1[1].plot(u, v, '-', label='Original population, $(u,v)$',
             color=(0/256, 68/256, 27/256), linewidth=3.0)
ax_1[1].plot(u_transformed, v_transformed, '-',
             label='Transformed population, $(u,\\hat{v})$', color=(35/256,139/256,69/256), linewidth=3.0)
ax_1[1].plot(asarray([u[v_indices[0]]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[0]),
             '--', label="$\\Gamma_{2,\\epsilon}^{\\mathrm{LV},v}$", color=(0, 0, 0), linewidth=2.0)
for index, v_index in enumerate(list(v_indices)):
    ax_1[1].plot(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),
                 asarray(v_sym[index]), '--', color=(0, 0, 0), linewidth=2.0)
ax_1[1].grid()
ax_1[1].legend(loc='best', prop={"size": 20})
#ax_1[1].legend(loc='best', prop={"size": 20})
# Set fontsize of labels
ax_1[0].set_xlabel(xlabel='Rabbits, $u(t)$', fontsize=25)
ax_1[0].set_ylabel(ylabel='Foxes, $v(t)$', fontsize=25)
ax_1[1].set_xlabel(xlabel='Rabbits, $u(t)$', fontsize=25)
ax_1[1].set_ylabel(ylabel='Foxes, $v(t)$', fontsize=25)
# Change the size of the ticks
ax_1[0].tick_params(axis='both', which='major', labelsize=20)
ax_1[0].tick_params(axis='both', which='minor', labelsize=20)
ax_1[1].tick_params(axis='both', which='major', labelsize=20)
ax_1[1].tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
f1.suptitle('Phase plane symmetries of the Lotka-Volterra model',
            fontsize=30, weight='bold')
f1.savefig('../Figures/phase_plane_symmetries_LV.png')
#plt.show()
# =================================================================================
# =================================================================================
# Plot phase plane symmetries in LaTeX as well...
# =================================================================================
# =================================================================================
# --------------------------------------------------------------------------------
# u-directional
# --------------------------------------------------------------------------------
# Plot the symmetries in the phase plane
for index, u_index in enumerate(list(u_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray(u_sym[index]), asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex", "color=black,->,>=latex,densely dashed,line width=1.0pt", "$\\Gamma^{\\mathrm{LV},u}_{2,\\epsilon}$")
    else:
        plot_LaTeX_2D(asarray(u_sym[index]), asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex", "color=black,->,>=latex,densely dashed,line width=1.0pt", [])
# Plot the solutions in the phase plane
plot_LaTeX_2D(u, v, "../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex","color=phase_1,line width=1.5pt,", "$(u,v)$")
plot_LaTeX_2D(u_transformed_2, v_transformed_2, "../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex","color=phase_2,line width=1.5pt,", "$(\\hat{u},\\hat{v})$")

# --------------------------------------------------------------------------------
# v-directional
# --------------------------------------------------------------------------------
# Plot the symmetries in the phase plane
for index, v_index in enumerate(list(v_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[index]),
                      "../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex", "color=black,->,>=latex,densely dashed,line width=1.0pt", "$\\Gamma^{\\mathrm{LV},v}_{2,\\epsilon}$")
    else:
        plot_LaTeX_2D(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[index]),
                      "../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex", "color=black,->,>=latex,densely dashed,line width=1.0pt", [])
# Plot the solutions in the phase plane
plot_LaTeX_2D(u, v, "../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex","color=phase_1,line width=1.5pt,", "$(u,v)$")
plot_LaTeX_2D(u_transformed, v_transformed, "../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex","color=phase_2,line width=1.5pt,", "$(\\hat{u},\\hat{v})$")
# =================================================================================
# =================================================================================
# Plot the lifted symmetries as well
# =================================================================================
# =================================================================================
# Re-define the original solution curve
# Dimensionless LV parameter
t = linspace(0, 10, 1000)              # time
# initials conditions: 10 rabbits and 5 foxes
X0 = array([1.00, 0.10])
X1, infodict = integrate.odeint(dX_dt_LV, X0, t, args=(alpha,), full_output=True)
infodict['message']                     # >>> 'Integration successful.'
#!python
u, v = X1.T
# Testing to find ICs
# Calculate internal energy
H = v[0] + alpha*u[0] - log((u[0]**alpha)*v[0])
# The transformation parameter
epsilon_vec = linspace(0,epsilon,200)
epsilon_vec_dense = linspace(0,epsilon,5000)
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# u-directional symmetry of the LV model
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# epsilon
# We know that Lambertz w satisfies the following (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lambertw.html):
# The principal branch (k=0) is real if z>-1/e,
# The branch (k>1) is real if -1/e < z < 0.
# Let's make an experiment.
# Magical index
magical_index = 200
# Take a point on the original solution curve which we transform initially
X0 = array([t[magical_index], u[magical_index], v[magical_index]])
# Define our lower limit
# Calculate the max and min values
umin = -lambertw(-exp((1-H)/alpha),0).real
umax = -lambertw(-exp((1-H)/alpha),-1).real
u0 = umin

# Loop over the dense epsilon vector and iteratively solve our ODE system for the symmetry
# Calculate our lovely branch index
if v[magical_index]>1:
    branch_index = -1
else:
    branch_index = 0
# Try to solve the ODE at hand
Gamma_epsilon, infodict = integrate.odeint(dX_deps_LV_u, X0, epsilon_vec_dense, args = (alpha,branch_index),full_output=True)
# Split the solution into its component parts
Gamma_u_t, Gamma_u_u, Gamma_u_v = Gamma_epsilon.T
# Define a new time vector 
t_2 = linspace(Gamma_u_t[-1], t[-1], 20000)
# We need to integrate backwards in time as well to start at 0
t_2_start = linspace(Gamma_u_t[-1], 0, 20000)
# Define new initial conditions for the transformed solutions
X0_2 = array([Gamma_u_u[-1], Gamma_u_v[-1]])  
# Solve the ODE at hand with the new initial conditions
# Integrate backwards to get to 0
X2_start, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2_start, args = (alpha,),full_output=True)
# Find the rest of the solution
X2, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2, args = (alpha,),full_output=True)
#infodict['message'] # >>> 'Integration successful.'
# Split the solution into its component parts
u_2, v_2 = X2.T
# Split the start bit as well
u_2_start, v_2_start = X2_start.T
# Concatenate to get the full solution curves and the full time
u_2_u = concatenate((flip(u_2_start,0), u_2), axis=0)
v_2_u = concatenate((flip(v_2_start,0), v_2), axis=0)
t_2_u = concatenate((flip(t_2_start,0), t_2), axis=0)
# Plot the symmetry again for some other point on the solution curves
# #---------------------------------------------------------------------------------
# Full design
magical_indices = [25, 75, 125, 175, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, 600, 650, 700, 775, 800, 850, 875, 900]
delta_t_forward_vec = [2, 2, 2, 2, 2, 0.2, 0, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0.2, 0.2, 2, 2, 2, 2]
delta_t_backward_vec = [2, 2, 2, 2, 2, 2, 2.0, 0.2, 0, 0, 0, 0, 0, 2, 2, 2, 2.5, 2, 0, 0, 0, 0]
# #---------------------------------------------------------------------------------
# Allocate memory for our lovely symmetry
Gamma_u_t_vec = []
Gamma_u_u_vec = []
Gamma_u_v_vec = []
# Loop over the indices and plot the symmetry transformation
for useful_index,new_magical_index in enumerate(magical_indices):
    # Extract the delta_t backwards and forwards
    delta_t_forward = delta_t_forward_vec[useful_index]
    delta_t_backward = delta_t_backward_vec[useful_index]    
    # Take a point on the original solution curve which we transform initially
    X0 = array([t[new_magical_index], u[new_magical_index], v[new_magical_index]])    
    # Calculate our lovely branch index
    if v[new_magical_index]>1:
        branch_index = -1
    else:
        branch_index = 0
    # Try to solve the ODE at hand
    Gamma_epsilon, infodict = integrate.odeint(dX_deps_LV_u, X0, epsilon_vec_dense, args = (alpha,branch_index),full_output=True)
    # Split the solution into its component parts
    Gamma_u_t_temp, Gamma_u_u_temp, Gamma_u_v_temp = Gamma_epsilon.T
    # Find the desired time
    v_val = v[new_magical_index]
    # Tak a lower time value
    t_lower = where(t_2_u==find_nearest(t_2_u, value=t[new_magical_index]-delta_t_backward))[0][0]
    t_upper = where(t_2_u==find_nearest(t_2_u, value=t[new_magical_index]+delta_t_forward))[0][0]
    t_temp = t_2_u[t_lower:t_upper]
    v_temp = v_2_u[t_lower:t_upper]    
    # Find the nearest v_val
    nearest_v_val = find_nearest(v_temp, value=v_val)    
    special_index = where(v_temp==nearest_v_val)[0][0]    
    # Find the desired time point on the transformed curve
    if special_index != len(v_temp)-1:
        t_desired = interp(array([v_val]), array([v_temp[special_index-1], v_temp[special_index], v_temp[special_index+1]]), array([t_temp[special_index-1], t_temp[special_index], t_temp[special_index+1]]))        
    else:
        t_desired = interp(array([v_val]), array([v_temp[special_index-1], v_temp[special_index]]), array([t_temp[special_index-1], t_temp[special_index]]))        
    # Calculate a time constant
    time_constant = ((Gamma_u_t_temp[-1] - t_desired) / (epsilon))
    # Update Gamma_u_temp
    Gamma_u_t_temp = reshape(array([val - time_constant*epsilon_vec_dense[index] for index,val in enumerate(Gamma_u_t_temp)]),Gamma_u_t_temp.size)    
    # Reduce the sparsity of these vectors which are very dense
    Gamma_u_t_sparse = interp(epsilon_vec, epsilon_vec_dense, Gamma_u_t_temp)
    Gamma_u_u_sparse = interp(epsilon_vec, epsilon_vec_dense, Gamma_u_u_temp)
    Gamma_u_v_sparse = interp(epsilon_vec, epsilon_vec_dense, Gamma_u_v_temp)                
    # Save the solution for the current point
    Gamma_u_t_vec.append(Gamma_u_t_sparse)
    Gamma_u_u_vec.append(Gamma_u_u_sparse)
    Gamma_u_v_vec.append(Gamma_u_v_sparse)
# Concatenate to get the full solution curves and the full time
t_2_u_sparse = linspace(t_2_u[0],t_2_u[-1],200)
u_2_u = interp(t_2_u_sparse, t_2_u, u_2_u)
v_2_u = interp(t_2_u_sparse, t_2_u, v_2_u)
t_2_u = t_2_u_sparse
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# v-directional symmetry of the LV model
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# Choose a magical index
magical_index = 100
# Take a point
X0 = array([t[magical_index], u[magical_index], v[magical_index]])
# Calculate the max and min values
vmin = -lambertw(-exp(alpha-H),0).real
vmax = -lambertw(-exp(alpha-H),-1).real
# Loop over the dense epsilon vector and iteratively solve our ODE system for the symmetry
# Calculate our lovely branch index
if u[magical_index]>1:
    branch_index = -1
else:
    branch_index = 0
# Try to solve the ODE at hand
Gamma_epsilon, infodict = integrate.odeint(dX_deps_LV_v, X0, epsilon_vec_dense, args = (alpha,branch_index),full_output=True)
# Split the solution into its component parts
Gamma_v_t, Gamma_v_u, Gamma_v_v = Gamma_epsilon.T
# Define a new time vector 
t_2 = linspace(Gamma_v_t[-1], 10, 20000)
# We need to integrate backwards in time as well to start at 0
t_2_start = linspace(Gamma_v_t[-1], 0, 20000)
# Define new initial conditions for the transformed solutions
X02 = array([Gamma_v_t[-1], Gamma_v_u[-1], Gamma_v_v[-1]])  
X0_2 = array([Gamma_v_u[-1], Gamma_v_v[-1]])  
# Solve the ODE at hand with the new initial conditions
# Integrate backwards to get to 0
X2_start, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2_start, args = (alpha,),full_output=True)
# Find the rest of the solution
X2, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2, args = (alpha,),full_output=True)
#infodict['message'] # >>> 'Integration successful.'
# Split the solution into its component parts
u_2, v_2 = X2.T
# Split the start bit as well
u_2_start, v_2_start = X2_start.T
# Concatenate to get the full solution curves and the full time
u_2_v = concatenate((flip(u_2_start,0), u_2), axis=0)
v_2_v = concatenate((flip(v_2_start,0), v_2), axis=0)
t_2_v = concatenate((flip(t_2_start,0), t_2), axis=0)
#---------------------------------------------------------------------------------
# Plot the symmetry again for some other point on the solution curves
magical_indices_v = [100, 115, 128, 213, 220, 240, 275, 300, 325, 645, 700, 750, 800, 850, 900,-5, -1]
# Illustrative example 1
delta_t_forward_vec = [2, 0.4, 0.4, 0.0, 0.0, 0.0, 0.2, 0.4, 0.4, 0, 0, 0, 0, 0, 0, 0, 0] 
delta_t_backward_vec = [0, 0, 0, 0.4, 0.4, 0.2, 0.2, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1] 
#---------------------------------------------------------------------------------
# Allocate memory for our lovely symmetry
Gamma_v_t_vec = []
Gamma_v_u_vec = []
Gamma_v_v_vec = []
# Loop over the indices and plot the symmetry transformation
for useful_index,new_magical_index in enumerate(magical_indices_v):
    # Extract the delta_t backwards and forwards
    delta_t_forward = delta_t_forward_vec[useful_index]
    delta_t_backward = delta_t_backward_vec[useful_index]    
    # Take a point on the original solution curve which we transform initially
    X0 = array([t[new_magical_index], u[new_magical_index], v[new_magical_index]])    
    # Calculate our lovely branch index
    if u[new_magical_index]>1:
        branch_index = -1
    else:
        branch_index = 0
    # Try to solve the ODE at hand
    Gamma_epsilon, infodict = integrate.odeint(dX_deps_LV_v, X0, epsilon_vec_dense, args = (alpha,branch_index),full_output=True)
    # Split the solution into its component parts
    Gamma_v_t_temp, Gamma_v_u_temp, Gamma_v_v_temp = Gamma_epsilon.T
    # Find the desired time
    u_val = u[new_magical_index]
    # Tak a lower time value
    t_lower = where(t_2_v==find_nearest(t_2_v, value=t[new_magical_index]-delta_t_backward))[0][0]
    t_upper = where(t_2_v==find_nearest(t_2_v, value=t[new_magical_index]+delta_t_forward))[0][0]
    t_temp = t_2_v[t_lower:t_upper]
    u_temp = u_2_v[t_lower:t_upper]    
    # Find the nearest u_val
    nearest_u_val = find_nearest(u_temp, value=u_val)
    special_index = where(u_temp==nearest_u_val)[0][0]    
    if special_index != len(u_temp)-1:
        t_desired = interp(array([u_val]), array([u_temp[special_index-1], u_temp[special_index], u_temp[special_index+1]]), array([t_temp[special_index-1], t_temp[special_index], t_temp[special_index+1]])) 
    else:
        t_desired = interp(array([u_val]), array([u_temp[special_index-1], u_temp[special_index]]), array([t_temp[special_index-1], t_temp[special_index]]))
    # Calculate a time constant
    time_constant = ((Gamma_v_t_temp[-1] - t_desired) / (epsilon))
    # Update Gamma_u_temp
    Gamma_v_t_temp = reshape(array([val - time_constant*epsilon_vec_dense[index] for index,val in enumerate(Gamma_v_t_temp)]),Gamma_v_t_temp.size)
    # Reduce the sparsity of these vectors which are very dense
    Gamma_v_t_sparse = interp(epsilon_vec, epsilon_vec_dense, Gamma_v_t_temp)
    Gamma_v_u_sparse = interp(epsilon_vec, epsilon_vec_dense, Gamma_v_u_temp)
    Gamma_v_v_sparse = interp(epsilon_vec, epsilon_vec_dense, Gamma_v_v_temp)                
    # Save the solution for the current point
    Gamma_v_t_vec.append(Gamma_v_t_sparse)
    Gamma_v_u_vec.append(Gamma_v_u_sparse)
    Gamma_v_v_vec.append(Gamma_v_v_sparse)        
 # Concatenate to get the full solution curves and the full time
t_2_v_sparse = linspace(t_2_v[0],t_2_v[-1],200)
u_2_v = interp(t_2_v_sparse, t_2_v, u_2_v)
v_2_v = interp(t_2_v_sparse, t_2_v, v_2_v)
t_2_v = t_2_v_sparse
#=====================================================================================================
# Se if we can plot our symmetry?
# Define the first figure
f2, ax_2 = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
#-----------------------------------------------------------------------------------------------------
# Subplot 1 out of 2: u-directional symmetry of the LV model
# Prey
ax_2[0].plot(t, u, '-', label="Prey, $u(t)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
# Prey transformed
ax_2[0].plot(t_2_u, u_2_u, '-', label="Prey, $\\hat{u}(t)$" ,color=(35/256,139/256,69/256),linewidth=3.0)
# Predator
ax_2[0].plot(t, v  , '-', label='Predator, $v(t)$',color=(77/256,0/256,75/256),linewidth=3.0)
# Predator transformed
ax_2[0].plot(t_2_u, v_2_u, '-', label="Predator, $\\hat{v}(t)$" ,color=(136/256,65/256,157/256),linewidth=3.0)
# The symmetry without legend
for index in range(len(Gamma_u_t_vec)):
    if index == 0:
        ax_2[0].plot(Gamma_u_t_vec[index],Gamma_u_u_vec[index], '-', label="Symmetry, $\\Gamma_{3,\\epsilon}^{\\mathrm{LV},u}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
        ax_2[0].plot(Gamma_u_t_vec[index],Gamma_u_v_vec[index], '-',color=(0/256,0/256,0/256),linewidth=3.0)
    else:
        ax_2[0].plot(Gamma_u_t_vec[index],Gamma_u_u_vec[index], '-',color=(0/256,0/256,0/256),linewidth=3.0)
        ax_2[0].plot(Gamma_u_t_vec[index],Gamma_u_v_vec[index], '-',color=(0/256,0/256,0/256),linewidth=3.0)
# Grid, legends and other stuff
ax_2[0].grid()
ax_2[0].legend(loc='best',prop={"size":20})
ax_2[0].set_xlabel(xlabel='Time, $t$',fontsize=25)
ax_2[0].set_ylabel(ylabel='Population density',fontsize=25)
ax_2[0].set_yticks([0,1,2,3,4,5]) 
# Change the size of the ticks
ax_2[0].tick_params(axis='both', which='major', labelsize=20)
ax_2[0].tick_params(axis='both', which='minor', labelsize=20)
#-----------------------------------------------------------------------------------------------------
# Subplot 1 out of 2: u-directional symmetry of the LV model
# Prey
ax_2[1].plot(t, u, '-', label="Prey, $u(t)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
# Prey transformed
ax_2[1].plot(t_2_v, u_2_v, '-', label="Prey, $\\hat{u}(t)$" ,color=(35/256,139/256,69/256),linewidth=3.0)
# Predator
ax_2[1].plot(t, v  , '-', label='Predator, $v(t)$',color=(77/256,0/256,75/256),linewidth=3.0)
# Predator transformed
ax_2[1].plot(t_2_v, v_2_v, '-', label="Predator, $\\hat{v}(t)$" ,color=(136/256,65/256,157/256),linewidth=3.0)
# The symmetry without legend
for index in range(len(Gamma_v_t_vec)):
    if index == 0:
        ax_2[1].plot(Gamma_v_t_vec[index],Gamma_v_u_vec[index], '-', label="Symmetry, $\\Gamma_{3,\\epsilon}^{\\mathrm{LV},u}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
        ax_2[1].plot(Gamma_v_t_vec[index],Gamma_v_v_vec[index], '-',color=(0/256,0/256,0/256),linewidth=3.0)
    else:
        ax_2[1].plot(Gamma_v_t_vec[index],Gamma_v_u_vec[index], '-',color=(0/256,0/256,0/256),linewidth=3.0)
        ax_2[1].plot(Gamma_v_t_vec[index],Gamma_v_v_vec[index], '-',color=(0/256,0/256,0/256),linewidth=3.0)
# Grid, legends and other stuff
ax_2[1].grid()
ax_2[1].legend(loc='best',prop={"size":20})
ax_2[1].set_xlabel(xlabel='Time, $t$',fontsize=25)
ax_2[1].set_ylabel(ylabel='Population density',fontsize=25)
ax_2[1].set_yticks([0,1,2,3,4,5]) 
# Change the size of the ticks
ax_2[1].tick_params(axis='both', which='major', labelsize=20)
ax_2[1].tick_params(axis='both', which='minor', labelsize=20)
plt.show()
#-----------------------------------------------------------------------------------------------------
# Title and saving the figure
#-----------------------------------------------------------------------------------------------------
# Title and saving the figure
f2.suptitle('LV symmetries in the time domain',fontsize=30,weight='bold');
f2.savefig('../Figures/time_domain_symmetries_LV.png')
plt.show()
#=================================================================================
#=================================================================================
# Plot time domain symmetries in LaTeX as well...
#=================================================================================
#=================================================================================
#--------------------------------------------------------------------------------
# u-directional
#--------------------------------------------------------------------------------
for index in range(len(Gamma_u_t_vec)):
    if index == 0:
        plot_LaTeX_2D(Gamma_u_t_vec[index],Gamma_u_u_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\\mathrm{LV},u}_{3,\\epsilon}$")
        plot_LaTeX_2D(Gamma_u_t_vec[index],Gamma_u_v_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
    else:
        plot_LaTeX_2D(Gamma_u_t_vec[index],Gamma_u_u_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
        plot_LaTeX_2D(Gamma_u_t_vec[index],Gamma_u_v_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])                
# Plot the solutions as well
plot_LaTeX_2D(t,u,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=r_1,line width=1.5pt,","$u(t)$")
plot_LaTeX_2D(t_2_u,u_2_u,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=r_2,line width=1.5pt,","$\\hat{u}(t)$")
plot_LaTeX_2D(t,v,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=r_3,line width=1.5pt,","$v(t)$")
plot_LaTeX_2D(t_2_u,v_2_u,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=r_4,line width=1.5pt,","$\\hat{v}(t)$")
#--------------------------------------------------------------------------------
# v-directional
#--------------------------------------------------------------------------------
for index in range(len(Gamma_v_t_vec)):
    if index == 0:
        plot_LaTeX_2D(Gamma_v_t_vec[index],Gamma_v_u_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\\mathrm{LV},v}_{3,\\epsilon}$")
        plot_LaTeX_2D(Gamma_v_t_vec[index],Gamma_v_v_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
    else:
        plot_LaTeX_2D(Gamma_v_t_vec[index],Gamma_v_v_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
        plot_LaTeX_2D(Gamma_v_t_vec[index],Gamma_v_u_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])            
# Plot the solutions as well
plot_LaTeX_2D(t,u,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_time.tex","color=r_1,line width=1.5pt,","$u(t)$")
plot_LaTeX_2D(t_2_v,u_2_v,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_time.tex","color=r_2,line width=1.5pt,","$\\hat{u}(t)$")
plot_LaTeX_2D(t,v,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_time.tex","color=r_3,line width=1.5pt,","$v(t)$")
plot_LaTeX_2D(t_2_v,v_2_v,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_time.tex","color=r_4,line width=1.5pt,","$\\hat{v}(t)$")
                      
