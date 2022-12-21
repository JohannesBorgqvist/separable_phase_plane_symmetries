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
# Transformation parameter
epsilon = 0.5
# Dimensionless LV parameter
alpha = 1
t = linspace(0, 10, 500)              # time
# initials conditions: 10 rabbits and 5 foxes
X0 = array([1, 0.10])
X1, infodict = integrate.odeint(dX_dt_LV, X0, t, args=(alpha,), full_output=True)
infodict['message']                     # >>> 'Integration successful.'
#!python
u, v = X1.T
# Testing to find ICs
# Calculate internal energy
H = v[0] + alpha*u[0] - log((u[0]**alpha)*v[0])
# Define the u-coordinate in the new initial condition
u0 = X0[0]
v0 = X0[0]
# Find the other v-coordinate in the transformed initial condition
v0_new = transform_IC_LV(u0, v0, alpha, H-alpha*epsilon)
v0_new_2 = transform_IC_LV(u0, v0, alpha, H+alpha*epsilon)
# Print our new initial conditions
print("The original IC:\t(u0,v0)\t=\t(%0.3f,%0.3f)" % (u0, v0))
print("The original energy:\tH\t=\t%0.3f" % (H))
print("The transformed IC:\t(u0,v0)\t=\t(%0.3f,%0.3f)" % (u0, v0_new))
print("The transformed energy:\tH\t=\t%0.3f" % (H-alpha*epsilon))
# Try to solve the LV-model again with these new ICs
X2, infodict = integrate.odeint(dX_dt_LV, array(
    [u0, v0_new]), t, args=(alpha,), full_output=True)
u_transformed, v_transformed = X2.T
X3, infodict = integrate.odeint(dX_dt_LV, array(
    [u0, v0_new_2]), t, args=(1,), full_output=True)
u_transformed_2, v_transformed_2 = X3.T
# Transformed solutions
t_trans = asarray([t_temp + epsilon for t_temp in list(t)])
# Transformations time translation
t_sym = []
epsilon_vec = arange(0, epsilon-0.005, epsilon/50)
t_indices = arange(60, 78, 3)
for t_index in list(t_indices):
    trans_vec = [t[t_index]+epsilon_temp for epsilon_temp in list(epsilon_vec)]
    t_sym.append(trans_vec)
# Transformations u
u_sym = []
u_indices = arange(70, 125, 3)
for u_index in list(u_indices):
    trans_vec = [u_transf(u[u_index], epsilon_temp, alpha)
                 for epsilon_temp in list(epsilon_vec)]
    u_sym.append(trans_vec)
# Transformations v
v_sym = []
v_indices = arange(107, 150, 2)
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
             label='Transformed population, $(\\hat{u},v)$', color=(77/256, 0/256, 75/256), linewidth=3.0)
ax_1[0].plot(asarray(u_sym[0]), asarray([v[u_indices[0]]*index/index for index in range(len(epsilon_vec))]), '--',
             label="$\\left.\\Gamma_{\\epsilon}^{u}\\right|_{\\epsilon=" + str(epsilon) + "}$", color=(0, 0, 0), linewidth=2.0)
for index, u_index in enumerate(list(u_indices)):
    ax_1[0].plot(asarray(u_sym[index]), asarray([v[u_index]*((index+1)/(index+1))
                                                 for index in range(len(epsilon_vec))]), '--', color=(0, 0, 0), linewidth=2.0)
ax_1[0].grid()
ax_1[0].legend(loc='best', prop={"size": 20})
# v-directional symmetry
ax_1[1].plot(u, v, '-', label='Original population, $(u,v)$',
             color=(0/256, 68/256, 27/256), linewidth=3.0)
ax_1[1].plot(u_transformed, v_transformed, '-',
             label='Transformed population, $(u,\\hat{v})$', color=(77/256, 0/256, 75/256), linewidth=3.0)
ax_1[1].plot(asarray([u[v_indices[0]]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[0]),
             '--', label="$\\left.\\Gamma_{\\epsilon}^{v}\\right|_{\\epsilon=" + str(epsilon) + "}$", color=(0, 0, 0), linewidth=2.0)
for index, v_index in enumerate(list(v_indices)):
    ax_1[1].plot(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),
                 asarray(v_sym[index]), '--', color=(0, 0, 0), linewidth=2.0)
ax_1[1].grid()
ax_1[1].legend(loc='best', prop={"size": 20})
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
# Plot the solutions in the phase plane
plot_LaTeX_2D(u, v, "../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex","color=phase_1,line width=1.5pt,", "$(u,v)$")
plot_LaTeX_2D(u_transformed_2, v_transformed_2, "../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex","color=phase_2,line width=1.5pt,", "$(\\hat{u},\\hat{v})$")
# Plot the symmetries in the phase plane
for index, u_index in enumerate(list(u_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray(u_sym[index]), asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex", "color=black,->,>=latex,densely dashed,line width=1.0pt", "$\\Gamma^{\\mathrm{LV},u}_{2,\\epsilon}$")
    else:
        plot_LaTeX_2D(asarray(u_sym[index]), asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex", "color=black,->,>=latex,densely dashed,line width=1.0pt", [])
# --------------------------------------------------------------------------------
# v-directional
# --------------------------------------------------------------------------------
# Plot the solutions in the phase plane
plot_LaTeX_2D(u, v, "../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex","color=phase_1,line width=1.5pt,", "$(u,v)$")
plot_LaTeX_2D(u_transformed, v_transformed, "../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex","color=phase_2,line width=1.5pt,", "$(\\hat{u},\\hat{v})$")
# Plot the symmetries in the phase plane
for index, v_index in enumerate(list(v_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[index]),
                      "../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex", "color=black,->,>=latex,densely dashed,line width=1.0pt", "$\\Gamma^{\\mathrm{LV},v}_{2,\\epsilon}$")
    else:
        plot_LaTeX_2D(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[index]),
                      "../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex", "color=black,->,>=latex,densely dashed,line width=1.0pt", [])
# =================================================================================
# =================================================================================
# Plot the lifted symmetries as well
# =================================================================================
# =================================================================================
# Epsilon value
epsilon = 1.0
# The transformation parameter
epsilon_vec = linspace(0,epsilon,200)              # epsilon
# We know that Lambertz w satisfies the following (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lambertw.html):
# The principal branch (k=0) is real if z>-1/e,
# The branch (k>1) is real if -1/e < z < 0.
# Let's make an experiment.
# Magical index
magical_index = 25
#magical_index = 89
# Calculate the minimal u value
u_min = min(u)
u_max = max(u)
# Set a lower integration limit
u0 = 1.0
# Prompt to the user
print("\n\n")
print("u_min\t=\t%0.3f"%(u_min))
print("u_max\t=\t%0.3f"%(u_max))
print("u_0\t=\t%0.3f"%(u0))
print("\n\n")
# Take a point
X0 = array([t[magical_index], u[magical_index], v[magical_index]])
# Try to solve the ODE at hand
#Gamma_epsilon, infodict = integrate.odeint(dX_deps_LV_u, X0, epsilon_vec, args = (alpha,-1,u0),full_output=True)
Gamma_epsilon, infodict = integrate.odeint(dX_deps_LV_u, X0, epsilon_vec, args = (alpha,0,u0),full_output=True)
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
X2_start, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2_start, args = (alpha,),full_output=True)
# Find the rest of the solution
X2, infodict = integrate.odeint(dX_dt_LV, X0_2, t_2, args = (alpha,),full_output=True)
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
magical_indices = [25, 35, 45, 55, 65, 75, 85, 100, 105, 110, 115, 120, 130]
#branch_indices = [-1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0]
branch_indices = [0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1]
# Allocate memory for our lovely symmetry
Gamma_u_t_vec = []
Gamma_u_u_vec = []
Gamma_u_v_vec = []
# Loop over the indices and plot the symmetry transformation
for branch_index,new_magical_index in enumerate(magical_indices):
    # Take a point
    X0 = array([t[new_magical_index], u[new_magical_index], v[new_magical_index]])  
    # Try to solve the ODE at hand
    Gamma_epsilon_temp, infodict = integrate.odeint(dX_deps_LV_u, X0, epsilon_vec, args = (alpha,branch_indices[branch_index],u0),full_output=True)    
    # Split the solution into its component parts
    Gamma_u_t_temp, Gamma_u_u_temp, Gamma_u_v_temp = Gamma_epsilon_temp.T    
    # Append our solutions
    Gamma_u_t_vec.append(Gamma_u_t_temp)    
    Gamma_u_u_vec.append(Gamma_u_u_temp)    
    Gamma_u_v_vec.append(Gamma_u_v_temp)    
#=====================================================================================================
# Se if we can plot our symmetry?
fig_2 = plt.figure(constrained_layout=True, figsize=(20, 8))
# Prey
plt.plot(t, u, '-', label="Prey, $u(t)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
# Prey transformed
plt.plot(t_2, u_2, '-', label="Prey, $\\hat{u}(t)$" ,color=(35/256,139/256,69/256),linewidth=3.0)
# Predator
plt.plot(t, v  , '-', label='Predator, $v(t)$',color=(77/256,0/256,75/256),linewidth=3.0)
# Prey transformed
plt.plot(t_2, v_2, '-', label="Predator, $\\hat{v}(t)$" ,color=(136/256,65/256,157/256),linewidth=3.0)
# The symmetry without legend
for index in range(len(Gamma_u_t_vec)):
    if index == 0:
        plt.plot(Gamma_u_t_vec[index],Gamma_u_u_vec[index], '--', label="Symmetry, $\\Gamma_{3,\\epsilon}^{\\mathrm{LV},u}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
        plt.plot(Gamma_u_t_vec[index],Gamma_u_v_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
    else:
        plt.plot(Gamma_u_t_vec[index],Gamma_u_u_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
        plt.plot(Gamma_u_t_vec[index],Gamma_u_v_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
plt.grid()
plt.legend(loc='best',prop={"size":20})
plt.xlabel(xlabel='Time, $t$',fontsize=25)
plt.ylabel(ylabel='Population density',fontsize=25)
# Change the size of the ticks
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
plt.title('LV symmetries in the time domain',fontsize=30,weight='bold');
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
        if magical_indices[index]<100 or magical_indices[index]>115:
            plot_LaTeX_2D(Gamma_u_t_vec[index],Gamma_u_v_vec[index],"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])                
# Plot the solutions as well
plot_LaTeX_2D(t,u,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=r_1,line width=1.5pt,","$u(t)$")
plot_LaTeX_2D(t_2,u_2,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=r_2,line width=1.5pt,","$\\hat{u}(t)$")
plot_LaTeX_2D(t,v,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=r_3,line width=1.5pt,","$v(t)$")
plot_LaTeX_2D(t_2,v_2,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_time.tex","color=r_4,line width=1.5pt,","$\\hat{v}(t)$")

                      
