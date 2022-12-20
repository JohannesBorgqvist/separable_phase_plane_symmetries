#=================================================================================
#=================================================================================
# Script:"LV"
# Date: 2022-12-20
# Implemented by: Johannes Borgqvist
# Description:
# The script plots the phase plane symmetries as well as the corresponding symmetries of the Lotka-Volterra model in the time domain.
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
# Solving the ODEs for the LV model
#=================================================================================
#=================================================================================
t = linspace(0, 10, 500)              # time
X0 = array([1, 0.10])                     # initials conditions: 10 rabbits and 5 foxes
X1, infodict = integrate.odeint(dX_dt_LV, X0, t, args = (1,),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
#!python
u, v = X1.T
epsilon=0.5
alpha = 1
# Testing to find ICs
# Calculate internal energy
H= v[0] + alpha*u[0] - log((u[0]**alpha)*v[0])
# Define the u-coordinate in the new initial condition
u0 = X0[0]
v0 = X0[0]
# Find the other v-coordinate in the transformed initial condition
v0_new = transform_IC_LV(u0, v0, alpha, H-alpha*epsilon)
v0_new_2 = transform_IC_LV(u0, v0, alpha, H+alpha*epsilon)
# Print our new initial conditions
print("The original IC:\t(u0,v0)\t=\t(%0.3f,%0.3f)"%(u0,v0))
print("The original energy:\tH\t=\t%0.3f"%(H))
print("The transformed IC:\t(u0,v0)\t=\t(%0.3f,%0.3f)"%(u0,v0_new))
print("The transformed energy:\tH\t=\t%0.3f"%(H-alpha*epsilon))
# Try to solve the LV-model again with these new ICs
X2, infodict = integrate.odeint(dX_dt_LV, array([u0,v0_new]), t, args = (1,),full_output=True)
u_transformed, v_transformed = X2.T
X3, infodict = integrate.odeint(dX_dt_LV, array([u0,v0_new_2]), t, args = (1,),full_output=True)
u_transformed_2, v_transformed_2 = X3.T
# Transformed solutions
t_trans = asarray([t_temp + epsilon for t_temp in list(t)])
# Transformations time translation
t_sym = []
epsilon_vec = arange(0,epsilon-0.005,epsilon/50)
t_indices = arange(60,78,3)
for t_index in list(t_indices):
    trans_vec = [t[t_index]+epsilon_temp for epsilon_temp in list(epsilon_vec)]
    t_sym.append(trans_vec)
# Transformations u 
u_sym = []
u_indices = arange(70,125,3)
for u_index in list(u_indices):
    trans_vec = [u_transf(u[u_index],epsilon_temp,alpha) for epsilon_temp in list(epsilon_vec)]
    u_sym.append(trans_vec)
# Transformations v
v_sym = []
v_indices = arange(107,150,2)
for v_index in list(v_indices):
    trans_vec = [v_transf(v[v_index],epsilon_temp,alpha) for epsilon_temp in list(epsilon_vec)]
    v_sym.append(trans_vec)
#=================================================================================
#=================================================================================
# Plotting the solutions
#=================================================================================
#=================================================================================
#Define the first figure
f1, ax_1 = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
# u-directional symmetry
ax_1[0].plot(u, v,'-', label='Original population, $(u,v)$',color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[0].plot(u_transformed_2,v_transformed_2, '-', label='Transformed population, $(\\hat{u},v)$',color=(77/256,0/256,75/256),linewidth=3.0)
ax_1[0].plot(asarray(u_sym[0]),asarray([v[u_indices[0]]*index/index for index in range(len(epsilon_vec))]), '--', label="$\\left.\\Gamma_{\\epsilon}^{u}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,u_index in enumerate(list(u_indices)):
    ax_1[0].plot(asarray(u_sym[index]),asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_1[0].grid()
ax_1[0].legend(loc='best',prop={"size":20})
# v-directional symmetry
ax_1[1].plot(u, v,'-', label='Original population, $(u,v)$',color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[1].plot(u_transformed,v_transformed, '-', label='Transformed population, $(u,\\hat{v})$',color=(77/256,0/256,75/256),linewidth=3.0)
ax_1[1].plot(asarray([u[v_indices[0]]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[0]), '--', label="$\\left.\\Gamma_{\\epsilon}^{v}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,v_index in enumerate(list(v_indices)):
    ax_1[1].plot(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),asarray(v_sym[index]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_1[1].grid()
ax_1[1].legend(loc='best',prop={"size":20})
# Set fontsize of labels
ax_1[0].set_xlabel(xlabel='Rabbits, $u(t)$',fontsize=25)
ax_1[0].set_ylabel(ylabel='Foxes, $v(t)$',fontsize=25)
ax_1[1].set_xlabel(xlabel='Rabbits, $u(t)$',fontsize=25)
ax_1[1].set_ylabel(ylabel='Foxes, $v(t)$',fontsize=25)
# Change the size of the ticks
ax_1[0].tick_params(axis='both', which='major', labelsize=20)
ax_1[0].tick_params(axis='both', which='minor', labelsize=20)
ax_1[1].tick_params(axis='both', which='major', labelsize=20)
ax_1[1].tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
f1.suptitle('Phase plane symmetries of the Lotka-Volterra model',fontsize=30,weight='bold')
f1.savefig('../Figures/phase_plane_symmetries_LV.png')
plt.show()
#=================================================================================
#=================================================================================
# Plot phase plane symmetries in LaTeX as well...
#=================================================================================
#=================================================================================
#--------------------------------------------------------------------------------
# u-directional
#--------------------------------------------------------------------------------
# Plot the solutions in the phase plane
plot_LaTeX_2D(u, v,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex","color=phase_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_transformed_2,v_transformed_2,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex","color=phase_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")
# Plot the symmetries in the phase plane
for index,u_index in enumerate(list(u_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray(u_sym[index]),asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\\mathrm{LV},u}_{2,\\epsilon}$")                
    else:
        plot_LaTeX_2D(asarray(u_sym[index]),asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/LV_symmetries/Input/LV_u_phase.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])                
#--------------------------------------------------------------------------------
# v-directional
#--------------------------------------------------------------------------------
# Plot the solutions in the phase plane
plot_LaTeX_2D(u, v,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex","color=phase_1,line width=1.5pt,","$(u,v)$")
plot_LaTeX_2D(u_transformed,v_transformed,"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex","color=phase_2,line width=1.5pt,","$(\\hat{u},\\hat{v})$")
# Plot the symmetries in the phase plane
for index,v_index in enumerate(list(v_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),asarray(v_sym[index]),"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\\mathrm{LV},v}_{2,\\epsilon}$")                
    else:
        plot_LaTeX_2D(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),asarray(v_sym[index]),"../Figures/LaTeX_figures/LV_symmetries/Input/LV_v_phase.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
