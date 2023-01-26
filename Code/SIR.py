# =================================================================================
# =================================================================================
# Script:"SIR"
# Date: 2022-12-21
# Implemented by: Johannes Borgqvist
# Description:
# The script plots the phase plane symmetries as well as the corresponding symmetries of the SIR model in the time domain.
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
# Import the lambert W function from scipy
from scipy.special import lambertw
# Set all parameters to tex
plt.rcParams['text.usetex'] = True
# =================================================================================
# =================================================================================
# The phase plane symmetries of the SIR model
# =================================================================================
# =================================================================================
#=================================================================================
#=================================================================================
# Plotting the SIR model
#=================================================================================
#=================================================================================
# DEFINE PARAMETER VALUES
r = 1 # In units per days
p = 1 # In units days (defined as p=a/r)
a = r*p # In units per individuals per days
N = 5 # Total population size in units individuals
S0 = N-0.1 # Initial condition for S in units individuals
I0 = N-S0 # Initial condition for I in units individuals
delta = N/50
# CALCULATE PROPERTIES OF THE SIR MODEL
# Define our constant for the solution trajectory
C = I0+S0-p*log(S0)
# Define I_max
I_max = N-p+p*log(p/S0)
# Define a large epsilon
epsilon = 0.5
# Allocate memory for the S-vector
S = linspace(0,N,500)
# Solve for the I-curve
I = asarray([C-S_temp +p*log(S_temp) for S_temp in S])
# Calculate the I directional symmetry as well
I_trans_I_dir = asarray([(C+epsilon)-S_temp +p*log(S_temp) for S_temp in S])
# Calculate the transformed solution
I_trans_S_dir = asarray([(C-epsilon)-S_temp +p*log(S_temp) for S_temp in S])
#Plot the action of the symmetry
epsilon_vec = linspace(0,epsilon,num=200,endpoint=True)
# Transformations S
S_sym = []
#S_indices = [5, 7, 10, 13, 16, 19, 22, 25, 28,30,235,260,290, 310, 340, 373, 420, 455]
S_indices = [13, 16, 19, 22, 25, 28,30,235,260,290, 310, 340, 373, 420, 455]
for S_index in list(S_indices):
    trans_vec = [S_transf(S[S_index],epsilon_temp,a,r) for epsilon_temp in list(epsilon_vec)]
    S_sym.append(list(trans_vec))
# Transformations I
I_sym = []
I_indices = arange(25,500,50)
for I_index in list(I_indices):
    trans_vec = [I[I_index]+epsilon_temp for epsilon_temp in list(epsilon_vec)]
    I_sym.append(list(trans_vec))
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Plot of symmetries for the SIR model
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Plot the solutions
f1, ax_1 = plt.subplots(1,2,constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# S-directional symmetry
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
ax_1[0].plot(S,I, '-', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[0].plot(S,I_trans_S_dir, '-', label="Transformed solution, $(\hat{S},I)$" ,color=(35/256,139/256,69/256),linewidth=3.0)
for index,S_index in enumerate(list(S_indices)):
    if index == 0:
        ax_1[0].plot(asarray(S_sym[index]),asarray([I[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),label="$\\Gamma^{\\mathrm{SIR},S}_{2,\\epsilon}$",linewidth=2.0)
    else:
        ax_1[0].plot(asarray(S_sym[index]),asarray([I[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_1[0].legend(loc='best',prop={"size":20})
ax_1[0].set_xlabel(xlabel='Susceptibles, $S(t)$',fontsize=25)
ax_1[0].set_ylabel(ylabel='Infected, $I(t)$',fontsize=25)
ax_1[0].tick_params(axis='both', which='major', labelsize=20)
ax_1[0].tick_params(axis='both', which='minor', labelsize=20)
ax_1[0].grid()
ax_1[0].set_xlim([0, S[-1]])
ax_1[0].set_ylim([0, max(I_trans_I_dir)+delta])
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# I-directional symmetry
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
ax_1[1].plot(S,I, '-', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[1].plot(S,I_trans_I_dir, '-', label="Transformed solution, $(S,\hat{I})$" ,color=(102/256,194/256,164/256),linewidth=3.0)
for index,I_index in enumerate(list(I_indices)):
    if index == 0:
        ax_1[1].plot(asarray([S[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]), '--' ,color=(0,0,0),label="$\\Gamma^{\\mathrm{SIR},I}_{2,\\epsilon}$",linewidth=2.0)
    else:
        ax_1[1].plot(asarray([S[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]), '--' ,color=(0,0,0),linewidth=2.0)
ax_1[1].tick_params(axis='both', which='major', labelsize=20)
ax_1[1].tick_params(axis='both', which='minor', labelsize=20)
ax_1[1].grid()
ax_1[1].legend(loc='best',prop={"size":20})
ax_1[1].set_xlabel(xlabel='Susceptibles, $S(t)$',fontsize=25)
ax_1[1].set_ylabel(ylabel='Infected, $I(t)$',fontsize=25)
ax_1[1].set_xlim([0, S[-1]])
ax_1[1].set_ylim([0, max(I_trans_I_dir)+delta])
# Title and saving the figure
f1.suptitle('Phase plane symmetries of the SIR- model',
            fontsize=30, weight='bold')
f1.savefig('../Figures/phase_plane_symmetries_SIR.png')
#plt.show()


#=================================================================================
#=================================================================================
# Plotting the solutions in LaTeX
#=================================================================================
#=================================================================================
# S transformation
#=================================================================================
#=================================================================================
# We add the solutions last
plot_LaTeX_2D(S,I,"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=clr_1,line width=2.0pt,","$(S,I)$")
plot_LaTeX_2D(S,I_trans_S_dir,"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=clr_2,line width=2.0pt,","$(\\hat{S},\\hat{I})$")
# The symmetry
for index,S_index in enumerate(list(S_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray(S_sym[index]),asarray([((I[S_index]*(temp_index+1))/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\mathrm{SIR},S}_{2,\\epsilon}$")
    else:
        plot_LaTeX_2D(asarray(S_sym[index]),asarray([((I[S_index]*(temp_index+1))/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,->,>=latex,densely dashed",[])

#=================================================================================
#=================================================================================
# I transformation
#=================================================================================
#=================================================================================
plot_LaTeX_2D(S,I,"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=clr_1,line width=2.0pt,","$(S,I)$")
plot_LaTeX_2D(S,I_trans_I_dir,"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=clr_2,line width=2.0pt,","$(\\hat{S},\\hat{I})$")
# Plot the transformation
for index,I_index in enumerate(list(I_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray([S[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]),"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\mathrm{SIR},I}_{2,\\epsilon}$")
    else:
        plot_LaTeX_2D(asarray([S[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]),"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=black,->,>=latex,densely dashed",[])

#=================================================================================
#=================================================================================
# Plot the lifted symmetries in the time domain
#=================================================================================
#=================================================================================
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# S-directional symmetry of the SIR model
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# Define a large epsilon
epsilon = 0.32
# Plot the action of the symmetry
epsilon_vec_dense = linspace(0,epsilon,num=20000,endpoint=True)
epsilon_vec = linspace(0,epsilon,num=100,endpoint=True)
# Calculate the initial condition (tryouts)
S0_time = 1
t = linspace(0, 3, 500)              # time
# initials conditions: Number of infectives I0 and number of susceptibles S0.
X0 = array([S0, I0])
X1, infodict = integrate.odeint(dX_dt_SIR, X0, t, args=(a,r), full_output=True)
infodict['message']                     # >>> 'Integration successful.'
#!python
S, I = X1.T
# Magical index (tryouts)
magical_index = 20
# Calculate the H-value
H_SIR = I[magical_index] + S[magical_index] - p*log(S[magical_index]) 
# Take a point
X0 = array([t[magical_index], S[magical_index], I[magical_index]])
# Try to solve the ODE at hand
Gamma_epsilon, infodict = integrate.odeint(dX_deps_SIR_S, X0, epsilon_vec_dense, args = (H_SIR,r,p,S0_time),full_output=True)
# Split the solution into its component parts
Gamma_S_t, Gamma_S_S, Gamma_S_I = Gamma_epsilon.T
# Define a new time vector
t_2 = linspace(Gamma_S_t[-1], t[-1], 250)
# We need to integrate backwards in time as well to start at 0
t_2_start = linspace(Gamma_S_t[-1], 0, 250)
# Define new initial conditions for the transformed solutions
X02 = array([Gamma_S_t[-1], Gamma_S_S[-1], Gamma_S_I[-1]])  
X0_2 = array([Gamma_S_S[-1], Gamma_S_I[-1]])  
# Solve the ODE at hand with the new initial conditions
# Integrate backwards to get to 0
X2_start, infodict = integrate.odeint(dX_dt_SIR, X0_2, t_2_start, args=(a,r), full_output=True)
# Find the rest of the solution
X2, infodict = integrate.odeint(dX_dt_SIR, X0_2, t_2, args=(a,r), full_output=True)
infodict['message'] # >>> 'Integration successful.'
# Split the solution into its component parts
S_2, I_2 = X2.T
# Split the start bit as well
S_2_start, I_2_start = X2_start.T
# Concatenate to get the full solution curves and the full time
S_2_S = concatenate((flip(S_2_start,0), S_2), axis=0)
I_2_S = concatenate((flip(I_2_start,0), I_2), axis=0)
t_2_S = concatenate((flip(t_2_start,0), t_2), axis=0)
# Plot the symmetry again for some other point on the solution curves
magical_indices = [20, 40, 60, 80, 100, 120, 130, 140, 150, 160, 167, 173, 295, 305, 315, 325, 335, 345]
# Allocate memory for our lovely symmetry
Gamma_S_t_vec = []
Gamma_S_S_vec = []
Gamma_S_I_vec = []
# Loop over the indices and plot the symmetry transformation
for index,new_magical_index in enumerate(magical_indices):
    # Take a point
    X0 = array([t[new_magical_index], S[new_magical_index], I[new_magical_index]])  
    # Try to solve the ODE at hand
    Gamma_epsilon_temp, infodict = integrate.odeint(dX_deps_SIR_S, X0, epsilon_vec_dense, args = (H_SIR,r,p,S0_time),full_output=True)    
    # Split the solution into its component parts
    Gamma_S_t_temp, Gamma_S_S_temp, Gamma_S_I_temp = Gamma_epsilon_temp.T
    # Make these vector more sparse than they are currently by means of interpolation
    Gamma_S_t_temp = interp(epsilon_vec, epsilon_vec_dense, Gamma_S_t_temp)
    Gamma_S_S_temp = interp(epsilon_vec, epsilon_vec_dense, Gamma_S_S_temp)
    Gamma_S_I_temp = interp(epsilon_vec, epsilon_vec_dense, Gamma_S_I_temp)
    # Append our solutions
    Gamma_S_t_vec.append(Gamma_S_t_temp)    
    Gamma_S_S_vec.append(Gamma_S_S_temp)    
    Gamma_S_I_vec.append(Gamma_S_I_temp)
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# I-directional symmetry of the SIR model
# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# Define a large epsilon
epsilon = 0.07
#Plot the action of the symmetry
epsilon_vec = linspace(0,epsilon,num=200,endpoint=True)
epsilon_vec_dense = linspace(0,epsilon,num=20000,endpoint=True)
# Magical index
magical_index = 170
branch_index_temp = -1
# Make a test value for the lower limit in the integral
I0_test = 0.1
# Take a point
X0 = array([t[magical_index], S[magical_index], I[magical_index]])
# Try to solve the ODE at hand
Gamma_epsilon, infodict = integrate.odeint(dX_deps_SIR_I, X0, epsilon_vec_dense, args = (H_SIR,r,p,I0_test,branch_index_temp),full_output=True)
# Split the solution into its component parts
Gamma_I_t, Gamma_I_S, Gamma_I_I = Gamma_epsilon.T
# Define a new time vector
t_2 = linspace(Gamma_I_t[-1], t[-1], 250)
# We need to integrate backwards in time as well to start at 0
t_2_start = linspace(Gamma_I_t[-1], 0, 250)
# Define new initial conditions for the transformed solutions
X02 = array([Gamma_I_t[-1], Gamma_I_S[-1], Gamma_I_I[-1]])  
X0_2 = array([Gamma_I_S[-1], Gamma_I_I[-1]])  
# Solve the ODE at hand with the new initial conditions
# Integrate backwards to get to 0
X2_start, infodict = integrate.odeint(dX_dt_SIR, X0_2, t_2_start, args=(a,r), full_output=True)
# Find the rest of the solution
X2, infodict = integrate.odeint(dX_dt_SIR, X0_2, t_2, args=(a,r), full_output=True)
infodict['message'] # >>> 'Integration successful.'
# Split the solution into its component parts
S_2, I_2 = X2.T
# Split the start bit as well
S_2_start, I_2_start = X2_start.T
# Concatenate to get the full solution curves and the full time
S_2_I = concatenate((flip(S_2_start,0), S_2), axis=0)
I_2_I = concatenate((flip(I_2_start,0), I_2), axis=0)
t_2_I = concatenate((flip(t_2_start,0), t_2), axis=0)
# Plot the symmetry again for some other point on the solution curves
#magical_indices_I = [20,25,30,35,65, 75, 85, 90]
#branches_I = [0,0,0,0,0, 0, 0, 0]
magical_indices_I = [120, 150,175,195, 280, 300, 320, 340]
branches_I = [-1, -1,-1, -1,0, 0, 0, 0]
# Allocate memory for our lovely symmetry
Gamma_I_t_vec = []
Gamma_I_S_vec = []
Gamma_I_I_vec = []
# Loop over the indices and plot the symmetry transformation
for index,new_magical_index in enumerate(magical_indices_I):
    # Take a point
    X0 = array([t[new_magical_index], S[new_magical_index], I[new_magical_index]])  
    # Try to solve the ODE at hand
    Gamma_epsilon_temp, infodict = integrate.odeint(dX_deps_SIR_I, X0, epsilon_vec_dense, args = (H_SIR,r,p,I0_test,branches_I[index]),full_output=True)    
    # Split the solution into its component parts
    Gamma_I_t_temp, Gamma_I_S_temp, Gamma_I_I_temp = Gamma_epsilon_temp.T    
    # Make these vector more sparse than they are currently by means of interpolation
    Gamma_I_t_temp = interp(epsilon_vec, epsilon_vec_dense, Gamma_I_t_temp)
    Gamma_I_S_temp = interp(epsilon_vec, epsilon_vec_dense, Gamma_I_S_temp)
    Gamma_I_I_temp = interp(epsilon_vec, epsilon_vec_dense, Gamma_I_I_temp)
    # Append our solutions
    Gamma_I_t_vec.append(Gamma_I_t_temp)    
    Gamma_I_S_vec.append(Gamma_I_S_temp)    
    Gamma_I_I_vec.append(Gamma_I_I_temp)
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Plot of symmetries for the SIR model in the time domain
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Plot the solutions
f2, ax_2 = plt.subplots(1,2,constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# S-directional symmetry
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# Original solution S
ax_2[0].plot(t,S, '-', label="$S$" ,color=(0/256,68/256,27/256),linewidth=3.0)
# Transformed solution S
ax_2[0].plot(t_2_S,S_2_S, '-', label="$\\hat{S}$" ,color=(35/256,139/256,69/256),linewidth=3.0)
# Original solution I
ax_2[0].plot(t,I, '-', label="$I$" ,color=(77/256,0/256,75/256),linewidth=3.0)
# # Transformed solution I
ax_2[0].plot(t_2_S,I_2_S, '-', label="$\\hat{I}$" ,color=(136/256,65/256,157/256),linewidth=3.0)
# The symmetry without legend
for index in range(len(Gamma_S_t_vec)):
    if index == 0:
        ax_2[0].plot(Gamma_S_t_vec[index],Gamma_S_S_vec[index], '--', label="Symmetry, $\\Gamma_{3,\\epsilon}^{\\mathrm{SIR},S}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
        ax_2[0].plot(Gamma_S_t_vec[index],Gamma_S_I_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
    else:
        ax_2[0].plot(Gamma_S_t_vec[index],Gamma_S_S_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
        ax_2[0].plot(Gamma_S_t_vec[index],Gamma_S_I_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
# General axes properties
ax_2[0].legend(loc='best',prop={"size":20})
ax_2[0].set_xlabel(xlabel='Time, $t$ [days]',fontsize=25)
ax_2[0].set_ylabel(ylabel='Population size [\\# individuals]',fontsize=25)
ax_2[0].tick_params(axis='both', which='major', labelsize=20)
ax_2[0].tick_params(axis='both', which='minor', labelsize=20)
ax_2[0].grid()
# #---------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------
# # I-directional symmetry
# #---------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------
# Solution S
ax_2[1].plot(t,S, '-', label="$S$" ,color=(0/256,68/256,27/256),linewidth=3.0)
# # Transformed solution S
ax_2[1].plot(t_2_I,S_2_I, '-', label="$\\hat{S}$" ,color=(35/256,139/256,69/256),linewidth=3.0)
# Solution I
ax_2[1].plot(t,I, '-', label="$I$" ,color=(77/256,0/256,75/256),linewidth=3.0)
# Transformed solution I
ax_2[1].plot(t_2_I,I_2_I, '-', label="$\\hat{I}$" ,color=(136/256,65/256,157/256),linewidth=3.0)
#The symmetry without legend
for index in range(len(Gamma_I_t_vec)):
    if index == 0:
        ax_2[1].plot(Gamma_I_t_vec[index],Gamma_I_S_vec[index], '--', label="Symmetry, $\\Gamma_{3,\\epsilon}^{\\mathrm{SIR},I}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
        ax_2[1].plot(Gamma_I_t_vec[index],Gamma_I_I_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
    else:
        ax_2[1].plot(Gamma_I_t_vec[index],Gamma_I_S_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
        ax_2[1].plot(Gamma_I_t_vec[index],Gamma_I_I_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
#General ticks settings etc.
ax_2[1].tick_params(axis='both', which='major', labelsize=20)
ax_2[1].tick_params(axis='both', which='minor', labelsize=20)
ax_2[1].grid()
ax_2[1].legend(loc='best',prop={"size":20})
ax_2[1].set_xlabel(xlabel='Time, $t$ [days]',fontsize=25)
ax_2[1].set_ylabel(ylabel='Population size [\\# individuals]',fontsize=25)
#-----------------------------------------------------------------------------------------------------
# Title and saving the figure
#-----------------------------------------------------------------------------------------------------
# Title and saving the figure
f2.suptitle('SIR symmetries in the time domain',fontsize=30,weight='bold');
f2.savefig('../Figures/time_domain_symmetries_SIR.png')


I_max = H_SIR - (p-p*log(p))
I_temp = linspace(0.01,I_max)
S_temp = array([-p*lambertw(-exp(((I_val-H_SIR)/(p)))/p,-1) for I_val in I_temp])
S_temp_2 = array([-p*lambertw(-exp(((I_val-H_SIR)/(p)))/p,0) for I_val in I_temp])

f3 = plt.figure(figsize =(5, 4)) 
ax_3 = f3.add_axes([0.1, 0.1, 0.8, 0.8])
ax_3.plot(S_temp, I_temp,label="Branch -1") 
ax_3.plot(S_temp_2, I_temp,label="Branch 0")
ax_3.grid()
ax_3.legend(loc='best',prop={"size":20})
ax_3.set_xlabel(xlabel='Susceptibles, $S(t)$',fontsize=25)
ax_3.set_ylabel(ylabel='Infected, $I(t)$',fontsize=25)
f3.suptitle('matplotlib.figure.Figure() class Example\n\n',  
             fontweight ="bold") 

plt.show()

#=================================================================================
#=================================================================================
# Plot time domain symmetries in LaTeX as well...
#=================================================================================
#=================================================================================
#--------------------------------------------------------------------------------
# S-directional
#--------------------------------------------------------------------------------
# Plot the symmetries
for index in range(len(Gamma_S_t_vec)):
    if index == 0:
        plot_LaTeX_2D(Gamma_S_t_vec[index],Gamma_S_S_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\\mathrm{SIR},S}_{3,\\epsilon}$")
        plot_LaTeX_2D(Gamma_S_t_vec[index],Gamma_S_I_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
    else:
        if index < 6 or index >12:
            plot_LaTeX_2D(Gamma_S_t_vec[index],Gamma_S_S_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
        if index>7:
            plot_LaTeX_2D(Gamma_S_t_vec[index],Gamma_S_I_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])                
# Plot the solutions as well
plot_LaTeX_2D(t,S,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=r_1,line width=1.5pt,","$S(t)$")
plot_LaTeX_2D(t_2_S,S_2_S,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=r_2,line width=1.5pt,","$\\hat{S}(t)$")
plot_LaTeX_2D(t,I,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=r_3,line width=1.5pt,","$I(t)$")
plot_LaTeX_2D(t_2_S,I_2_S,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=r_4,line width=1.5pt,","$\\hat{I}(t)$")


#--------------------------------------------------------------------------------
# I-directional
#--------------------------------------------------------------------------------
for index in range(len(Gamma_I_t_vec)):
    if index == 0:
        plot_LaTeX_2D(Gamma_I_t_vec[index],Gamma_I_S_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\\mathrm{SIR},I}_{3,\\epsilon}$")
        plot_LaTeX_2D(Gamma_I_t_vec[index],Gamma_I_I_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
    else:
        plot_LaTeX_2D(Gamma_I_t_vec[index],Gamma_I_S_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
        plot_LaTeX_2D(Gamma_I_t_vec[index],Gamma_I_I_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])                
# Plot the solutions as well
plot_LaTeX_2D(t,S,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=r_1,line width=1.5pt,","$S(t)$")
plot_LaTeX_2D(t_2_I,S_2_I,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=r_2,line width=1.5pt,","$\\hat{S}(t)$")
plot_LaTeX_2D(t,I,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=r_3,line width=1.5pt,","$I(t)$")
plot_LaTeX_2D(t_2_I,I_2_I,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=r_4,line width=1.5pt,","$\\hat{I}(t)$")
