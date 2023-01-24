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
r = 0.00218 # In units per days
p = 202 # In units days (defined as p=a/r)
a = r*p # In units per individuals per days
#S0_ori = 300 # Initial condition for S in units individuals
N = 763 # Total population size in units individuals
S0 = N-1 # Initial condition for S in units individuals
I0 = N-S0 # Initial condition for I in units individuals
# CALCULATE PROPERTIES OF THE SIR MODEL
# Define our constant for the solution trajectory
C = I0+S0-p*log(S0)
# Define I_max
I_max = N-p+p*log(p/S0)
# Define a large epsilon
epsilon = 100
# Allocate memory for the S-vector
S = linspace(0,S0,500)
I = asarray([C-S_temp +p*log(S_temp) for S_temp in S])
#I_trans = asarray([(C+epsilon)-S_temp +p*log(S_temp) for S_temp in S])
#I_trans_2 = asarray([(C+2*epsilon)-S_temp +p*log(S_temp) for S_temp in S])
I_trans = asarray([(C-epsilon)-S_temp +p*log(S_temp) for S_temp in S])
I_trans_2 = asarray([(C-2*epsilon)-S_temp +p*log(S_temp) for S_temp in S])
# Calculate the I directional symmetry as well
I_trans_I_dir = asarray([(C+epsilon)-S_temp +p*log(S_temp) for S_temp in S])
# Calculate the transformed solution
S_transformed = asarray([S_transf(S_temp, epsilon,a,r) for S_temp in list(S)])
# Modify S and I so that we only plot up until the threshold N
S_new = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I[index])<=N and I[index]>0])
I_new = asarray([I[index] for index,S_temp in enumerate(S) if (S_temp+I[index])<=N and I[index]>0])
S_trans_I_dir_new = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I_trans_I_dir[index])<=N and I_trans_I_dir[index]>0])
I_trans_I_dir_new = asarray([I_trans_I_dir[index] for index,S_temp in enumerate(S) if (S_temp+I_trans_I_dir[index])<=N and I_trans_I_dir[index]>0])
# The same goes for the transformed solution
S_trans_new = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I_trans[index])<=N and I_trans[index]>0])
I_trans_new = asarray([I_trans[index] for index,S_temp in enumerate(S) if (S_temp+I_trans[index])<=N and I_trans[index]>0])
S_trans_new_2 = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I_trans_2[index])<=N and I_trans_2[index]>0])
I_trans_new_2 = asarray([I_trans_2[index] for index,S_temp in enumerate(S) if (S_temp+I_trans_2[index])<=N and I_trans_2[index]>0])
#I_trans = asarray([I[index] for index,S_temp in enumerate(S_new) if (S_temp+I[index])<=N and I[index]>0])
# Plot the action of the symmetry
epsilon_vec = linspace(0,epsilon,num=100,endpoint=True)
# Transformations S
S_sym = []
S_indices = [3, 10, 17, 24, 300, 350, 400, 450]
for S_index in list(S_indices):
    trans_vec = [S_transf(S_new[S_index],epsilon_temp,a,r) for epsilon_temp in list(epsilon_vec)]
    S_sym.append(list(trans_vec))
# Add the other ones as well
S_sym_2 = []
S_indices_2 = arange(21,30,2)
for S_index in list(S_indices_2):
    trans_vec = [S_transf(S_trans_new[S_index],epsilon_temp,a,r) for epsilon_temp in list(epsilon_vec)]
    S_sym_2.append(list(trans_vec))
# Transformations I
I_sym = []
I_indices = arange(25,300,50)
for I_index in list(I_indices):
    trans_vec = [I_new[I_index]+epsilon_temp for epsilon_temp in list(epsilon_vec)]
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
#ax_1.plot(S,I, '--', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[0].plot(S_new,I_new, '-', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[0].plot(S_trans_new,I_trans_new, '-', label="Transformed solution, $(\hat{S},I)$" ,color=(35/256,139/256,69/256),linewidth=3.0)
for index,S_index in enumerate(list(S_indices)):
    if index == 0:
        ax_1[0].plot(asarray(S_sym[index]),asarray([I_new[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),label="$\\Gamma^{\\mathrm{SIR},S}_{\\epsilon}$",linewidth=2.0)
    else:
        ax_1[0].plot(asarray(S_sym[index]),asarray([I_new[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_1[0].plot(asarray([N, 0]),asarray([0, N]),label="Total population size, $N=I+S$",color=(0,0,0),linewidth=3.0)
ax_1[0].plot(asarray([0, p]),asarray([I_max, I_max]), '--' ,color=(0,0,0),linewidth=3.0)
#ax_1[0].plot(asarray([0, p]),asarray([I_max+epsilon, I_max+epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[0].plot(asarray([0, p]),asarray([I_max-epsilon, I_max-epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[0].plot(asarray([p, p]),asarray([0, I_max]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[0].legend(loc='best',prop={"size":20})
ax_1[0].set_xlabel(xlabel='Susceptibles, $S(t)$',fontsize=25)
ax_1[0].set_ylabel(ylabel='Infected, $I(t)$',fontsize=25)
ax_1[0].tick_params(axis='both', which='major', labelsize=20)
ax_1[0].tick_params(axis='both', which='minor', labelsize=20)
ax_1[0].grid()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# I-directional symmetry
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
ax_1[1].plot(S_new,I_new, '-', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[1].plot(S_trans_I_dir_new,I_trans_I_dir_new, '-', label="Transformed solution, $(\hat{\hat{S}},I)$" ,color=(102/256,194/256,164/256),linewidth=3.0)
for index,I_index in enumerate(list(I_indices)):
    if index == 0:
        ax_1[1].plot(asarray([S_new[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]), '--' ,color=(0,0,0),label="$\\Gamma^{\\mathrm{SIR},I}_{\\epsilon}$",linewidth=2.0)
    else:
        ax_1[1].plot(asarray([S_new[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]), '--' ,color=(0,0,0),linewidth=2.0)
ax_1[1].plot(asarray([N, 0]),asarray([0, N]),label="Total population size, $N=I+S$",color=(0,0,0),linewidth=3.0)
ax_1[1].plot(asarray([0, p]),asarray([I_max, I_max]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[1].plot(asarray([0, p]),asarray([I_max+epsilon, I_max+epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[1].plot(asarray([p, p]),asarray([0, I_max]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[1].plot(asarray([p, p]),asarray([0, I_max+epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[1].tick_params(axis='both', which='major', labelsize=20)
ax_1[1].tick_params(axis='both', which='minor', labelsize=20)
ax_1[1].grid()
ax_1[1].legend(loc='best',prop={"size":20})
ax_1[1].set_xlabel(xlabel='Susceptibles, $S(t)$',fontsize=25)
ax_1[1].set_ylabel(ylabel='Infected, $I(t)$',fontsize=25)
# Title and saving the figure
f1.suptitle('Phase plane symmetries of the SIR- model',
            fontsize=30, weight='bold')
f1.savefig('../Figures/phase_plane_symmetries_SIR.png')
plt.show()


# #=================================================================================
# #=================================================================================
# # Plotting the solutions in LaTeX
# #=================================================================================
# #=================================================================================
# # S transformation
# #=================================================================================
# #=================================================================================
# # The symmetry
# for index,S_index in enumerate(list(S_indices)):
#     if index == 0:
#         plot_LaTeX_2D(asarray(S_sym[index]),asarray([((I_new[S_index]*(temp_index+1))/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\mathrm{SIR},S}_{2,\\epsilon}$")
#     elif index<5:
#         plot_LaTeX_2D(asarray(S_sym[index]),asarray([((I_new[S_index]*(temp_index+1))/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,->,>=latex,densely dashed",[])
# # We add the solutions last
# plot_LaTeX_2D(S_new,I_new,"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=clr_1,line width=2.0pt,","$I(S)$")
# plot_LaTeX_2D(S_trans_new,I_trans_new,"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=clr_2,line width=2.0pt,","$\\hat{I}(S;\\epsilon)$")
# # for index in range(len(S_sym_2)):
# #     plot_LaTeX_2D(asarray(S_sym_2[index]),asarray([I_trans_new[S_indices_2[index]]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
# plot_LaTeX_2D([N*x for x in linspace(0,1,50)],[N*(1-x) for x in linspace(0,1,50)],"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,mark=diamond*,only marks,line width=0.75pt,","$N=I+S$")
# plot_LaTeX_2D([0, p],[I_max,I_max],"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
# plot_LaTeX_2D([0, p],[I_max-epsilon,I_max-epsilon],"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
# #plot_LaTeX_2D([0, p],[I_max+2*epsilon,I_max+2*epsilon],"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
# plot_LaTeX_2D([p, p],[0,I_max],"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
# #plot_LaTeX_2D(S_trans_new_2,I_trans_new_2,"../Figures/LaTeX_figures/SIR_symmetries/Input/S_trans.tex","color=clr_3,line width=2.0pt,","$\\hat{\\hat{I}}(S;\\epsilon)=\\hat{I}(S;2\\epsilon)$")
# #=================================================================================
# #=================================================================================
# # I transformation
# #=================================================================================
# #=================================================================================
# # Plot the transformation
# for index,I_index in enumerate(list(I_indices)):
#     if index == 0:
#         plot_LaTeX_2D(asarray([S_new[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]),"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\mathrm{SIR},I}_{2,\\epsilon}$")
#     elif index<5:
#         plot_LaTeX_2D(asarray([S_new[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]),"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=black,->,>=latex,densely dashed",[])
# plot_LaTeX_2D(S_new,I_new,"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=clr_1,line width=2.0pt,","$I(S)$")
# plot_LaTeX_2D(S_trans_I_dir_new,I_trans_I_dir_new,"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=clr_2,line width=2.0pt,","$\\hat{I}(S;\\epsilon)$")
# plot_LaTeX_2D([N*x for x in linspace(0,1,50)],[N*(1-x) for x in linspace(0,1,50)],"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=black,mark=diamond*,only marks,line width=0.75pt,","$N=I+S$")
# plot_LaTeX_2D([0, p],[I_max,I_max],"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
# plot_LaTeX_2D([0, p],[I_max+epsilon,I_max+epsilon],"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
# plot_LaTeX_2D([p, p],[0,I_max+epsilon],"../Figures/LaTeX_figures/SIR_symmetries/Input/I_trans.tex","color=black,densely dashed,line width=1.0pt,",[])        


# #=================================================================================
# #=================================================================================
# # Plot the lifted symmetries in the time domain
# #=================================================================================
# #=================================================================================
# # ---------------------------------------------------------------------------------
# # ---------------------------------------------------------------------------------
# # S-directional symmetry of the SIR model
# # ---------------------------------------------------------------------------------
# # ---------------------------------------------------------------------------------
# t = linspace(0, 15, 500)              # time
# # initials conditions: Number of infectives I0 and number of susceptibles S0.
# X0 = array([S0, I0])
# X1, infodict = integrate.odeint(dX_dt_SIR, X0, t, args=(a,r), full_output=True)
# infodict['message']                     # >>> 'Integration successful.'
# #!python
# S, I = X1.T
# # Magical index
# magical_index = 300
# # Calculate the H-value
# H_SIR = I[magical_index] + S[magical_index] - p*log(S[magical_index])
# H_SIR_2 = I[200] + S[200] - p*log(S[200])
# # Take a point
# X0 = array([t[magical_index], S[magical_index], I[magical_index]])
# # Try to solve the ODE at hand
# Gamma_epsilon, infodict = integrate.odeint(dX_deps_SIR_S, X0, epsilon_vec, args = (H_SIR,r,p,round(S0/10)),full_output=True)
# # Split the solution into its component parts
# Gamma_S_t, Gamma_S_S, Gamma_S_I = Gamma_epsilon.T
# # Define a new time vector
# t_2 = linspace(Gamma_S_t[-1], t[-1], 250)
# # We need to integrate backwards in time as well to start at 0
# t_2_start = linspace(Gamma_S_t[-1], 0, 250)
# # Define new initial conditions for the transformed solutions
# X02 = array([Gamma_S_t[-1], Gamma_S_S[-1], Gamma_S_I[-1]])  
# X0_2 = array([Gamma_S_S[-1], Gamma_S_I[-1]])  
# # Solve the ODE at hand with the new initial conditions
# # Integrate backwards to get to 0
# X2_start, infodict = integrate.odeint(dX_dt_SIR, X0_2, t_2_start, args=(a,r), full_output=True)
# # Find the rest of the solution
# X2, infodict = integrate.odeint(dX_dt_SIR, X0_2, t_2, args=(a,r), full_output=True)
# infodict['message'] # >>> 'Integration successful.'
# # Split the solution into its component parts
# S_2, I_2 = X2.T
# # Split the start bit as well
# S_2_start, I_2_start = X2_start.T
# # Concatenate to get the full solution curves and the full time
# S_2_S = concatenate((flip(S_2_start,0), S_2), axis=0)
# I_2_S = concatenate((flip(I_2_start,0), I_2), axis=0)
# t_2_S = concatenate((flip(t_2_start,0), t_2), axis=0)
# # Plot the symmetry again for some other point on the solution curves
# magical_indices = [290, 300, 310, 320, 330, 340, 350, 360, 370]
# # Allocate memory for our lovely symmetry
# Gamma_S_t_vec = []
# Gamma_S_S_vec = []
# Gamma_S_I_vec = []
# # Loop over the indices and plot the symmetry transformation
# for index,new_magical_index in enumerate(magical_indices):
#     # Take a point
#     X0 = array([t[new_magical_index], S[new_magical_index], I[new_magical_index]])  
#     # Try to solve the ODE at hand
#     Gamma_epsilon_temp, infodict = integrate.odeint(dX_deps_SIR_S, X0, epsilon_vec, args = (H_SIR,r,p,round(S0/10)),full_output=True)    
#     # Split the solution into its component parts
#     Gamma_S_t_temp, Gamma_S_S_temp, Gamma_S_I_temp = Gamma_epsilon_temp.T    
#     # Append our solutions
#     Gamma_S_t_vec.append(Gamma_S_t_temp)    
#     Gamma_S_S_vec.append(Gamma_S_S_temp)    
#     Gamma_S_I_vec.append(Gamma_S_I_temp)
# # ---------------------------------------------------------------------------------
# # ---------------------------------------------------------------------------------
# # I-directional symmetry of the SIR model
# # ---------------------------------------------------------------------------------
# # ---------------------------------------------------------------------------------
# # Magical index
# magical_index = 300
# branch_index_temp = 0
# # Make a test value for the lower limit in the integral
# I0_test = 100
# # Take a point
# X0 = array([t[magical_index], S[magical_index], I[magical_index]])
# # Try to solve the ODE at hand
# Gamma_epsilon, infodict = integrate.odeint(dX_deps_SIR_I, X0, epsilon_vec, args = (H_SIR,r,p,I0_test,branch_index_temp),full_output=True)
# # Split the solution into its component parts
# Gamma_I_t, Gamma_I_S, Gamma_I_I = Gamma_epsilon.T
# # Define a new time vector
# t_2 = linspace(Gamma_I_t[-1], t[-1], 250)
# # We need to integrate backwards in time as well to start at 0
# t_2_start = linspace(Gamma_I_t[-1], 0, 250)
# # Define new initial conditions for the transformed solutions
# X02 = array([Gamma_I_t[-1], Gamma_I_S[-1], Gamma_I_I[-1]])  
# X0_2 = array([Gamma_I_S[-1], Gamma_I_I[-1]])  
# # Solve the ODE at hand with the new initial conditions
# # Integrate backwards to get to 0
# X2_start, infodict = integrate.odeint(dX_dt_SIR, X0_2, t_2_start, args=(a,r), full_output=True)
# # Find the rest of the solution
# X2, infodict = integrate.odeint(dX_dt_SIR, X0_2, t_2, args=(a,r), full_output=True)
# infodict['message'] # >>> 'Integration successful.'
# # Split the solution into its component parts
# S_2, I_2 = X2.T
# # Split the start bit as well
# S_2_start, I_2_start = X2_start.T
# # Concatenate to get the full solution curves and the full time
# S_2_I = concatenate((flip(S_2_start,0), S_2), axis=0)
# I_2_I = concatenate((flip(I_2_start,0), I_2), axis=0)
# t_2_I = concatenate((flip(t_2_start,0), t_2), axis=0)
# # Plot the symmetry again for some other point on the solution curves
# magical_indices_I = [240,250,260,270,280,290]
# branches_I = [0,0,0,0,0,0]
# #branches_I = [-1,-1,-1,-1,-1]
# # Allocate memory for our lovely symmetry
# Gamma_I_t_vec = []
# Gamma_I_S_vec = []
# Gamma_I_I_vec = []
# # Loop over the indices and plot the symmetry transformation
# for index,new_magical_index in enumerate(magical_indices_I):
#     # Take a point
#     X0 = array([t[new_magical_index], S[new_magical_index], I[new_magical_index]])  
#     # Try to solve the ODE at hand
#     Gamma_epsilon_temp, infodict = integrate.odeint(dX_deps_SIR_I, X0, epsilon_vec, args = (H_SIR,r,p,I0_test,branches_I[index]),full_output=True)    
#     # Split the solution into its component parts
#     Gamma_I_t_temp, Gamma_I_S_temp, Gamma_I_I_temp = Gamma_epsilon_temp.T    
#     # Append our solutions
#     Gamma_I_t_vec.append(Gamma_I_t_temp)    
#     Gamma_I_S_vec.append(Gamma_I_S_temp)    
#     Gamma_I_I_vec.append(Gamma_I_I_temp)
# #---------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------
# # Plot of symmetries for the SIR model in the time domain
# #---------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------
# # Plot the solutions
# f2, ax_2 = plt.subplots(1,2,constrained_layout=True, figsize=(20, 8))
# #---------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------
# # S-directional symmetry
# #---------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------
# # Original solution S
# ax_2[0].plot(t,S, '-', label="$S$" ,color=(0/256,68/256,27/256),linewidth=3.0)
# # Transformed solution S
# ax_2[0].plot(t_2_S,S_2_S, '-', label="$\\hat{S}$" ,color=(35/256,139/256,69/256),linewidth=3.0)
# # Original solution I
# ax_2[0].plot(t,I, '-', label="$I$" ,color=(77/256,0/256,75/256),linewidth=3.0)
# # Transformed solution I
# ax_2[0].plot(t_2_S,I_2_S, '-', label="$\\hat{I}$" ,color=(136/256,65/256,157/256),linewidth=3.0)
# # The symmetry without legend
# for index in range(len(Gamma_S_t_vec)):
#     if index == 0:
#         ax_2[0].plot(Gamma_S_t_vec[index],Gamma_S_S_vec[index], '--', label="Symmetry, $\\Gamma_{3,\\epsilon}^{\\mathrm{SIR},S}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
#         ax_2[0].plot(Gamma_S_t_vec[index],Gamma_S_I_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
#     else:
#         ax_2[0].plot(Gamma_S_t_vec[index],Gamma_S_S_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
#         ax_2[0].plot(Gamma_S_t_vec[index],Gamma_S_I_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
# # General axes properties
# ax_2[0].legend(loc='best',prop={"size":20})
# ax_2[0].set_xlabel(xlabel='Time, $t$ [days]',fontsize=25)
# ax_2[0].set_ylabel(ylabel='Population size [\\# individuals]',fontsize=25)
# ax_2[0].tick_params(axis='both', which='major', labelsize=20)
# ax_2[0].tick_params(axis='both', which='minor', labelsize=20)
# ax_2[0].grid()
# #---------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------
# # I-directional symmetry
# #---------------------------------------------------------------------------------------------------
# #---------------------------------------------------------------------------------------------------
# # Solution S
# ax_2[1].plot(t,S, '-', label="$S$" ,color=(0/256,68/256,27/256),linewidth=3.0)
# # Transformed solution S
# ax_2[1].plot(t_2_I,S_2_I, '-', label="$\\hat{S}$" ,color=(35/256,139/256,69/256),linewidth=3.0)
# # Solution I
# ax_2[1].plot(t,I, '-', label="$I$" ,color=(77/256,0/256,75/256),linewidth=3.0)
# # Transformed solution I
# ax_2[1].plot(t_2_I,I_2_I, '-', label="$\\hat{I}$" ,color=(136/256,65/256,157/256),linewidth=3.0)
# # The symmetry without legend
# for index in range(len(Gamma_I_t_vec)):
#     if index == 0:
#         ax_2[1].plot(Gamma_I_t_vec[index],Gamma_I_S_vec[index], '--', label="Symmetry, $\\Gamma_{3,\\epsilon}^{\\mathrm{SIR},I}$" ,color=(0/256,0/256,0/256),linewidth=3.0)
#         ax_2[1].plot(Gamma_I_t_vec[index],Gamma_I_I_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
#     else:
#         ax_2[1].plot(Gamma_I_t_vec[index],Gamma_I_S_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
#         ax_2[1].plot(Gamma_I_t_vec[index],Gamma_I_I_vec[index], '--',color=(0/256,0/256,0/256),linewidth=3.0)
# # General ticks settings etc.
# ax_2[1].tick_params(axis='both', which='major', labelsize=20)
# ax_2[1].tick_params(axis='both', which='minor', labelsize=20)
# ax_2[1].grid()
# ax_2[1].legend(loc='best',prop={"size":20})
# ax_2[1].set_xlabel(xlabel='Time, $t$ [days]',fontsize=25)
# ax_2[1].set_ylabel(ylabel='Population size [\\# individuals]',fontsize=25)
# #-----------------------------------------------------------------------------------------------------
# # Title and saving the figure
# #-----------------------------------------------------------------------------------------------------
# # Title and saving the figure
# f2.suptitle('SIR symmetries in the time domain',fontsize=30,weight='bold');
# f2.savefig('../Figures/time_domain_symmetries_SIR.png')
# plt.show()
# #plt.show()
# #=================================================================================
# #=================================================================================
# # Plot time domain symmetries in LaTeX as well...
# #=================================================================================
# #=================================================================================
# #--------------------------------------------------------------------------------
# # S-directional
# #--------------------------------------------------------------------------------
# for index in range(len(Gamma_S_t_vec)):
#     if index == 0:
#         plot_LaTeX_2D(Gamma_S_t_vec[index],Gamma_S_S_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\\mathrm{SIR},S}_{3,\\epsilon}$")
#         plot_LaTeX_2D(Gamma_S_t_vec[index],Gamma_S_I_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
#     else:
#         plot_LaTeX_2D(Gamma_S_t_vec[index],Gamma_S_S_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
#         if magical_indices[index]>320:
#             plot_LaTeX_2D(Gamma_S_t_vec[index],Gamma_S_I_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])                
# # Plot the solutions as well
# plot_LaTeX_2D(t,S,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=r_1,line width=1.5pt,","$S(t)$")
# plot_LaTeX_2D(t_2_S,S_2_S,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=r_2,line width=1.5pt,","$\\hat{S}(t)$")
# plot_LaTeX_2D(t,I,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=r_3,line width=1.5pt,","$I(t)$")
# plot_LaTeX_2D(t_2_S,I_2_S,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_S_time.tex","color=r_4,line width=1.5pt,","$\\hat{I}(t)$")
# #--------------------------------------------------------------------------------
# # I-directional
# #--------------------------------------------------------------------------------
# for index in range(len(Gamma_I_t_vec)):
#     if index == 0:
#         plot_LaTeX_2D(Gamma_I_t_vec[index],Gamma_I_S_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\\mathrm{SIR},I}_{3,\\epsilon}$")
#         plot_LaTeX_2D(Gamma_I_t_vec[index],Gamma_I_I_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])        
#     else:
#         plot_LaTeX_2D(Gamma_I_t_vec[index],Gamma_I_S_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
#         plot_LaTeX_2D(Gamma_I_t_vec[index],Gamma_I_I_vec[index],"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])                
# # Plot the solutions as well
# plot_LaTeX_2D(t,S,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=r_1,line width=1.5pt,","$S(t)$")
# plot_LaTeX_2D(t_2_I,S_2_I,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=r_2,line width=1.5pt,","$\\hat{S}(t)$")
# plot_LaTeX_2D(t,I,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=r_3,line width=1.5pt,","$I(t)$")
# plot_LaTeX_2D(t_2_I,I_2_I,"../Figures/LaTeX_figures/SIR_symmetries/Input/SIR_I_time.tex","color=r_4,line width=1.5pt,","$\\hat{I}(t)$")
