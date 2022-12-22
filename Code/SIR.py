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
I_trans = asarray([(C+epsilon)-S_temp +p*log(S_temp) for S_temp in S])
I_trans_2 = asarray([(C+2*epsilon)-S_temp +p*log(S_temp) for S_temp in S])
# Calculate the I directional symmetry as well
I_trans_I_dir = asarray([(C-epsilon)-S_temp +p*log(S_temp) for S_temp in S])
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
epsilon_vec = arange(0,(epsilon)-0.005,(epsilon)/50)
# Transformations S
S_sym = []
S_indices = arange(35,65,5)
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
    trans_vec = [I_new[I_index]-epsilon_temp for epsilon_temp in list(epsilon_vec)]
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
#ax_1[0].plot(S_trans_new_2,I_trans_new_2, '-', label="Transformed solution, $(\hat{\hat{S}},I)$" ,color=(102/256,194/256,164/256),linewidth=3.0)
for index,S_index in enumerate(list(S_indices)):
    if index == 0:
        ax_1[0].plot(asarray(S_sym[index]),asarray([I_new[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),label="$\\Gamma^{\\mathrm{SIR},S}_{\\epsilon}$",linewidth=2.0)
    else:
        ax_1[0].plot(asarray(S_sym[index]),asarray([I_new[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)
#for index,S_val in enumerate(S_sym_2):
# for index in range(len(S_sym_2)):
#     ax_1[0].plot(asarray(S_sym_2[index]),asarray([I_trans_new[S_indices_2[index]]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)
    
ax_1[0].plot(asarray([N, 0]),asarray([0, N]),label="Total population size, $N=I+S$",color=(0,0,0),linewidth=3.0)
ax_1[0].plot(asarray([0, p]),asarray([I_max, I_max]), '--' ,color=(0,0,0),linewidth=3.0)
#ax_1[0].plot(asarray([0, p]),asarray([I_max+epsilon, I_max+epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[0].plot(asarray([0, p]),asarray([I_max+epsilon, I_max+epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[0].plot(asarray([p, p]),asarray([0, I_max+epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
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
ax_1[1].plot(asarray([0, p]),asarray([I_max-epsilon, I_max-epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[1].plot(asarray([p, p]),asarray([0, I_max]), '--' ,color=(0,0,0),linewidth=3.0)
ax_1[1].plot(asarray([p, p]),asarray([0, I_max-epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
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
#plt.show()
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
plot_LaTeX_2D(S_new,I_new,"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=clr_1,line width=2.0pt,","$I(S)$")
plot_LaTeX_2D(S_trans_new,I_trans_new,"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=clr_2,line width=2.0pt,","$\\hat{I}(S;\\epsilon)$")
# The symmetry
for index,S_index in enumerate(list(S_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray(S_sym[index]),asarray([((I_new[S_index]*(temp_index+1))/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\mathrm{SIR},S}_{2,\\epsilon}$")
    elif index<5:
        plot_LaTeX_2D(asarray(S_sym[index]),asarray([((I_new[S_index]*(temp_index+1))/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=black,->,>=latex,densely dashed",[])
# for index in range(len(S_sym_2)):
#     plot_LaTeX_2D(asarray(S_sym_2[index]),asarray([I_trans_new[S_indices_2[index]]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt",[])
plot_LaTeX_2D([N*x for x in linspace(0,1,50)],[N*(1-x) for x in linspace(0,1,50)],"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=black,mark=diamond*,only marks,line width=0.75pt,","$N=I+S$")
plot_LaTeX_2D([0, p],[I_max,I_max],"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
plot_LaTeX_2D([0, p],[I_max+epsilon,I_max+epsilon],"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
#plot_LaTeX_2D([0, p],[I_max+2*epsilon,I_max+2*epsilon],"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
plot_LaTeX_2D([p, p],[0,I_max+epsilon],"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
#plot_LaTeX_2D(S_trans_new_2,I_trans_new_2,"../Figures/LaTeX_figures/SIR_symmetry/Input/S_trans.tex","color=clr_3,line width=2.0pt,","$\\hat{\\hat{I}}(S;\\epsilon)=\\hat{I}(S;2\\epsilon)$")
#=================================================================================
#=================================================================================
# I transformation
#=================================================================================
#=================================================================================
plot_LaTeX_2D(S_new,I_new,"../Figures/LaTeX_figures/SIR_symmetry/Input/I_trans.tex","color=clr_1,line width=2.0pt,","$I(S)$")
plot_LaTeX_2D(S_trans_I_dir_new,I_trans_I_dir_new,"../Figures/LaTeX_figures/SIR_symmetry/Input/I_trans.tex","color=clr_2,line width=2.0pt,","$\\hat{I}(S;\\epsilon)$")
# Plot the transformation
for index,I_index in enumerate(list(I_indices)):
    if index == 0:
        plot_LaTeX_2D(asarray([S_new[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]),"../Figures/LaTeX_figures/SIR_symmetry/Input/I_trans.tex","color=black,->,>=latex,densely dashed,line width=1.0pt","$\\Gamma^{\mathrm{SIR},I}_{2,\\epsilon}$")
    elif index<5:
        plot_LaTeX_2D(asarray([S_new[I_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]),asarray(I_sym[index]),"../Figures/LaTeX_figures/SIR_symmetry/Input/I_trans.tex","color=black,->,>=latex,densely dashed",[])
plot_LaTeX_2D([N*x for x in linspace(0,1,50)],[N*(1-x) for x in linspace(0,1,50)],"../Figures/LaTeX_figures/SIR_symmetry/Input/I_trans.tex","color=black,mark=diamond*,only marks,line width=0.75pt,","$N=I+S$")
plot_LaTeX_2D([0, p],[I_max,I_max],"../Figures/LaTeX_figures/SIR_symmetry/Input/I_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
plot_LaTeX_2D([0, p],[I_max-epsilon,I_max-epsilon],"../Figures/LaTeX_figures/SIR_symmetry/Input/I_trans.tex","color=black,densely dashed,line width=1.0pt,",[])
plot_LaTeX_2D([p, p],[0,I_max],"../Figures/LaTeX_figures/SIR_symmetry/Input/I_trans.tex","color=black,densely dashed,line width=1.0pt,",[])        


#=================================================================================
#=================================================================================
# Plot the lifted symmetries in the time domain
#=================================================================================
#=================================================================================
t = linspace(0, 0.05, 500)              # time
# initials conditions: 10 rabbits and 5 foxes
X0 = array([S0, I0])
X1, infodict = integrate.odeint(dX_dt_SIR, X0, t, args=(a,r), full_output=True)
infodict['message']                     # >>> 'Integration successful.'
#!python
S, I = X1.T




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
ax_2[0].plot(t,S, '-', label="$S$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_2[0].plot(t,I, '-', label="$I$" ,color=(35/256,139/256,69/256),linewidth=3.0)
ax_2[0].legend(loc='best',prop={"size":20})
ax_2[0].set_xlabel(xlabel='Time, $t$',fontsize=25)
ax_2[0].set_ylabel(ylabel='State',fontsize=25)
ax_2[0].tick_params(axis='both', which='major', labelsize=20)
ax_2[0].tick_params(axis='both', which='minor', labelsize=20)
ax_2[0].grid()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# I-directional symmetry
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
ax_2[1].plot(t,S, '-', label="$S$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_2[1].plot(t,I, '-', label="$I$" ,color=(35/256,139/256,69/256),linewidth=3.0)
ax_2[1].tick_params(axis='both', which='major', labelsize=20)
ax_2[1].tick_params(axis='both', which='minor', labelsize=20)
ax_2[1].grid()
ax_2[1].legend(loc='best',prop={"size":20})
ax_2[1].set_xlabel(xlabel='Time, $t$',fontsize=25)
ax_2[1].set_ylabel(ylabel='State',fontsize=25)
#-----------------------------------------------------------------------------------------------------
# Title and saving the figure
#-----------------------------------------------------------------------------------------------------
# Title and saving the figure
f2.suptitle('SIR symmetries in the time domain',fontsize=30,weight='bold');
f2.savefig('../Figures/time_domain_symmetries_SIR.png')
plt.show()
#plt.show()
