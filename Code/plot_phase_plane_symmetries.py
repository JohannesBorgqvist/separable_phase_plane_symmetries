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
from symmetry_toolbox import *
#=================================================================================
#=================================================================================
# Solving the ODEs for the LV model
#=================================================================================
#=================================================================================
t = linspace(0, 10, 500)              # time
X0 = array([1, 0.10])                     # initials conditions: 10 rabbits and 5 foxes
#X0_2 = array([20, 10])                     # initials conditions: 20 rabbits and 10 foxes
X1, infodict = integrate.odeint(dX_dt_LV, X0, t, args = (1,),full_output=True)
infodict['message']                     # >>> 'Integration successful.'

#!python
#X2, infodict = integrate.odeint(dX_dt_LV, X0, t, args = (2,),full_output=True)
#infodict['message']                     # >>> 'Integration successful.'

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
#u_transformed = asarray([u_transf(u_temp, epsilon) for u_temp in list(u)])
#v_transformed = asarray([v_transf(v_temp, epsilon,alpha) for v_temp in list(v)])
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
    trans_vec = [u_transf(u[u_index],epsilon_temp) for epsilon_temp in list(epsilon_vec)]
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
f1, ax_1 = p.subplots(1, 3, constrained_layout=True, figsize=(20, 8))
ax_1[0].plot(asarray(t_sym[0]),asarray([u[t_indices[0]]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--', label="$\\left.\\Gamma_{\\epsilon}^{\\tau}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,t_index in enumerate(list(t_indices)):
    ax_1[0].plot(asarray(t_sym[index]),asarray([u[t_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)
    ax_1[0].plot(asarray(t_sym[index]),asarray([v[t_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_1[0].plot(t, u, '-', label="$u(\\tau)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[0].plot(t, v  , '-', label='$v(\\tau)$',color=(153/256,216/256,201/256),linewidth=3.0)
ax_1[0].plot(t_trans, u, '-', label='$u(\\hat{\\tau})$',color=(77/256,0/256,75/256),linewidth=3.0)
ax_1[0].plot(t_trans, v  , '-', label='$v(\\hat{\\tau})$',color=(140/256,150/256,198/256),linewidth=3.0)
ax_1[0].grid()
ax_1[0].legend(loc='best',prop={"size":20})
ax_1[1].plot(u, v,'-', label='Original population, $(u,v)$',color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[1].plot(u_transformed_2,v_transformed_2, '-', label='Transformed population, $(\\hat{u},v)$',color=(77/256,0/256,75/256),linewidth=3.0)
ax_1[1].plot(asarray(u_sym[0]),asarray([v[u_indices[0]]*index/index for index in range(len(epsilon_vec))]), '--', label="$\\left.\\Gamma_{\\epsilon}^{u}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,u_index in enumerate(list(u_indices)):
    ax_1[1].plot(asarray(u_sym[index]),asarray([v[u_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_1[1].grid()
ax_1[1].legend(loc='best',prop={"size":20})
ax_1[2].plot(u, v,'-', label='Original population, $(u,v)$',color=(0/256,68/256,27/256),linewidth=3.0)
ax_1[2].plot(u_transformed,v_transformed, '-', label='Transformed population, $(u,\\hat{v})$',color=(77/256,0/256,75/256),linewidth=3.0)
ax_1[2].plot(asarray([u[v_indices[0]]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]), asarray(v_sym[0]), '--', label="$\\left.\\Gamma_{\\epsilon}^{v}\\right|_{\\epsilon=" +str(epsilon) + "}$" ,color=(0,0,0),linewidth=2.0)
for index,v_index in enumerate(list(v_indices)):
    ax_1[2].plot(asarray([u[v_index]*((index+1)/(index+1)) for index in range(len(epsilon_vec))]),asarray(v_sym[index]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_1[2].grid()
ax_1[2].legend(loc='best',prop={"size":20})
# Set fontsize of labels
ax_1[0].set_xlabel(xlabel='Time, $\\tau$',fontsize=25)
ax_1[0].set_ylabel(ylabel='Population size',fontsize=25)
ax_1[1].set_xlabel(xlabel='Rabbits, $u(\\tau)$',fontsize=25)
ax_1[1].set_ylabel(ylabel='Foxes, $v(\\tau)$',fontsize=25)
ax_1[2].set_xlabel(xlabel='Rabbits, $u(\\tau)$',fontsize=25)
ax_1[2].set_ylabel(ylabel='Foxes, $v(\\tau)$',fontsize=25)
# Change the size of the ticks
ax_1[0].tick_params(axis='both', which='major', labelsize=20)
ax_1[0].tick_params(axis='both', which='minor', labelsize=20)
ax_1[1].tick_params(axis='both', which='major', labelsize=20)
ax_1[1].tick_params(axis='both', which='minor', labelsize=20)
ax_1[2].tick_params(axis='both', which='major', labelsize=20)
ax_1[2].tick_params(axis='both', which='minor', labelsize=20)
# Title and saving the figure
f1.suptitle('Symmetries of the Lotka-Volterra model',fontsize=30,weight='bold')
f1.savefig('../Figures/LV_symmetries.png')
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
# Calculate the transformed solution
S_transformed = asarray([S_transf(S_temp, epsilon,a,r) for S_temp in list(S)])
# Modify S and I so that we only plot up until the threshold N
S_new = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I[index])<=N and I[index]>0])
I_new = asarray([I[index] for index,S_temp in enumerate(S) if (S_temp+I[index])<=N and I[index]>0])
# The same goes for the transformed solution
S_trans_new = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I_trans[index])<=N and I_trans[index]>0])
I_trans_new = asarray([I_trans[index] for index,S_temp in enumerate(S) if (S_temp+I_trans[index])<=N and I_trans[index]>0])
S_trans_new_2 = asarray([S_temp for index,S_temp in enumerate(S) if (S_temp+I_trans_2[index])<=N and I_trans_2[index]>0])
I_trans_new_2 = asarray([I_trans_2[index] for index,S_temp in enumerate(S) if (S_temp+I_trans_2[index])<=N and I_trans_2[index]>0])
#I_trans = asarray([I[index] for index,S_temp in enumerate(S_new) if (S_temp+I[index])<=N and I[index]>0])
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon-0.005,epsilon/50)
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
# Plot the solutions
f2, ax_2 = plt.subplots(1,1,constrained_layout=True, figsize=(20, 8))
#ax_2.plot(S,I, '--', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_2.plot(S_new,I_new, '-', label="Original solution, $(S,I)$" ,color=(0/256,68/256,27/256),linewidth=3.0)
ax_2.plot(S_trans_new,I_trans_new, '-', label="Transformed solution, $(\hat{S},I)$" ,color=(35/256,139/256,69/256),linewidth=3.0)
ax_2.plot(S_trans_new_2,I_trans_new_2, '-', label="Transformed solution, $(\hat{\hat{S}},I)$" ,color=(102/256,194/256,164/256),linewidth=3.0)
for index,S_index in enumerate(list(S_indices)):
    if index == 0:
        ax_2.plot(asarray(S_sym[index]),asarray([I_new[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),label="$\\Gamma^{\\mathrm{SIR},S}_{\\epsilon}$",linewidth=2.0)
    else:
        ax_2.plot(asarray(S_sym[index]),asarray([I_new[S_index]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)
#for index,S_val in enumerate(S_sym_2):
for index in range(len(S_sym_2)):
    ax_2.plot(asarray(S_sym_2[index]),asarray([I_trans_new[S_indices_2[index]]*((temp_index+1)/(temp_index+1)) for temp_index in range(len(epsilon_vec))]), '--' ,color=(0,0,0),linewidth=2.0)    
ax_2.plot(asarray([N, 0]),asarray([0, N]),label="Total population size, $N=I+S$",color=(0,0,0),linewidth=3.0)
ax_2.plot(asarray([0, p]),asarray([I_max, I_max]), '--' ,color=(0,0,0),linewidth=3.0)
ax_2.plot(asarray([0, p]),asarray([I_max+epsilon, I_max+epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_2.plot(asarray([0, p]),asarray([I_max+2*epsilon, I_max+2*epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_2.plot(asarray([p, p]),asarray([0, I_max+2*epsilon]), '--' ,color=(0,0,0),linewidth=3.0)
ax_2.tick_params(axis='both', which='major', labelsize=20)
ax_2.tick_params(axis='both', which='minor', labelsize=20)
ax_2.grid()
ax_2.legend(loc='best',prop={"size":20})
ax_2.set_xlabel(xlabel='Susceptibles, $S(t)$',fontsize=25)
ax_2.set_ylabel(ylabel='Infected, $I(t)$',fontsize=25)
f2.suptitle('Symmetries of the SIR model',fontsize=30,weight='bold')
f2.savefig('../Figures/SIR_symmetries.png')
#=================================================================================
#=================================================================================
# Transforming the linear model with our symmetries
#=================================================================================
#=================================================================================
# ---------------------------------------------------------------------------
# Overall properties (a,b,c,d)=(-1,0,0,-2) Stable node
# ---------------------------------------------------------------------------
# Define the parameters for the stable point
a_stable = -1
b_stable = 0
c_stable = 0
d_stable = -2
# Transformation parameters
epsilon_r=0.33
epsilon_theta = pi/6
# Define our initial condition
u0 = 1
v0 = 3
# Define the time vector
t = linspace(0, 2, 500)              # Time
X0 = array([u0, v0])                  # ICs
X1, infodict = integrate.odeint(dX_dt_linear, X0, t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the original solution with the defined parameters and initial conditions
u_1, v_1 = X1.T
# The indices we wish to transform
lin_indices = arange(0,len(t),15)
# ---------------------------------------------------------------------------
# Stable node (a,b,c,d)=(-1,0,0,-2) radial symmetry
# ---------------------------------------------------------------------------
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_r-0.005,epsilon_r/50)
# Allocate our empty lists
Gamma_trans_r_u_1 = []
Gamma_trans_r_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_ODE, array([u_1[lin_index],v_1[lin_index]]), epsilon_vec, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_r_u_temp, Gamma_trans_r_v_temp = X_r.T
    # Save our transformations
    Gamma_trans_r_u_1.append(Gamma_trans_r_u_temp)
    Gamma_trans_r_v_1.append(Gamma_trans_r_v_temp)    
# Solve to get the transformed angle solution
X1_r, infodict = integrate.odeint(dX_dt_linear, array([Gamma_trans_r_u_1[0][-1], Gamma_trans_r_v_1[0][-1]]), t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
u_r_1, v_r_1 = X1_r.T
# Allocate our empty lists
Gamma_trans_r_u_2 = []
Gamma_trans_r_v_2 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_ODE, array([u_r_1[lin_index],v_r_1[lin_index]]), epsilon_vec, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_r_u_temp, Gamma_trans_r_v_temp = X_r.T
    # Save our transformations
    Gamma_trans_r_u_2.append(Gamma_trans_r_u_temp)
    Gamma_trans_r_v_2.append(Gamma_trans_r_v_temp)    
# Solve to get the transformed angle solution
X2_r, infodict = integrate.odeint(dX_dt_linear, array([Gamma_trans_r_u_2[0][-1], Gamma_trans_r_v_2[0][-1]]), t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
u_r_2, v_r_2 = X2_r.T
# ---------------------------------------------------------------------------
# Stable node (a,b,c,d)=(-1,0,0,-2) angle symmetry
# ---------------------------------------------------------------------------
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_theta-0.005,epsilon_theta/50)
# Allocate our empty lists
Gamma_trans_theta_u_1 = []
Gamma_trans_theta_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_ODE, array([u_1[lin_index],v_1[lin_index]]), epsilon_vec, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_theta_u_temp, Gamma_trans_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_trans_theta_u_1.append(Gamma_trans_theta_u_temp)
    Gamma_trans_theta_v_1.append(Gamma_trans_theta_v_temp)    
# Solve to get the transformed angle solution
X1_angle, infodict = integrate.odeint(dX_dt_linear, array([Gamma_trans_theta_u_1[0][-1], Gamma_trans_theta_v_1[0][-1]]), t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
u_angle_1, v_angle_1 = X1_angle.T
# Allocate our empty lists
Gamma_trans_theta_u_2 = []
Gamma_trans_theta_v_2 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_ODE, array([u_angle_1[lin_index],v_angle_1[lin_index]]), epsilon_vec, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_theta_u_temp, Gamma_trans_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_trans_theta_u_2.append(Gamma_trans_theta_u_temp)
    Gamma_trans_theta_v_2.append(Gamma_trans_theta_v_temp)    
# Solve to get the transformed angle solution
X2_angle, infodict = integrate.odeint(dX_dt_linear, array([Gamma_trans_theta_u_2[0][-1], Gamma_trans_theta_v_2[0][-1]]), t, args = (a_stable,b_stable,c_stable,d_stable),full_output=True)
u_angle_2, v_angle_2 = X2_angle.T
# ---------------------------------------------------------------------------
# Overall properties (a,b,c,d)=(1,0,0,-2) Saddle point
# ---------------------------------------------------------------------------
# Define the parameters for the stable point
a_saddle = 1
b_saddle = 0
c_saddle = 0
d_saddle = -2
# Solve the system for the saddle parameters
X2, infodict = integrate.odeint(dX_dt_linear, X0, t, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
infodict['message']                     # >>> 'Integration successful.'
# Extract the original solution with the defined parameters and initial conditions
u_2, v_2 = X2.T
# ---------------------------------------------------------------------------
# Saddle point (a,b,c,d)=(1,0,0,-2) radial symmetry
# ---------------------------------------------------------------------------
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_r-0.005,epsilon_r/50)
# Allocate our empty lists
Gamma_saddle_r_u_1 = []
Gamma_saddle_r_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_ODE, array([u_2[lin_index],v_2[lin_index]]), epsilon_vec, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_saddle_r_u_temp, Gamma_saddle_r_v_temp = X_r.T
    # Save our transformations
    Gamma_saddle_r_u_1.append(Gamma_saddle_r_u_temp)
    Gamma_saddle_r_v_1.append(Gamma_saddle_r_v_temp)    
# Solve to get the transformed angle solution
X1_r, infodict = integrate.odeint(dX_dt_linear, array([Gamma_saddle_r_u_1[0][-1], Gamma_saddle_r_v_1[0][-1]]), t, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
u_saddle_r_1, v_saddle_r_1 = X1_r.T
# Allocate our empty lists
Gamma_saddle_r_u_2 = []
Gamma_saddle_r_v_2 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_r, infodict = integrate.odeint(Gamma_r_ODE, array([u_saddle_r_1[lin_index],v_saddle_r_1[lin_index]]), epsilon_vec, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_r_u_temp, Gamma_trans_r_v_temp = X_r.T
    # Save our transformations
    Gamma_saddle_r_u_2.append(Gamma_trans_r_u_temp)
    Gamma_saddle_r_v_2.append(Gamma_trans_r_v_temp)    
# Solve to get the transformed angle solution
X2_r, infodict = integrate.odeint(dX_dt_linear, array([Gamma_saddle_r_u_2[0][-1], Gamma_saddle_r_v_2[0][-1]]), t, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
u_saddle_r_2, v_saddle_r_2 = X2_r.T
# ---------------------------------------------------------------------------
# Stable node (a,b,c,d)=(-1,0,0,-2) angle symmetry
# ---------------------------------------------------------------------------
# Plot the action of the symmetry
epsilon_vec = arange(0,epsilon_theta-0.005,epsilon_theta/50)
# Allocate our empty lists
Gamma_saddle_theta_u_1 = []
Gamma_saddle_theta_v_1 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_ODE, array([u_2[lin_index],v_2[lin_index]]), epsilon_vec, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_theta_u_temp, Gamma_trans_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_saddle_theta_u_1.append(Gamma_trans_theta_u_temp)
    Gamma_saddle_theta_v_1.append(Gamma_trans_theta_v_temp)    
# Define a temporary time vector since we need to solve this ODE for a longer time
t_trans_1 = linspace(0, 2.75, 500)              # Time
# Solve to get the transformed angle solution
X1_angle, infodict = integrate.odeint(dX_dt_linear, array([Gamma_saddle_theta_u_1[0][-1], Gamma_saddle_theta_v_1[0][-1]]), t_trans_1, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
# Extract the solution
u_saddle_angle_1, v_saddle_angle_1 = X1_angle.T
# Allocate our empty lists
Gamma_saddle_theta_u_2 = []
Gamma_saddle_theta_v_2 = []
# Loop over our indices and find the transformations
for lin_index in lin_indices:
    # Get our transformations
    X_theta, infodict = integrate.odeint(Gamma_theta_ODE, array([u_saddle_angle_1[lin_index],v_saddle_angle_1[lin_index]]), epsilon_vec, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    # Extract transformations
    Gamma_trans_theta_u_temp, Gamma_trans_theta_v_temp = X_theta.T
    # Save our transformations
    Gamma_saddle_theta_u_2.append(Gamma_trans_theta_u_temp)
    Gamma_saddle_theta_v_2.append(Gamma_trans_theta_v_temp)    
# Define a temporary time vector since we need to solve this ODE for a longer time
t_trans_2 = linspace(0, 3.5, 500)              # Time
# Solve to get the transformed angle solution
X2_angle, infodict = integrate.odeint(dX_dt_linear, array([Gamma_saddle_theta_u_2[0][-1], Gamma_saddle_theta_v_2[0][-1]]), t_trans_2, args = (a_saddle,b_saddle,c_saddle,d_saddle),full_output=True)
# Extract the solution
u_saddle_angle_2, v_saddle_angle_2 = X2_angle.T
#=================================================================================
#=================================================================================
# Plot the linear system
#=================================================================================
#=================================================================================
#Define the first figure
f3, ax_3 = plt.subplots(2, 2, constrained_layout=True, figsize=(20, 8))
# Plot 1: Stable node with radial symmetry
ax_3[0][0].plot(u_1,v_1,color=(0/256,68/256,27/256),label="$(u,v)$",linewidth=4.0)
ax_3[0][0].plot(u_r_1,v_r_1,color=(0/256,109/256,44/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
ax_3[0][0].plot(u_r_2,v_r_2,color=(102/256,194/256,164/256),label="$(\\hat{\\hat{u}},\\hat{\\hat{v}})$",linewidth=4.0)
for index in range(len(Gamma_trans_r_u_1)):
    if index == 0:
        ax_3[0][0].plot(Gamma_trans_r_u_1[index], Gamma_trans_r_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{r}$",linewidth=0.5)
    else:
        ax_3[0][0].plot(Gamma_trans_r_u_1[index], Gamma_trans_r_v_1[index],'--',color=(0,0,0),linewidth=0.5)
for index in range(len(Gamma_trans_r_u_2)):
    ax_3[0][0].plot(Gamma_trans_r_u_2[index], Gamma_trans_r_v_2[index],'--',color=(0,0,0),linewidth=0.5)
ax_3[0][0].grid()
ax_3[0][0].legend(loc='best',prop={"size":20})
ax_3[0][0].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
ax_3[0][0].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
ax_3[0][0].tick_params(axis='both', which='major', labelsize=20)
ax_3[0][0].tick_params(axis='both', which='minor', labelsize=20)
# Plot 2: Stable node with angular symmetry
ax_3[0][1].plot(u_1,v_1,color=(103/256,0/256,31/256),label="$(u,v)$",linewidth=4.0)
ax_3[0][1].plot(u_angle_1,v_angle_1,color=(206/256,18/256,86/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
ax_3[0][1].plot(u_angle_2,v_angle_2,color=(223/256,101/256,176/256),label="$(\\hat{\\hat{u}},\\hat{\\hat{v}})$",linewidth=4.0)
for index in range(len(Gamma_trans_theta_u_1)):
    if index == 0:
        ax_3[0][1].plot(Gamma_trans_theta_u_1[index], Gamma_trans_theta_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{\\theta}$",linewidth=0.5)
    else:
        ax_3[0][1].plot(Gamma_trans_theta_u_1[index], Gamma_trans_theta_v_1[index],'--',color=(0,0,0),linewidth=0.5)
for index in range(len(Gamma_trans_theta_u_2)):
    ax_3[0][1].plot(Gamma_trans_theta_u_2[index], Gamma_trans_theta_v_2[index],'--',color=(0,0,0),linewidth=0.5)
ax_3[0][1].grid()
ax_3[0][1].legend(loc='best',prop={"size":20})
ax_3[0][1].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
ax_3[0][1].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
ax_3[0][1].tick_params(axis='both', which='major', labelsize=20)
ax_3[0][1].tick_params(axis='both', which='minor', labelsize=20)
# Plot 3: Saddle point with radial symmetry
ax_3[1][0].plot(u_2,v_2,color=(2/256,56/256,88/256),label="$(u,v)$",linewidth=4.0)
ax_3[1][0].plot(u_saddle_r_1,v_saddle_r_1,color=(5/256,112/256,176/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
ax_3[1][0].plot(u_saddle_r_2,v_saddle_r_2,color=(116/256,169/256,207/256),label="$(\\hat{\\hat{u}},\\hat{\\hat{v}})$",linewidth=4.0)
for index in range(len(Gamma_saddle_r_u_1)):
    if index == 0:
        ax_3[1][0].plot(Gamma_saddle_r_u_1[index], Gamma_saddle_r_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{r}$",linewidth=0.5)
    else:
        ax_3[1][0].plot(Gamma_saddle_r_u_1[index], Gamma_saddle_r_v_1[index],'--',color=(0,0,0),linewidth=0.5)
for index in range(len(Gamma_saddle_r_u_2)):
    ax_3[1][0].plot(Gamma_saddle_r_u_2[index], Gamma_saddle_r_v_2[index],'--',color=(0,0,0),linewidth=0.5)
ax_3[1][0].grid()
ax_3[1][0].legend(loc='best',prop={"size":20})
ax_3[1][0].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
ax_3[1][0].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
ax_3[1][0].tick_params(axis='both', which='major', labelsize=20)
ax_3[1][0].tick_params(axis='both', which='minor', labelsize=20)
# Plot 4: Saddle point with angular symmetry
ax_3[1][1].plot(u_2,v_2,color=(102/256,37/256,6/256),label="$(u,v)$",linewidth=4.0)
ax_3[1][1].plot(u_saddle_angle_1,v_saddle_angle_1,color=(204/256,76/256,2/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
ax_3[1][1].plot(u_saddle_angle_2,v_saddle_angle_2,color=(254/256,153/256,41/256),label="$(\\hat{\\hat{u}},\\hat{\\hat{v}})$",linewidth=4.0)
for index in range(len(Gamma_saddle_theta_u_1)):
    if index == 0:
        ax_3[1][1].plot(Gamma_saddle_theta_u_1[index], Gamma_saddle_theta_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{\\theta}$",linewidth=0.5)
    else:
        ax_3[1][1].plot(Gamma_saddle_theta_u_1[index], Gamma_saddle_theta_v_1[index],'--',color=(0,0,0),linewidth=0.5)
for index in range(len(Gamma_saddle_theta_u_2)):
    ax_3[1][1].plot(Gamma_saddle_theta_u_2[index], Gamma_saddle_theta_v_2[index],'--',color=(0,0,0),linewidth=0.5)
ax_3[1][1].grid()
ax_3[1][1].legend(loc='best',prop={"size":20})
ax_3[1][1].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
ax_3[1][1].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
ax_3[1][1].tick_params(axis='both', which='major', labelsize=20)
ax_3[1][1].tick_params(axis='both', which='minor', labelsize=20)
# We have a title of this figure as well
f3.suptitle('Symmetries of the linear model',fontsize=30,weight='bold')
f3.savefig('../Figures/symmetries_linear_model.png')


#=================================================================================
#=================================================================================
# Transforming the biological oscillator
#=================================================================================
#=================================================================================
# ---------------------------------------------------------------------------
# Overall properties
# ---------------------------------------------------------------------------
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
# Plot the oscillatory system
#=================================================================================
#=================================================================================
#Define the first figure
f4, ax_4 = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
# Plot 1: Radial symmetry on the biological oscillator
ax_4[0].plot(u_3,v_3,color=(0/256,68/256,27/256),label="$(u,v)$",linewidth=4.0)
ax_4[0].plot(u_os_r_2,v_os_r_2,color=(0/256,109/256,44/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
for index in range(len(Gamma_os_r_u_1)):
    if index == 0:
        ax_4[0].plot(Gamma_os_r_u_1[index], Gamma_os_r_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{r}$",linewidth=0.5)
    else:
        ax_4[0].plot(Gamma_os_r_u_1[index], Gamma_os_r_v_1[index],'--',color=(0,0,0),linewidth=0.5)
ax_4[0].grid()
ax_4[0].legend(loc='best',prop={"size":20})
ax_4[0].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
ax_4[0].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
ax_4[0].tick_params(axis='both', which='major', labelsize=20)
ax_4[0].tick_params(axis='both', which='minor', labelsize=20)
# Plot 2: Angular symmetry on the biological oscillator
ax_4[1].plot(u_3,v_3,color=(103/256,0/256,31/256),label="$(u,v)$",linewidth=4.0)
ax_4[1].plot(u_os_theta_2,v_os_theta_2,color=(206/256,18/256,86/256),label="$(\\hat{u},\\hat{v})$",linewidth=4.0)
for index in range(len(Gamma_os_theta_u_1)):
    if index == 0:
        ax_4[1].plot(Gamma_os_theta_u_1[index], Gamma_os_theta_v_1[index],'--',color=(0,0,0),label="$\\Gamma_{\\epsilon}^{\\theta}$",linewidth=0.5)
    else:
        ax_4[1].plot(Gamma_os_theta_u_1[index], Gamma_os_theta_v_1[index],'--',color=(0,0,0),linewidth=0.5)
ax_4[1].grid()
ax_4[1].legend(loc='best',prop={"size":20})
ax_4[1].set_xlabel(xlabel="$u(\\tau)$",fontsize=25)
ax_4[1].set_ylabel(ylabel="$v(\\tau)$",fontsize=25)
ax_4[1].tick_params(axis='both', which='major', labelsize=20)
ax_4[1].tick_params(axis='both', which='minor', labelsize=20)
# We have a title of this figure as well
f4.suptitle('Symmetries of the biological oscillator',fontsize=30,weight='bold')
f4.savefig('../Figures/symmetries_biological_oscillator.png')
# Show the plot in the end
plt.show()


