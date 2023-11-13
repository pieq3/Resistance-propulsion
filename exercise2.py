"""
Corso di Architetttura Navale
Analisi delle prove di rimorchio secondo il
Metodo di trasferimento modello-vero: ITTC'78
Last rev.:04/11/2023
"""

# Libraries
#---------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.interpolate
from scipy.stats import linregress

import os                                               # to place the workspace in the script's directory
path = os.getcwd()                                      # ""
os.chdir('../')                                         # ""
os.chdir(os.path.dirname(os.path.abspath(__file__)))    # ""


# 01 - Ship's characteristics
# --------------------------------------------------------------------------------------------------------------
L_S = 142.00            # length at waterline, (m)
S = 2949.5              # wetted surface area, (m^2)
B = 18.9                # breadth, (m)
T = 6.16                # draft, (m)
Nabla = 8425.4          # hull volume, (m^3)
Cb = Nabla/(L_S*B*T)    # block coefficient, (-)


# 02 - Loading experimental data
# --------------------------------------------------------------------------------------------------------------
# Loading the data from basin test in the "vasca" matrix 
vasca = np.loadtxt('prove_rimorchio_insean2340.txt', delimiter=' ')

# In the first column there are the recorded velocities, V_M
# in the second column are the corresponding resistance values, RT_M

V_M = vasca[:,0]
RT_M = vasca[:,1]       # resistance is given in kgf


# 03 - MODEL data processing
# --------------------------------------------------------------------------------------------------------------
Lambda = 24.824         # model-ship ratio
t_M = 18.3              # average temperature during testing, (degrees Celsius)
rho_M = 998.543         # density; WARNING! Make sure to update viscosity and density values whenever temperature changes
nu_M = 0.0000010462     # kinematic viscosity, (m^2/s); WARNING! Same as previous line
g = 9.81                # gravity acceleration, (m/s^2)
RT_M = RT_M*g           # resistance expressed in (N)

L_M = L_S/Lambda
S_M = S/(Lambda**2)

Fn_M = V_M/(np.sqrt(L_M*g))
Rn_M = V_M*L_M/nu_M
CT_M = RT_M/(0.5*rho_M*S_M*vasca[:,0]**2)
CF_M = 0.075/((np.log10(Rn_M)-2)**2)


# 03 - Prohaska's method
# --------------------------------------------------------------------------------------------------------------
Fn_aux = []              # preallocating vector for Fn_M <= 0.2 (as ITTC recommends)
CF_aux = []
CT_aux = []
j = 0

for i in range(0,len(Fn_M)):
    if (Fn_M[i] >= 0.12) and (Fn_M[i] <= 0.2):
        Fn_aux.append(Fn_M[i])
        CF_aux.append(CF_M[i])
        CT_aux.append(CT_M[i])
        j = j+1

Fn_aux = np.array(Fn_aux)
CF_aux = np.array(CF_aux) 
CT_aux = np.array(CT_aux)

x_axis = Fn_aux**4/CF_aux
y_axis = CT_aux/CF_aux
c, oneplusk, r_value, p_value, std_err = linregress(x_axis, y_axis)

x = np.linspace(min(x_axis), max(x_axis))
fig = plt.figure(figsize=(8,4))              
ax = fig.add_subplot(111)
ax.plot(x, c*x + oneplusk, color=(0.8500, 0.3250, 0.0980))
#ax.plot(x, c_alt*x+1+k_alt, color=(0, 0.4470, 0.7410))
plt.scatter(x_axis, y_axis, color=(0.8500, 0.3250, 0.0980), marker='x')
ax.grid(True)

plt.savefig('oneplusk.png', dpi=800, transparent=True)
#plt.show()


# 03 - Alternatives to Prohaska's method
# --------------------------------------------------------------------------------------------------------------
# empirical formulae, alternative to Prohaska method
# Wattanabe's formula
k_Wattanabe = -0.0095+25.6*((Nabla)/(L_S*B*T))/((L_S/B)**2*np.sqrt(B/T))

# Granwille's formula
k_Granwille = -0.03 + 32.8 * (Cb**2 / (L_S/B)**2 * (B/T))

# Prohaska's statistical formula
k_Prohaska = 0.11 + 0.128*(B/T) - 0.0157*(B/T)**2 -3.10*(Cb/(L_S/B)) +28.8*( (Cb/(L_S/B)))**2


# 04 - SHIP data processing
# --------------------------------------------------------------------------------------------------------------
T_S = 15                # average operating temperature, (degrees Celsius)
rho_S = 1026.021        # seawater density, (kg/m^3)
nu_S = 0.0000011892     # kinematic viscosity, (m^2/s)
A_T = 400               # transversal windage area, (m^2)
k_MAA = 150*10**(-6)    # recommended ITTC coefficient, (m)

V_S = V_M*np.sqrt(Lambda)
Fn_S = V_S/(np.sqrt(L_S*g))
Rn_S = V_S*L_S/nu_S
CF_S = 0.075/((np.log10(Rn_S)-2)**2)
CR = CT_M - oneplusk*CF_M
DeltaCF = (105*(k_MAA/L_S)**(1/3)-0.64)*10**(-3)
C_AA = 0.001*A_T/S
CT_S = oneplusk*CF_S + CR + DeltaCF + C_AA
CT_S_Wattanabe = (1+k_Wattanabe)*CF_S + CR + DeltaCF + C_AA
CT_S_Granville = (1+k_Granwille)*CF_S + CR + DeltaCF + C_AA
CT_S_Prohaska = (1+k_Prohaska)*CF_S + CR + DeltaCF + C_AA

RT_S = CT_S * 0.5 * rho_S * S * V_S**2
RT_S_Wattanabe = CT_S_Wattanabe * 0.5 * rho_S * S * V_S**2
RT_S_Granville = CT_S_Granville * 0.5 * rho_S * S * V_S**2
RT_S_Prohaska = CT_S_Prohaska * 0.5 * rho_S * S * V_S**2
P_E = RT_S * V_S


# 05 - Diagrams plotting
# --------------------------------------------------------------------------------------------------------------
# diagram 1 - MODEL data processing: CF_M,CT_M and (1+k)CF_M as functions of Rn_M
# diagram 2 - MODEL data processing: RT_M as a function of V_M
# diagram 3 - SHIP data processing: CF_S,CT_S and (1+k)CF_S as functions of Rn_S
# diagram 4 - SHIP data processing: RT_S as a function of V_S
# diagram 5 - SHIP data processing: P_E as a function of V_S
# diagram 6 - SHIP data processing: RT_S calculated in different ways


# diagram 1 - MODEL data processing: CF_M and CT_M as functions of Rn_M
fig = plt.figure(figsize=(8,4))              #figsize=(8,4)
ax = fig.add_subplot(111)
ax.plot(Rn_M, CT_M, color=(0.8500, 0.3250, 0.0980), label='$C_{T_M}$')
ax.plot(Rn_M, CF_M, color=(0, 0.4470, 0.7410), label='$C_{F_M}=\dfrac{0.075}{(log_{10}Rn_M-2)^2}$')
ax.plot(Rn_M, oneplusk*CF_M, color=(0.051 , 0.251 , 0.051), label='$(1+k)C_{F_M}$')
#plt.scatter(Rn_M, CT_M, color=(0.8500, 0.3250, 0.0980), marker='x')                                 #label='Punti',

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('MODEL data processing: $C_{F_M}$ and $C_{T_M}$ as functions of $Rn_M$')
ax.set_ylabel('$C_{T_M}$, $C_{F_M}$, $(1+k)C_{F_M}$ (-)', fontsize='12')
ax.set_xlabel('$Rn_M$ (-)', fontsize='12')
ax.grid(True)

plt.legend()
plt.savefig('CTMCFM.png', dpi=800, transparent=True)
#plt.show()


# diagram 2 - MODEL data processing: RT_M as a function of V_M
fig2 = plt.figure(figsize=(8,4))              
ax = fig2.add_subplot(111)
ax.plot(V_M,RT_M/1000, color=(0.8500, 0.3250, 0.0980), label='$R_{T_M}$')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('MODEL data processing: $R_{T_M}$ as a function of $V_M$')
ax.set_ylabel('$R_{T_M}$ (kN)', fontsize='12')
ax.set_xlabel('$V_M$ (m/s)', fontsize='12')
ax.grid(True)

#plt.legend()
plt.savefig('RTM.png', dpi=800, transparent=True)
#plt.show()


# diagram 3 - SHIP data processing: CF_S and CT_S as functions of Rn_S
fig = plt.figure(figsize=(8,4))              #figsize=(8,4)
ax = fig.add_subplot(111)
ax.plot(Rn_S,CT_S, color=(0.8500, 0.3250, 0.0980), label='$C_{T_S}$')
ax.plot(Rn_S,CF_S, color=(0, 0.4470, 0.7410), label='$C_{F_S}=\dfrac{0.075}{(log_{10}Rn_S-2)^2}$')
ax.plot(Rn_S, oneplusk*CF_S, color=(0.051 , 0.251 , 0.051), label='$(1+k)C_{F_M}$')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('SHIP data processing: $C_{F_S}$ and $C_{T_S}$ as functions of $Rn_S$')
ax.set_ylabel('$C_{T_S}$, $C_{F_S}$, $(1+k)C_{F_M}$ (-)', fontsize='12')
ax.set_xlabel('$Rn_S$ (-)', fontsize='12')
ax.grid(True)

plt.legend()
plt.savefig('CTSCFS.png', dpi=800, transparent=True)
#plt.show()


# diagram 4 - SHIP data processing: RT_S as a function of V_S
fig2 = plt.figure(figsize=(8,4))              
ax = fig2.add_subplot(111)
ax.plot(V_S*1.94384,RT_S/1000, color=(0.8500, 0.3250, 0.0980), label='$R_{T_S}$')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('SHIP data processing: $R_{T_S}$ as a function of $V_S$')
ax.set_ylabel('$R_{T_S}$ (kN)', fontsize='12')
ax.set_xlabel('$V_S$ (kn)', fontsize='12')
ax.grid(True)

#plt.legend()
plt.savefig('RTS.png', dpi=800, transparent=True)
#plt.show()


# diagram 5 - SHIP data processing: P_E as a function of V_S
fig2 = plt.figure(figsize=(8,4))              
ax = fig2.add_subplot(111)
ax.plot(V_S*1.94384,P_E, color=(0.8500, 0.3250, 0.0980), label='$P_E$')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('SHIP data processing: $R_{T_S}$ as a function of $V_S$')
ax.set_ylabel('$P_E$ (kW)', fontsize='12')
ax.set_xlabel('$V_S$ (kn)', fontsize='12')
ax.grid(True)

#plt.legend()
plt.savefig('P_E.png', dpi=800, transparent=True)
#plt.show()


# diagram 6 - SHIP data processing: RT_S calculated in different ways
fig2 = plt.figure(figsize=(10,6))              
ax = fig2.add_subplot(111)
ax.plot(V_S,RT_S/1000, color=(0.8500, 0.3250, 0.0980), label='$R_{T_S}$')
ax.plot(V_S,RT_S_Wattanabe/1000, color=(0, 0.4470, 0.7410), label='$R_{T_SWattanabe}$')
ax.plot(V_S,RT_S_Granville/1000, color=(0.051 , 0.251 , 0.051), label='$R_{T_SGranville}$')
ax.plot(V_S,RT_S_Prohaska/1000, color=(1 , 0.839 , 0.071), label='$R_{T_SProhaska}$')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('SHIP data processing: RT_S calculated in different ways')
ax.set_ylabel('$R_T$ (kN)', fontsize='12')
ax.set_xlabel('$V_S$ (kn)', fontsize='12')
#ax.set_yscale('log')
ax.grid(True)

plt.legend()
plt.savefig('RTSconfronto.png', dpi=800, transparent=True)
plt.show()


# 06 - DATA EXPORT
# --------------------------------------------------------------------------------------------------------------
# data export section, making them ready to be shown in a report; LaTeX optimised
auxmatrix_M = np.empty((len(V_M), 5))
auxmatrix_S = np.empty((len(V_S), 8))
auxmatrix_RT = np.empty((len(V_S), 4))


# For cycle to create a matrix which contains the most significant data
for i in range(0,len(V_M)):
    auxmatrix_M[i] = np.mat([V_M[i] , Fn_M[i] , Rn_M[i] , CF_M[i] , CT_M[i]])
    auxmatrix_S[i] = np.matrix([V_S[i] , Fn_S[i] , Rn_S[i] , CF_S[i] , CR[i] , CT_S[i] , RT_S[i] , P_E[i]])
    auxmatrix_RT[i] = np.mat([RT_S[i] , RT_S_Wattanabe[i] , RT_S_Granville[i] , RT_S_Prohaska[i] ])


# Export files
with open('eserc1Model.txt' , 'w') as f:                                        # creating the new file, watch for the build directory 
    for i in range(0,len(V_M)):
        f.write(np.array2string( auxmatrix_M[i],                                # data source
                                formatter={'float_kind':'{:13.3e}'.format},     # data precision, scientific notation
                                floatmode='fixed',                              # imposing same precision for all data
                                separator=' & '))                               # '&' separator, ready-to-use for a LaTeX table
        f.write('\n')

with open('eserc1Ship.txt' , 'w') as f:
    for i in range(0,len(V_M)):
        f.write(np.array2string( auxmatrix_S[i], 
                                formatter={'float_kind':'{:13.3e}'.format}, 
                                floatmode='fixed', 
                                separator=' & ')) 
        f.write('\n')

with open('eserc1ShipRes.txt' , 'w') as f:
    for i in range(0,len(V_M)):
        f.write(np.array2string( auxmatrix_RT[i], 
                                formatter={'float_kind':'{:13.3e}'.format}, 
                                floatmode='fixed', 
                                separator=' & ')) 
        f.write('\n')
