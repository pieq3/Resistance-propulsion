"""
Corso di Architetttura Navale
Analisi delle prove di rimorchio secondo il
Metodo di trasferimento modello-vero: ITTC'57
Last rev.:22/10/2023
"""

# Libraries
#---------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os                                               # to place the workspace in the script's directory
path = os.getcwd()                                      # ""
os.chdir('../')                                         # ""
os.chdir(os.path.dirname(os.path.abspath(__file__)))    # ""


# 01 - Ship's characteristics
# --------------------------------------------------------------------------------------------------------------
L_S = 159.47        # length at waterline, [m]
S = 4605            # wetted surface area, [m^2]


# 02 - Loading experimental data
# --------------------------------------------------------------------------------------------------------------
# Loading the data from basin test in the "vasca" matrix 
vasca = np.loadtxt('prove_rimorchio.txt', delimiter=' ')

'''#imposto percorso assoluto, così il file viene riconosciuto
import os

# Ottengo il percorso assoluto della directory corrente
current_directory = os.path.abspath(os.path.dirname(__file__))

# Costruisco il percorso completo al file
file_path = os.path.join(current_directory, "prove_rimorchio.txt")

# Carico i dati dal file
vasca = np.loadtxt(file_path, delimiter=' ')'''

# In the first column there are the recorded velocities, V_M
# in the second column are the corresponding resistance values, RT_M

V_M = vasca[:,0]
RT_M = vasca[:,1]


# 03 - MODEL data processing
# --------------------------------------------------------------------------------------------------------------
Lambda = 24.175         # model-ship ratio
T_M = 18.3              # average temperature during testing, (degrees Celsius)
rho_M = 998.543         # density; WARNING! Make sure to update viscosity and density values whenever temperature changes
nu_M = 0.0000010462     # kinematic viscosity, (m^2/s); WARNING! Same as previous line
g = 9.81                # gravity acceleration, (m/s^2)

L_M = L_S/Lambda
S_M = S/(Lambda**2)

Fn_M = V_M/(np.sqrt(L_M*g))
Rn_M = V_M*L_M/nu_M
CT_M = RT_M/(0.5*rho_M*S_M*vasca[:,0]**2)
CF_M = 0.075/((np.log10(Rn_M)-2)**2)
CR_M = CT_M - CF_M


# 04 - SHIP data processing
# --------------------------------------------------------------------------------------------------------------
T_S = 15                # average operating temperature, (degrees Celsius)
rho_S = 1026.021        # seawater density, (kg/m^3)
nu_S = 0.0000011892     # kinematic viscosity, (m^2/s)

V_S = V_M*np.sqrt(Lambda)
Fn_S = V_S/(np.sqrt(L_S*g))
Rn_S = V_S*L_S/nu_S
CF_S = 0.075/((np.log10(Rn_S)-2)**2)
CR_S = CR_M
DeltaCF = 0.0002
CT_S = CF_S + CR_S + DeltaCF
RT_S = CT_S * 0.5 * rho_S * S * V_S**2
P_E = RT_S * V_S


# 05 - Diagrams plotting
# --------------------------------------------------------------------------------------------------------------
# diagram 1 - MODEL data processing: CF_M and CT_M as functions of Rn_M
# diagram 2 - MODEL data processing: RT_M as a function of V_M
# diagram 3 - SHIP data processing: CF_S and CT_S as functions of Rn_S
# diagram 4 - SHIP data processing: RT_S as a function of V_S
# diagram 5 - SHIP data processing: P_E as a function of V_S


# diagram 1 - MODEL data processing: CF_M and CT_M as functions of Rn_M
fig = plt.figure(figsize=(8,4))              #figsize=(8,4)
ax = fig.add_subplot(111)
ax.plot(Rn_M,CT_M, color=(0.8500, 0.3250, 0.0980), label='$C_{T_M}$')
ax.plot(Rn_M,CF_M, color=(0, 0.4470, 0.7410), label='$C_{F_M}=\dfrac{0.075}{(log_{10}Rn_M-2)^2}$')
#plt.scatter(Rn_M, CT_M, color=(0.8500, 0.3250, 0.0980), marker='x')                                 #label='Punti',

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('MODEL data processing: $C_{F_M}$ and $C_{T_M}$ as functions of $Rn_M$')
ax.set_ylabel('$C_{T_M}$, $C_{F_M}$ (-)', fontsize='12')
ax.set_xlabel('$Rn_M$ (-)', fontsize='12')
ax.grid(True)

plt.legend()
plt.savefig('CTMCFM.png', dpi=1400, transparent=True)
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
plt.savefig('RTM.png', dpi=1400, transparent=True)
#plt.show()


# diagram 3 - SHIP data processing: CF_S and CT_S as functions of Rn_S
fig = plt.figure(figsize=(8,4))              #figsize=(8,4)
ax = fig.add_subplot(111)
ax.plot(Rn_S,CT_S, color=(0.8500, 0.3250, 0.0980), label='$C_{T_S}$')
ax.plot(Rn_S,CF_S, color=(0, 0.4470, 0.7410), label='$C_{F_S}=\dfrac{0.075}{(log_{10}Rn_S-2)^2}$')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('SHIP data processing: $C_{F_S}$ and $C_{T_S}$ as functions of $Rn_S$')
ax.set_ylabel('$C_{T_S}$, $C_{F_S}$ (-)', fontsize='12')
ax.set_xlabel('$Rn_S$ (-)', fontsize='12')
ax.grid(True)

plt.legend()
plt.savefig('CTSCFS.png', dpi=1400, transparent=True)
#plt.show()


# diagram 4 - SHIP data processing: RT_S as a function of V_S
fig2 = plt.figure(figsize=(8,4))              
ax = fig2.add_subplot(111)
ax.plot(V_S,RT_S/1000, color=(0.8500, 0.3250, 0.0980), label='$R_{T_S}$')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('SHIP data processing: $R_{T_S}$ as a function of $V_S$')
ax.set_ylabel('$R_{T_S}$ (kN)', fontsize='12')
ax.set_xlabel('$V_S$ m/s', fontsize='12')
ax.grid(True)

#plt.legend()
plt.savefig('RTS.png', dpi=1400, transparent=True)
#plt.show()


# diagram 5 - SHIP data processing: P_E as a function of V_S
fig2 = plt.figure(figsize=(8,4))              
ax = fig2.add_subplot(111)
ax.plot(V_S,P_E, color=(0.8500, 0.3250, 0.0980), label='$P_E$')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('SHIP data processing: $R_{T_S}$ as a function of $V_S$')
ax.set_ylabel('$P_E$ (kW)', fontsize='12')
ax.set_xlabel('$V_S$ m/s', fontsize='12')
ax.grid(True)

#plt.legend()
plt.savefig('P_E.png', dpi=1400, transparent=True)
#plt.show()


# 06 - DATA EXPORT
# --------------------------------------------------------------------------------------------------------------
# data export section, making them ready to be shown in a report; LaTeX optimised
auxmatrix_M = np.empty((len(V_M), 6))
auxmatrix_S = np.empty((len(V_S), 7))


# For cycle to create a matrix which contains the most significant data
for i in range(0,len(V_M)):
    auxmatrix_M[i] = np.mat([V_M[i] , Fn_M[i] , Rn_M[i] , CF_M[i] , CT_M[i] , CR_M[i]])
    auxmatrix_S[i] = np.matrix([V_S[i] , Fn_S[i] , Rn_S[i] , CF_S[i] , CT_S[i] , RT_S[i] , P_E[i]])


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


#in case file is missing from the desidered directory
'''# Ottengo il percorso assoluto della directory corrente
current_directory = os.path.abspath(os.path.dirname(__file__))
# Costruisco il percorso completo al file
file_path = os.path.join(current_directory, "eserc1.txt")
print(file_path)'''


































#bin
'''
fig, axs = plt.subplots(2, 4, figsize=(18,7.5))
fig.tight_layout()

# grafico 1 - analisi modello: CF_M in funzione di Rn_M
axs[0,0].plot(CF_M, Rn_M)
axs[0,0].set_title('Analisi modello: $C_{F_M}$ in funzione di $Rn_M$')

# grafico 2 - analisi modello: CT_M in funzione di Rn_M
axs[0, 1].plot(CT_M, Rn_M)
axs[0, 1].set_title('Analisi modello: $C_{T_M}$ in funzione di $Rn_M$')

# grafico 3 - analisi modello: RT_M in funzione di V_M
axs[0, 2].plot(RT_M, V_M)
axs[0, 2].set_title('Analisi modello: $R_{T_M}$ in funzione di $V_M$')

#empty "slot"
axs[0, 3].axis('off')

# grafico 4 - analisi nave: CF_S in funzione di Rn_S
axs[1, 0].plot(CF_S, Rn_S)
axs[1, 0].set_title('Analisi nave: $C_{F_S}$ in funzione di $Rn_S$')

# grafico 5 - analisi nave: CT_S in funzione di Rn_S
axs[1, 1].plot(CT_S, Rn_S)
axs[1, 1].set_title('Analisi nave: $C_{T_S}$ in funzione di $Rn_S$')

# grafico 6 - analisi nave: RT_S in funzione di V_S
axs[1, 2].plot(RT_S, V_S)
axs[1, 2].set_title('Analisi nave: $R_{T_S}$ in funzione di $V_S$')

# grafico 7 - analisi nave: P_E in funzione di V_S
axs[1, 3].plot(P_E, V_S)
axs[1, 3].set_title('Analisi nave: $P_E$ in funzione di $V_S$')
'''