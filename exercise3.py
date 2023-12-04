"""
Corso di Architetttura Navale
Onde di piccola ampiezza
Last rev.: 13/11/2023
"""

# Libraries
#---------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import cycle
from matplotlib import colormaps
import matplotlib.transforms as mtransforms

import os                                               # to place the workspace in the script's directory
path = os.getcwd()                                      # ""
os.chdir('../')                                         # ""
os.chdir(os.path.dirname(os.path.abspath(__file__)))    # ""


# 01 - Wave characteristics
# --------------------------------------------------------------------------------------------------------------
# wave characteristics as given
Lambda = np.linspace(1,1000,1000)    # wave length, (m)
d = [1, 5, 10, 20, 50, 100, 5000]    # seabed depth, (m)
g = 9.81                             # gravity acceleration, (m/s^2)


#'''
# 02 - Data processing
# --------------------------------------------------------------------------------------------------------------
# Calculating wave frequency (omega_w), wave celerity (c) and wave period (T)
k = 2*np.pi/Lambda        # wave number, (1/m)
omega_w = []              # wave frequency, (rad/s)

# for loop to find wave frequency (\omega) matrix 
for row in range(len(d)):    
    a = []
    for column in range(len(k)):   
        a.append(np.sqrt(g*k[column]*np.tanh(k[column]*d[row])))
    omega_w.append(np.array(a))


c = []                    # wave celerity, (m/s)

# for loop to find wave celerity
for row in range(len(d)):    
    a = []
    for column in range(len(k)):   
        a.append(omega_w[row][column]/k[column])
    c.append(np.array(a))


T = []                    # wave period, (s)

# for loop to find wave period
for row in range(len(d)):    
    a = []
    for column in range(len(k)):   
        a.append(2*np.pi/omega_w[row][column])
    T.append(np.array(a))


#'''
# 03 - Diagrams plotting
# --------------------------------------------------------------------------------------------------------------
# diagram 1 - Wave frequency at different depths
# diagram 2 - Wave celerity at different depths
# diagram 3 - Wave period at different depths

def get_cmap(n, name='tab20b'):             # alternative color maps:
                                            #'tab20c'  
                                            #'jet'
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)
cmap1 = get_cmap(len(d))                    # randomise len(d)-colors, RGB coded'''


# diagram 1 - Wave frequency at different depths
fig2 = plt.figure(figsize=(8,5))              
ax = fig2.add_subplot(111)
for i in range(len(d)):
    ax.plot(Lambda,omega_w[i], color=cmap1(i) , label = f"d = {d[i]}")

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('')
ax.set_ylabel('Wavelength, $\lambda$ (m)', fontsize='12')
ax.set_xlabel('Wave Frequency, $\omega$ (rad/s)', fontsize='12')
ax.set_xscale('log')
ax.grid(True)

plt.legend()
plt.savefig('omega.png', dpi=600, transparent=True)
plt.show()


# diagram 2 - Wave celerity at different depths
fig2 = plt.figure(figsize=(8,5))              
ax = fig2.add_subplot(111)
for i in range(len(d)):
    ax.plot(Lambda,c[i], c=cmap1(i), label = f"d = {d[i]}")
    plt.text(max(Lambda),max(c[i]),f" d = {d[i]}", verticalalignment = 'center')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('')
ax.set_ylabel('Wavelength, $\lambda$ (m)', fontsize='12')
ax.set_xlabel('Wave Celerity, $c$ (m/s)', fontsize='12')
ax.set_xscale('log')
ax.grid(True)

#plt.legend()
plt.savefig('c.png', dpi=600, transparent=True)
plt.show()


# diagram 3 - Wave period at different depths
fig2 = plt.figure(figsize=(8,6))              
ax = fig2.add_subplot(111)
for i in range(len(d)):
    ax.plot(Lambda,T[i], c=cmap1(i), label = f"d = {d[i]}") 
    plt.text(max(Lambda),max(T[i]),f" d = {d[i]}", verticalalignment = 'center')

ax.set_position( [0.15,0.15,0.8,0.8])
#ax.set_title('')
ax.set_ylabel('Wavelength, $\lambda$ (m)', fontsize='12')
ax.set_xlabel('Wave Period, $T$ (s)', fontsize='12')
ax.set_xscale('log')
ax.grid(True)

#plt.legend()
plt.savefig('T.png', dpi=600, transparent=True)
plt.show()
