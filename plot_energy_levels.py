import numpy as np
import math
import matplotlib.pyplot as plt
from numpy import linalg as LA

import constants.hartee_units as au
from constants.vib_energy_levels import I2

#---------------------------------------------#
#     CALCULATION OF ENERGY WAVEFUNCTIONS     #
#---------------------------------------------#

# Number of points to sample wavefuntion
N = 1000

# Coordinate system
interval = np.array([2.3, 5.8]) * au.Ang
r = np.linspace(interval[0],interval[1],N)
dr = (interval[1]-interval[0])/N

# get vibrational energy Hamiltonian in internuclear distance coordinates
X_H = I2["X"].hamiltonian(r)
B_H = I2["B"].hamiltonian(r)

# Find energy eigenstates
# eigh: special eigenvalue solver for hermetian matrices
E_X, psi_X = LA.eigh(X_H)
E_B, psi_B = LA.eigh(B_H)


#---------------------------------------------#
#     COMPARE CALCULATED ENERGY TO THEORY     #
#---------------------------------------------#

N = 50
v = np.array(range(0,N))

# Plotting
plt.plot(v, E_X[0:N] / au.cminv, "rx")
plt.plot(v, I2["X"].energy(v) / au.cminv, "k--")
plt.show()