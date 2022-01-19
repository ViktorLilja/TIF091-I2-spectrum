import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

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


#------------------#
#     PLOTTING     #
#------------------#

plot_energy = True

if plot_energy:
    # Probability distributions
    P_X = np.multiply(psi_X, np.conj(psi_X))
    P_B = np.multiply(psi_B, np.conj(psi_B))

    # Plot energy levels for ground state
    plt.plot(r/au.Ang, I2["X"].potential(r)/au.cminv)
    for n in range(0, 10):
        E = E_X[n] + 0.1*P_X[:,n]
        plt.plot(r/au.Ang, E/au.cminv , color="gray")

    # Plot energy levels for excited state
    plt.plot(r/au.Ang, I2["B"].potential(r)/au.cminv)
    for n in range(0, 10):
        E = E_B[n] + 0.1*P_B[:,n]
        plt.plot(r/au.Ang, E/au.cminv, color="gray")

    plt.xlabel("Atomic separation [Å]")
    plt.ylabel("Energy [cm^-1]")
    plt.show()

