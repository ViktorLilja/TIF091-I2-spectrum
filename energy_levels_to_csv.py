import module.hartee_units as au
from module.vib_energy_levels import I2
from numpy import linalg as LA
import numpy as np

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


# Save energy levels to textfile
with open('./output/calculated_energy_levels.csv', 'w') as f:

    # Header
    f.write("Electronic state, nu, energy [cm^-1]\n")

    for vbis in range(0,50):
        energy = E_X[vbis] / au.cminv
        f.write("X, " + str(vbis) + ", " + str(energy) + "\n")

    for vprime in range(0,50):
        energy = E_B[vprime] /au.cminv
        f.write("B, " + str(vprime) + ", " + str(energy) + "\n")
