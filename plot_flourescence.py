import numpy as np
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


#------------------#
#     PLOTTING     #
#------------------#

# Wavelength of laser: 612nm = 16339.87 cminv
# Roughly equivalent to energy difference between v´´= 0 and v´= 5

# Generate list of lines in emission spectrum
# line = (wavelength,intensity)

vprim = 28   # Select which level is excited by laser

wavelength = []
intensity = []
for vbis in range(0, 12):
    deltaE = I2["B"].energy(vprim) - I2["X"].energy(vbis)
    deltaE /= au.cminv                  # Convert from hartee to cminv
    wavelength.append(1e7 / deltaE)     # Convert from cminv to nm

    # Calculate relatve intensity using Franck-Condon factor
    FCfactor = np.sum(np.multiply(psi_X[:,vbis], psi_B[:,vprim])*dr)
    intensity.append(FCfactor * FCfactor)

plt.scatter(wavelength, intensity)
plt.show()