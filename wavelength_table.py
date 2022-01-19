import constants.hartee_units as au
from constants.vib_energy_levels import I2
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

vbis_max = 6
vprime_max = 30

# Write table of wavelengths to textfile in Latex format
with open('output/wavelength_table.txt', 'w') as f:

    # Header
    for vbis in range(0, vbis_max+1): f.write(" & " + str(vbis))
    f.write("\\\\ \\specialrule{0.8pt}{0pt}{0pt} \n")

    for vprime in range(0,vprime_max):

        f.write("%d & " % vprime)

        for vbis in range(0,vbis_max+1):
            energy_upper = I2["B"].energy(vprime) / au.cminv
            energy_lower = I2["X"].energy(vbis)   / au.cminv
            energy_diff = energy_upper - energy_lower

            # lambda [nm] = 10^7 / wavenumber [cm^-1]
            wavelength = 1e7 / energy_diff

            f.write(("%.2f" % wavelength).replace('.', ','))

            if vbis == vbis_max:    f.write(" \\\\")
            else:                   f.write(" & ")

        f.write("\n")