import numpy as np
import math
import matplotlib.pyplot as plt
from numpy import linalg as LA

from module import hartee_units as au
from module.vib_energy_levels import I2

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


#---------------------------------------------------#
#     Calculate position and intensity of peaks     #
#---------------------------------------------------#

wavelengths = []
intensities = []

vprimes = [28]
weights = [1]

# Iterate over excited states
for i in range(0, len(vprimes)):
    vprime = vprimes[i]
    weight = weights[i]

    # Iterate over ground states
    for vbis in range(0, 25):
        # Calculate wavelength of light from energy difference between levels
        deltaE = E_B[vprime] - E_X[vbis]
        deltaE /= au.cminv                  # Convert from hartee to cminv
        wavelength = 1e7 / deltaE           # Convert from cminv to nm
        wavelengths.append(wavelength)    

        # Calculate intensity using Franck-Condon factor
        FCfactor = np.sum(np.multiply(psi_X[:,vbis], psi_B[:,vprime])*dr)
        frequency = 3 / (wavelength * 10)   # Frequency from wavelegth in nm
        intensities.append(weight * FCfactor * FCfactor * frequency**3)

    print(wavelengths)

#---------------------------------------#
#     Generate spectrum from peaks      #
#---------------------------------------#

# Coordinate system
n_points = 1000          # Number of datapoints in interval
wl_inter = (530, 630)    # [nm] wavelength scan interval
lam = np.linspace(wl_inter[0], wl_inter[1], n_points)   # [nm] wavelength
spec = np.zeros_like(lam)                               # intensity (FC factor)

# Gaussian function as shape of peak
def gaussian(x, mean, std, height):
    normalization = 1/ (math.sqrt(2 * math.pi) * std)
    norm_dist = normalization * np.exp(-0.5*((lam-mean)/std)**2)
    return height * norm_dist

# Iterate over peaks and add them to the spectrum
for i in range(0, len(wavelengths)):
    peak_mid = wavelengths[i]
    peak_height = intensities[i]
    peak_width = 0.5

    peak = gaussian(wavelengths[i], peak_mid, peak_width, peak_height)
    spec += peak

# Normalize intensity to highest peak
spec /= np.max(spec)

plt.plot(lam, spec)
plt.show()

#----------------------#
#     Save to CSV      #
#----------------------#

with open("./output/simulated_f_spec.csv", "w") as f:
    # Header
    f.write("wavelength [nm], relative intensity\n")

    # Data
    for i in range(0, n_points):
        f.write("%.5f, %.8f\n" % (lam[i], spec[i]))
