import numpy as np
import math
import matplotlib.pyplot as plt
from numpy import linalg as LA

import module.hartee_units as au
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

# Select range of transisitons to simulate
vprime_range = range(0,500)
vbis_range   = range(0,50)

wavelengths = []
intensities = []
for vbis in vbis_range:
    for vprime in vprime_range:
        # Calculate wavelength of light from energy difference between levels
        E_excited = E_B[vprime]
        E_ground  = E_X[vbis]
        E_0       = E_X[0]
        deltaE = E_excited - E_ground
        deltaE /= au.cminv                  # Convert from hartee to cminv
        wavelength = 1e7 / deltaE           # Convert from cminv to nm
        wavelengths.append(wavelength)    
        
        # Calculate Franck-Condon factor
        intensity = np.sum(np.multiply(psi_X[:,vbis], psi_B[:,vprime])*dr)
        
        # Scale by probablility of hot band being excited
        # Assuming bolzmann distribution at temperature 100C
        # T = (100 + 273.15) * au.K
        # beta = 1 / (au.k * T)
        # intensity *= math.exp(- beta * (E_ground-E_0))

        frequency = 3 / (wavelength * 10)
        intensities.append(intensity * intensity * frequency**3)


#---------------------------------------#
#     Generate spectrum from peaks      #
#---------------------------------------#

n_points = 10000          # Number of datapoints in interval
wl_inter = (500, 750)    # [nm] wavelength scan interval

# Gaussian function as shape of peak
def gaussian(x, mean, std, height):
    normalization = 1/ (math.sqrt(2 * math.pi) * std)
    norm_dist = normalization * np.exp(-0.5*((lam-mean)/std)**2)
    return height * norm_dist

# Coordinate system
lam = np.linspace(wl_inter[0], wl_inter[1], n_points)   # [nm] wavelength
spec = np.zeros_like(lam)                               # absorbance

# Iterate over peaks and add them to the spectrum
for i in range(0, len(wavelengths)):
    peak_mid = wavelengths[i]
    peak_height = intensities[i]
    peak_width = 0.6

    peak = gaussian(wavelengths[i], peak_mid, peak_width, peak_height)
    spec += peak

# Normalize absorbance to highest peak
spec /= np.max(spec)
p = 2.0*spec
T = np.exp(-p)
A = 1-T

plt.plot(lam, A)
plt.xlim(520, 660)
plt.show()

#----------------------#
#     Save to CSV      #
#----------------------#

with open("./output/simulated_a_spec.csv", "w") as f:
    # Header
    f.write("Wavelength [nm], Absorption\n")

    # Data
    for i in range(0, n_points):
        f.write("%.5f, %.8f\n" % (lam[i], A[i]))
