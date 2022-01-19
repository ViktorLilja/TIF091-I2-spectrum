import constants.hartee_units as au
from constants.vib_energy_levels import I2

# Save energy levels to textfile
with open('output/calculated_energy_levels.csv', 'w') as f:

    # Header
    f.write("Electronic state, nu, energy [cm^-1]\n")

    for vbis in range(0,50):
        energy = I2["X"].energy(vbis)/au.cminv
        f.write("X, " + str(vbis) + ", " + str(energy) + "\n")

    for vprime in range(0,50):
        energy = I2["B"].energy(vprime)/au.cminv
        f.write("B, " + str(vprime) + ", " + str(energy) + "\n")
