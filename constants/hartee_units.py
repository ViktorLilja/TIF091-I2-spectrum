# Physical constants and conversion factors in Hartee atomic units
# https://en.wikipedia.org/wiki/Hartree_atomic_units

# Definitions of Hartee units
hbar    = 1                                 # Reduced Planck constant
me      = 1                                 # Electron rest mass
b       = 1                                 # Bohr radius
e       = 1                                 # Elementary charge

# Derived units
Eh      = 1                                 # Hartee energy unit
K       = 1                                 # Kelvin

# Conversion from SI
kg      = 1.0977683828808e3 * me            # kilogram
m       = 1.8897259885789e9 * b             # meter
J       = 2.293710448690e17 * Eh            # Joule
s       = 4.1341373e16      * hbar/Eh       # Second


# Conversion from various units used in molecular physics and chemistry
u       = 1822.888          * me            # atomic mass unit
Ang     = 1.8897259886      * b             # Angstrom
cminv   = 4.55633e-6        * Eh            # inverse cm as a unit of energy
eV      = 0.0367493088      * Eh            # electron Volt

# Physical constants
c       = 137.036           * b*Eh/hbar     # speed of light
k       = 1.380658e-23      * J/K           # Bolzmanns constant