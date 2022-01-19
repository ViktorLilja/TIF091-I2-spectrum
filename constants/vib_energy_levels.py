# Vibrational energy levels for diatomic molecule
# All parameters muste be given in Hartee atomic units

# mA:   mass of atom 1
# mB:   mass of atom 2
# we:   first order energy term
# wexe: second order energy term
# Te:   Potential well bottom relative to ground state bottom
# re:   internuclear distance with minimum potential energy

import numpy as np
import math
import constants.hartee_units as au

#------------------------#
#     VibLevel class     #
#------------------------#

class VibLevel:
    def __init__(self, mA, mB, we, wexe, Te, re):

        # Initialize parameters converted to hartee units
        self.mu     = mA*mB/(mA+mB) # Reduced mass
        self.we     = we            
        self.wexe   = wexe          
        self.Te     = Te            
        self.re     = re            

        # Parameters for Morse potential
        # V(r) = De [1 - exp(-beta * (r - re))]
        self.De     = self.we**2 / (4 * self.wexe)   
        self.beta   = we * math.sqrt(0.5 * self.mu / self.De)
        
        # Source of formulas:   "Molecular spectrum of Iodine Revisited"
        #                       Ian J. McNaught (1980)



    # Morse potential energy at r
    def potential(self, r):
        # Definition of Morse potential
        return self.Te + self.De * (1 - np.exp(-self.beta*(r-self.re)))**2



    # Hamiltonian matrix in position basis
    # r is a vector of equally spaced position values
    def hamiltonian(self, r):
        N = len(r)
        interval = r[-1]-r[0]
        dr = interval/N

        H = np.zeros((N,N))
        i,j = np.indices(H.shape)

        # Operator matrix for -1 * numerical second derivative with respect to x
        #   formula: d2y/dx2 = (y(k-1) - 2y(k) + y(k+1))/4dx2
        H[i==j] = 2/(dr*dr)
        H[abs(i-j)==1] = -1/(dr*dr)

        # Scale by hbar^2/2m, but hbar = 1 in atomic units
        # After this, H is the kinetic energy Hamiltonian
        H *= 0.5 / self.mu

        # Add potential to diagonal to form full H = T + V matrix
        H[i==j] += self.potential(r)

        return H



    # Energy of vibrational energy level with quantum number nu
    def energy(self, nu):
        E1 = self.we * (0.5 + nu)           # First order term
        E2 = self.wexe * (0.5 + nu)**2      # Second order term
        return self.Te - E1 + E2


#--------------------------------------------------------------#
#     Some vibrational energy levels of diatomic molecules     #
#--------------------------------------------------------------#

# Iodine
# Source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7553562&Mask=1000#Diatomic
I2 = {
    "X": VibLevel(  126.90447 * au.u, 
                    126.90447 * au.u, 
                    214.50    * au.cminv, 
                    0.614     * au.cminv, 
                    0         * au.cminv, 
                    2.666     * au.Ang),

    "B": VibLevel(  126.90447 * au.u, 
                    126.90447 * au.u, 
                    125.69    * au.cminv, 
                    0.764     * au.cminv, 
                    15769.01  * au.cminv, 
                    3.024     * au.Ang)
}