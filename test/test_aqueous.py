

import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from aqueous import gibbs_energy, AqueousSolution
import numpy as np

def test_waterdissociation():
    species = ['H+','OH-']
    m = np.array([1.0e-12,1.0e-12])
    s = AqueousSolution(species)
    T = 298
    P = 1
    s.equilibrate(m, T, P)
    print("Equilibrium concentrations of H+ and OH- in water:")
    for i,sp in enumerate(species):
        print("{:10}".format(sp)+"= "+"{:15}".format('%.10e'%m[i]))
    print()    
    
def test_enceladus():
    species = ['H2,aq','CO2,aq','CH4']
    m = np.array([1e-4, 7e-5, 3e-5]) # Table S11 in Waite et al. (2017)
    s = AqueousSolution(species)
    T = 273.15
    P = 1
    s.equilibrate(m, T, P)
    A = -(s.G_opt-s.G_init)*(1/3e-5)
    print("Chemical Affinity of 4H2 + CO2 => CH4 + 2H2O in Enceladus' Ocean = "+("%.1f"%A)+" J/(mol CH4)")
    print()

if __name__ == "__main__":
    test_waterdissociation()
    test_enceladus()