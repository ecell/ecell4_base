from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions

@species_attributes
def attributegen():
    A(f=off) | 1
    B() | 2
    C(f=off) | 0
    D() | 3
    E(f=off) | 0
    I() | 0

@reaction_rules
def rulegen():
    A(f=_1) + B() == C(f=_1) | (1,2)
    C(f=_1) + D() == E(f=_1) | (3,4)
    A(f=off) + I == A(f=on) | (5,6)
 
#begin parameters
#NA 6.02e23 # Acogadro's number( molecules/mol)
#f 0.1 # Fraction of the cell to simulate
#Vo f*1.0e-10 # Extracellular volume=1/cell_density (L)
#V f*3.0e-12 # Cytoplasmic volume (L)
## Initial concentrations (copies per cell)
#A_tot 10000
#B_tot 8000
#D_tot 50000
## Rate constants
## Divide by NA*V to convert bimolecular rate constants
## from /M/sec to /(molecule/cell)/sec
#kpAB 3.0e6/(NA*V)
#kmAB 0.06
#kpCD 1.0e6/(NA*V)
#kmCD 0.06
#kpI 1.0e7/(NA*V)
#kmI 0.1
#end parameters
#
#begin molecule types
#A(f~off~on)
#B()
#C(f~off~on)
#D()
#E(f~off~on)
#I()
#end molecule types
#
#begin seed species
#A(f~off) A_tot
#B() B_tot
#C(f~off) 0
#D() D_tot
#E(f~off) 0
#I() 0
#end seed species
#
#begin reaction rules
#A(f%1) + B() <-> C(f%1) kpAB, kmAB
#C(f%1) + D() <-> E(f%1) kpCD, kmCD
#A(f~off) + I <-> A(f~on) kpI, kmI
#end reaction rules
#
#generate_network({overwrite=>1});
