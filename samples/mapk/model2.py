#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

import math



V_str = sys.argv[1]
D_ratio_str = sys.argv[2]
D_mode = sys.argv[3]
ti_str = sys.argv[4]
mode = sys.argv[5]
seq = sys.argv[6]
T_str = sys.argv[7]

V = float( V_str )
D_ratio = float( D_ratio_str )
ti = float( ti_str )
T = float( T_str )

if ti == 0:
    ki = float( 'inf' )
else:
    ki = math.log( 2 ) / ti


D_ref = 1e-12

D_move = D_ref * D_ratio

if D_mode == 'normal':
    D_react = D_move
elif D_mode == 'fixed':
    D_react = D_ref



#V = 1e-15 #liter

L = math.pow( V * 1e-3, 1.0 / 3.0 )

s = EGFRDSimulator()
s.setWorldSize( L )

N = 200
matrixSize = min( max( 3, int( (3 * N) ** (1.0/3.0) ) ), 60 )
print 'matrixSize=', matrixSize
s.setMatrixSize( matrixSize )


box1 = CuboidalSurface( [0,0,0],[L,L,L] )
# not supported yet
#s.addSurface( box1 )

model='mapk1'

#D = 2e-12 # run1
D = 1e-12 # run2
#D = 5e-13 # run3
#D = 0.25e-12 # run4

radius = 5e-9

K = Species( 'K', D, radius )
s.addSpecies( K )
KK = Species( 'KK', D, radius )
s.addSpecies( KK )
P = Species( 'P', D, radius )
s.addSpecies( P )
Kp = Species( 'Kp', D, radius )
s.addSpecies( Kp )
Kpp = Species( 'Kpp', D, radius )
s.addSpecies( Kpp )
K_KK = Species( 'K_KK', D, radius )
s.addSpecies( K_KK )
Kp_KK = Species( 'Kp_KK', D, radius )
s.addSpecies( Kp_KK )
Kpp_KK = Species( 'Kpp_KK', D, radius )
s.addSpecies( Kpp_KK )
Kpp_P = Species( 'Kpp_P', D, radius )
s.addSpecies( Kpp_P )
Kp_P = Species( 'Kp_P', D, radius )
s.addSpecies( Kp_P )

# inactive forms
Kpi = Species( 'Kpi', D, radius )
s.addSpecies( Kpi )
Kppi = Species( 'Kppi', D, radius )
s.addSpecies( Kppi )
Ki = Species( 'Ki', D, radius )
s.addSpecies( Ki )



#  1 2   K + KK   <-> K_KK
#  3     K_KK       -> Kp + KK
#  4 5   Kp + KK  <-> Kp_KK
#  6     Kp_KK      -> Kpp + KK 
#  7 8   Kpp + P <-> Kpp_P
#  9     Kpp_P     -> Kp + P
# 10 11  Kp + P  <-> Kp_P
# 12     Kp_P      -> K + P
# 12     Kp_P      -> K + P
# 13     Kpp     -> Kp
# 14     Kp      -> K


sigma = radius * 2
kD = k_D( D_react * 2, sigma )


s.throwInParticles( K, C2N( 200e-9, V ), box1 )
s.throwInParticles( KK, C2N( 50e-9, V ), box1 )
s.throwInParticles( P, C2N( 50e-9, V ), box1 )

#endTime = .5
endTime = 0
while 1:
    s.step()
    nextTime = s.scheduler.getTopTime()
    if nextTime > endTime:
        s.stop( endTime )
        break

s.reset()

r1 = BindingReactionType( K, KK, K_KK, k_a( Mtom3( 0.02e9 ), kD ) )
s.addReactionType( r1 )
r2 = UnbindingReactionType( K_KK, K, KK, k_d( 1.0, Mtom3( 0.02e9 ), kD ) )
s.addReactionType( r2 )
#r3 = UnbindingReactionType( K_KK, Kp, KK, k_d( 1.5, Mtom3( 0.02e9 ), kD ) )
r3a = UnbindingReactionType( K_KK, Kpi, KK, 1.5 )
s.addReactionType( r3a )
r3b = UnimolecularReactionType( Kpi, Kp, ki )
s.addReactionType( r3b )

r4 = BindingReactionType( Kp, KK, Kp_KK, k_a( Mtom3( 0.032e9 ), kD ) )
s.addReactionType( r4 )
r5 = UnbindingReactionType( Kp_KK, Kp, KK, k_d( 1.0, Mtom3( 0.032e9 ), kD ) )
s.addReactionType( r5 )
#r6 = UnbindingReactionType( Kp_KK, Kpp, KK, k_d( 15.0, Mtom3( 0.032e9 ), kD ) )
r6a = UnbindingReactionType( Kp_KK, Kppi, KK, 15.0 )
s.addReactionType( r6a )
r6b = UnimolecularReactionType( Kppi, Kpp, ki )
s.addReactionType( r6b )

r7 = BindingReactionType( Kpp, P, Kpp_P, k_a( Mtom3( 0.02e9 ), kD ) )
s.addReactionType( r7 )
r8 = UnbindingReactionType( Kpp_P, Kpp, P, k_d( 1.0, Mtom3( 0.02e9 ), kD ) )
s.addReactionType( r8 )
#r9 = UnbindingReactionType( Kpp_P, Kp, P, k_d( 1.5, Mtom3( 0.02e9 ), kD ) )
r9a = UnbindingReactionType( Kpp_P, Kpi, P, 1.5 )
s.addReactionType( r9a )
r9b = UnimolecularReactionType( Kpi, Kp, ki )
s.addReactionType( r9b )

r10 = BindingReactionType( Kp, P, Kp_P, k_a( Mtom3( 0.032e9 ), kD ))
s.addReactionType( r10 )
r11 = UnbindingReactionType( Kp_P, Kp, P, k_d( 1.0, Mtom3( 0.032e9 ), kD ) )
s.addReactionType( r11 )
#r12 = UnbindingReactionType( Kp_P, K, P, k_d( 15.0, Mtom3( 0.032e9 ), kD ) )
r12a = UnbindingReactionType( Kp_P, Ki, P, 15.0 )
s.addReactionType( r12a )
r12b = UnimolecularReactionType( Ki, K, ki )
s.addReactionType( r12b )

#r13 = UnimolecularReactionType( Kpp, Kp, 1e-1 )
#s.addReactionType( r13 )
#r14 = UnimolecularReactionType( Kp, K, 1e-1 )
#s.addReactionType( r14 )


l = Logger( s, 
            logname = model + '_' + '_'.join( sys.argv[1:6] ) )



#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-0 )
l.log()

while s.t < 30:
    s.step()
    #s.dumpPopulation()
    #s.check()
    l.log()

