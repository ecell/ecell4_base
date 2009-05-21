#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

import math


# PROCESSIVE MODEL

model='mapk6'

V_str = sys.argv[1]
D_ratio_str = sys.argv[2]
N_KK_str = sys.argv[3]
N_P_str = sys.argv[4]
N_K_total_str = sys.argv[5]
ti_str = sys.argv[6]  # DUMMY
T_str = sys.argv[7]

V = float( V_str )
D_ratio = float( D_ratio_str )
ti = float( ti_str )

N_KK = int( N_KK_str )
N_P = int( N_P_str )
N_K_total = int( N_K_total_str )

# N_K == N_Kpp
N_K = N_K_total * .5
N_Kpp = N_K

T = float( T_str )

if ti == 0:
     ki = float( 'inf' )
else:
     ki = math.log( 2 ) / ti


D_ref = 1e-12

D_move = D_ref * D_ratio

D_mode = 'fixed'

if D_mode == 'normal':
    D_react = D_move
elif D_mode == 'fixed':
    D_react = D_ref
else:
    raise 'unexpected value'

# V in liter, L in meter
L = math.pow( V * 1e-3, 1.0 / 3.0 )

s = EGFRDSimulator()
s.setWorldSize( L )

N = 300
matrixSize = min( max( 3, int( (3 * N) ** (1.0/3.0) ) ), 60 )
print 'matrixSize=', matrixSize
s.setMatrixSize( matrixSize )


box1 = CuboidalSurface( [0,0,0],[L,L,L] )
# not supported yet
#s.addSurface( box1 )

radius = 2.5e-9

K = Species( 'K', D_move, radius )
s.addSpecies( K )
KK = Species( 'KK', D_move, radius )
s.addSpecies( KK )
P = Species( 'P', D_move, radius )
s.addSpecies( P )
Kp = Species( 'Kp', D_move, radius )
s.addSpecies( Kp )
Kpp = Species( 'Kpp', D_move, radius )
s.addSpecies( Kpp )
K_KK = Species( 'K_KK', D_move, radius )
s.addSpecies( K_KK )
Kp_KK = Species( 'Kp_KK', D_move, radius )
s.addSpecies( Kp_KK )
Kpp_KK = Species( 'Kpp_KK', D_move, radius )
s.addSpecies( Kpp_KK )
Kpp_P = Species( 'Kpp_P', D_move, radius )
s.addSpecies( Kpp_P )
Kp_P = Species( 'Kp_P', D_move, radius )
s.addSpecies( Kp_P )

# inactive forms
KKi = Species( 'KKi', D_move, radius )
s.addSpecies( KKi )
Pi = Species( 'Pi', D_move, radius )
s.addSpecies( Pi )



#  1 2   K + KK   <-> K_KK
#  3     K_KK       -> Kp_KK
#  4     Kp_KK       -> Kpp + KK
#  5 6   Kpp + P <-> Kpp_P
#  7     Kpp_P     -> Kp_P
#  8     Kp_P  <-> K + P


sigma = radius * 2
kD = k_D( D_react * 2, sigma )

s.throwInParticles( Kpp, N_Kpp, box1 )
s.throwInParticles( K, N_K, box1 )
s.throwInParticles( KK, N_KK, box1 )
s.throwInParticles( P, N_P, box1 )

# print kD
# print k_a( Mtom3( 0.02e9 ), kD )
# print k_a( Mtom3( 0.032e9 ), kD )
# sys.exit(0)

endTime = 5
while 1:
    s.step()
    nextTime = s.scheduler.getTopTime()
    if nextTime > endTime:
        s.stop( endTime )
        break

s.reset()
k1 = k_a( Mtom3( 0.02e9 ), kD )
k2 = k_d( 1.0, Mtom3( 0.02e9 ), kD )
k3 = 1.5
k4 = k_a( Mtom3( 0.032e9 ), kD )
k5 = k_d( 1.0, Mtom3( 0.032e9 ), kD )
k6 = 15.0

r1 = BindingReactionRule( K, KK, K_KK, k1 )
s.addReactionRule( r1 )
r2 = UnbindingReactionRule( K_KK, K, KK, k2 )
s.addReactionRule( r2 )
r3 = UnimolecularReactionRule( K_KK, Kp_KK, k3 )
s.addReactionRule( r3 )
r4 = UnbindingReactionRule( Kp_KK, Kpp, KK, k6 )
s.addReactionRule( r4 )


r5 = BindingReactionRule( Kpp, P, Kpp_P, k1 )
s.addReactionRule( r5 )
r6 = UnbindingReactionRule( Kpp_P, Kpp, P, k2 )
s.addReactionRule( r6 )
r7 = UnimolecularReactionRule( Kpp_P, Kp_P,k3 )
s.addReactionRule( r7 )
r8 = UnbindingReactionRule( Kp_P, K, P, k6 )
s.addReactionRule( r8 )


logname = model + '_' + '_'.join( sys.argv[1:7] )  + '_' +\
          os.environ[ 'SGE_TASK_ID' ]
l = Logger( s, 
            logname = logname,
            comment = '@ model=\'%s\'; D_move=%g; D_react=%g\n' %
            ( model, D_move, D_react ) +
            '#@ V=%s; N_K_total=%d; N_K=%d; N_Kpp=%d; N_KK=%d; N_P=%d;\n' % 
            ( V_str, N_K_total, N_K, N_Kpp, N_KK, N_P ) +
            '#@ k1=%g; k2=%g; k3=%g; k4=%g; k5=%g; k6=%g;\n' %
            ( k1, k2, k3, k4, k5, k6 ) +
            '#@ ti=%g; ki=%g;' %
            ( ti, ki ) )

rfile = open( 'data/' + logname + '_reactions.dat', 'w' )


#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-0 )
l.log()

while s.t < T:
    s.step()

    if s.lastReaction:
        r = s.lastReaction
        line = '( %18.18g,\t%s,\t%s )\n' % ( s.t, r.reactants, r.products )
        #print line
        rfile.write( line )
        rfile.flush()

        l.log()

