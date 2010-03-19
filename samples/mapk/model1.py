#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

import math



V_str = sys.argv[1]
D_ratio_str = sys.argv[2]
D_mode = sys.argv[3]
mode = sys.argv[4]
seq = sys.argv[5]
T_str = sys.argv[6]

V = float(V_str)
D_ratio = float(D_ratio_str)
T = float(T_str)

D_ref = 1e-12

D_move = D_ref * D_ratio

if D_mode == 'normal':
    D_react = D_move
elif D_mode == 'fixed':
    D_react = D_ref



#V = 1e-15 #liter

L = math.pow(V * 1e-3, 1.0 / 3.0)

s = EGFRDSimulator()
s.setWorldSize(L)

N = 200
matrixSize = min(max(3, int((3 * N) ** (1.0/3.0))), 60)
print 'matrixSize=', matrixSize
s.setMatrixSize(matrixSize)


box1 = CuboidalRegion([0,0,0],[L,L,L])
# not supported yet
#s.addSurface(box1)

model='mapk1'

#D = 2e-12 # run1
D = 1e-12 # run2
#D = 5e-13 # run3
#D = 0.25e-12 # run4

radius = 2.5e-9

m = ParticleModel()

K = m.new_species_type('K', D, radius)
KK = m.new_species_type('KK', D, radius)
P = m.new_species_type('P', D, radius)
Kp = m.new_species_type('Kp', D, radius)
Kpp = m.new_species_type('Kpp', D, radius)
K_KK = m.new_species_type('K_KK', D, radius)
Kp_KK = m.new_species_type('Kp_KK', D, radius)
Kpp_P = m.new_species_type('Kpp_P', D, radius)
Kp_P = m.new_species_type('Kp_P', D, radius)

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

s.setModel(m)

sigma = radius * 2
kD = k_D(D_react * 2, sigma)


s.throwInParticles(K, C2N(200e-9, V), box1)
s.throwInParticles(KK, C2N(50e-9, V), box1)
s.throwInParticles(P, C2N(50e-9, V), box1)

#endTime = .5
endTime = 0
while 1:
    s.step()
    nextTime = s.scheduler.getTopTime()
    if nextTime > endTime:
        s.stop(endTime)
        break

s.reset()

r1 = createBindingReactionRule(K, KK, K_KK, k_a(Mtom3(0.02e9), kD))
m.network_rules.add_reaction_rule(r1)
r2 = createUnbindingReactionRule(K_KK, K, KK, k_d(1.0, Mtom3(0.02e9), kD))
m.network_rules.add_reaction_rule(r2)
#r3 = createUnbindingReactionRule(K_KK, Kp, KK, k_d(1.5, Mtom3(0.02e9), kD))
r3 = createUnbindingReactionRule(K_KK, Kp, KK, 1.5)
m.network_rules.add_reaction_rule(r3)

r4 = createBindingReactionRule(Kp, KK, Kp_KK, k_a(Mtom3(0.032e9), kD))
m.network_rules.add_reaction_rule(r4)
r5 = createUnbindingReactionRule(Kp_KK, Kp, KK, k_d(1.0, Mtom3(0.032e9), kD))
m.network_rules.add_reaction_rule(r5)
#r6 = createUnbindingReactionRule(Kp_KK, Kpp, KK, k_d(15.0, Mtom3(0.032e9), kD))
r6 = createUnbindingReactionRule(Kp_KK, Kpp, KK, 15.0)
m.network_rules.add_reaction_rule(r6)

r7 = createBindingReactionRule(Kpp, P, Kpp_P, k_a(Mtom3(0.02e9), kD))
m.network_rules.add_reaction_rule(r7)
r8 = createUnbindingReactionRule(Kpp_P, Kpp, P, k_d(1.0, Mtom3(0.02e9), kD))
m.network_rules.add_reaction_rule(r8)
#r9 = createUnbindingReactionRule(Kpp_P, Kp, P, k_d(1.5, Mtom3(0.02e9), kD))
r9 = createUnbindingReactionRule(Kpp_P, Kp, P, 1.5)
m.network_rules.add_reaction_rule(r9)

r10 = createBindingReactionRule(Kp, P, Kp_P, k_a(Mtom3(0.032e9), kD))
m.network_rules.add_reaction_rule(r10)
r11 = createUnbindingReactionRule(Kp_P, Kp, P, k_d(1.0, Mtom3(0.032e9), kD))
m.network_rules.add_reaction_rule(r11)
#r12 = createUnbindingReactionRule(Kp_P, K, P, k_d(15.0, Mtom3(0.032e9), kD))
r12 = createUnbindingReactionRule(Kp_P, K, P, 15.0)
m.network_rules.add_reaction_rule(r12)

#r13 = createUnimolecularReactionRule(Kpp, Kp, 1e-1)
#m.network_rules.add_reaction_rule(r13)
#r14 = createUnimolecularReactionRule(Kp, K, 1e-1)
#m.network_rules.add_reaction_rule(r14)

s.setModel(m)


l = Logger(s, 
           logname = model + '_' + '_'.join(sys.argv[1:6]))



#l.setParticleOutput(('Ea','X','EaX','Xp','Xpp','EaI'))
l.setParticleOutInterval(1e-0)
l.log()

while s.t < 30:
    s.step()
    #s.dumpPopulation()
    #s.check()
    l.log()

