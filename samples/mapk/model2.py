#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

import math

model='mapk2'



V_str = sys.argv[1]
D_ratio_str = sys.argv[2]
D_mode = sys.argv[3]
ti_str = sys.argv[4]
mode = sys.argv[5]
seq = sys.argv[6]
T_str = sys.argv[7]

V = float(V_str)
D_ratio = float(D_ratio_str)
ti = float(ti_str)
T = float(T_str)

if ti == 0:
    ki = float('inf')
else:
    ki = math.log(2) / ti


D_ref = 1e-12

D_move = D_ref * D_ratio

if D_mode == 'normal':
    D_react = D_move
elif D_mode == 'fixed':
    D_react = D_ref



#V = 1e-15 #liter

L = math.pow(V * 1e-3, 1.0 / 3.0)

s = EGFRDSimulator()
s.set_world_size(L)

N = 200
matrix_size = min(max(3, int((3 * N) ** (1.0/3.0))), 60)
print 'matrix_size=', matrix_size
s.set_matrix_size(matrix_size)


box1 = CuboidalRegion([0,0,0],[L,L,L])
# not supported yet
#s.add_surface(box1)

radius = 2.5e-9

m = ParticleModel()

K = m.new_species_type('K', D_move, radius)
KK = m.new_species_type('KK', D_move, radius)
P = m.new_species_type('P', D_move, radius)
Kp = m.new_species_type('Kp', D_move, radius)
Kpp = m.new_species_type('Kpp', D_move, radius)
K_KK = m.new_species_type('K_KK', D_move, radius)
Kp_KK = m.new_species_type('Kp_KK', D_move, radius)
Kpp_P = m.new_species_type('Kpp_P', D_move, radius)
Kp_P = m.new_species_type('Kp_P', D_move, radius)

# inactive forms
Kpi = m.new_species_type('Kpi', D_move, radius)
Kppi = m.new_species_type('Kppi', D_move, radius)
Ki = m.new_species_type('Ki', D_move, radius)

s.set_model(m)

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
kD = k_D(D_react * 2, sigma)


s.throw_in_particles(K, C2N(200e-9, V), box1)
s.throw_in_particles(KK, C2N(50e-9, V), box1)
s.throw_in_particles(P, C2N(50e-9, V), box1)

# print kD
# print k_a(Mtom3(0.02e9), kD)
# print k_a(Mtom3(0.032e9), kD)
# sys.exit(0)

#end_time = .5
end_time = 10
while 1:
    s.step()
    next_time = s.scheduler.getTopTime()
    if next_time > end_time:
        s.stop(end_time)
        break

s.reset()

r1 = create_binding_reaction_rule(K, KK, K_KK, k_a(Mtom3(0.02e9), kD))
m.network_rules.add_reaction_rule(r1)
r2 = create_unbinding_reaction_rule(K_KK, K, KK, k_d(1.0, Mtom3(0.02e9), kD))
m.network_rules.add_reaction_rule(r2)
#r3 = create_unbinding_reaction_rule(K_KK, Kp, KK, k_d(1.5, Mtom3(0.02e9), kD))
r3a = create_unbinding_reaction_rule(K_KK, Kpi, KK, 1.5)
m.network_rules.add_reaction_rule(r3a)
r3b = create_unimolecular_reaction_rule(Kpi, Kp, ki)
m.network_rules.add_reaction_rule(r3b)

r4 = create_binding_reaction_rule(Kp, KK, Kp_KK, k_a(Mtom3(0.032e9), kD))
m.network_rules.add_reaction_rule(r4)
r5 = create_unbinding_reaction_rule(Kp_KK, Kp, KK, k_d(1.0, Mtom3(0.032e9), kD))
m.network_rules.add_reaction_rule(r5)
#r6 = create_unbinding_reaction_rule(Kp_KK, Kpp, KK, k_d(15.0, Mtom3(0.032e9), kD))
r6a = create_unbinding_reaction_rule(Kp_KK, Kppi, KK, 15.0)
m.network_rules.add_reaction_rule(r6a)
r6b = create_unimolecular_reaction_rule(Kppi, Kpp, ki)
m.network_rules.add_reaction_rule(r6b)

r7 = create_binding_reaction_rule(Kpp, P, Kpp_P, k_a(Mtom3(0.02e9), kD))
m.network_rules.add_reaction_rule(r7)
r8 = create_unbinding_reaction_rule(Kpp_P, Kpp, P, k_d(1.0, Mtom3(0.02e9), kD))
m.network_rules.add_reaction_rule(r8)
#r9 = create_unbinding_reaction_rule(Kpp_P, Kp, P, k_d(1.5, Mtom3(0.02e9), kD))
r9a = create_unbinding_reaction_rule(Kpp_P, Kpi, P, 1.5)
m.network_rules.add_reaction_rule(r9a)
# same as r3b
#r9b = create_unimolecular_reaction_rule(Kpi, Kp, ki)
#m.network_rules.add_reaction_rule(r9b)

r10 = create_binding_reaction_rule(Kp, P, Kp_P, k_a(Mtom3(0.032e9), kD))
m.network_rules.add_reaction_rule(r10)
r11 = create_unbinding_reaction_rule(Kp_P, Kp, P, k_d(1.0, Mtom3(0.032e9), kD))
m.network_rules.add_reaction_rule(r11)
#r12 = create_unbinding_reaction_rule(Kp_P, K, P, k_d(15.0, Mtom3(0.032e9), kD))
r12a = create_unbinding_reaction_rule(Kp_P, Ki, P, 15.0)
m.network_rules.add_reaction_rule(r12a)
r12b = create_unimolecular_reaction_rule(Ki, K, ki)
m.network_rules.add_reaction_rule(r12b)

#r13 = UnimolecularReactionRule(Kpp, Kp, 1e-1)
#s.add_reaction_rule(r13)
#r14 = UnimolecularReactionRule(Kp, K, 1e-1)
#s.add_reaction_rule(r14)

s.set_model(m)

l = Logger(logname=model + '_' + '_'.join(sys.argv[1:7]))
interrupter = FixedIntervalInterrupter(s, 1e-0, l)

l.start(s)
while s.t < T:
    interrupter.step()
