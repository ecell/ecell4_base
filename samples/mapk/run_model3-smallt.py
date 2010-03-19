#!/usr/bin/env python

'''
# D_factor ti T N

LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py 1 1 100 10

'''


from egfrd import *
from bd import *

def run(outfilename, D_factor, ti, T, N):
    print outfilename

    outfile_t = open(outfilename + '_t.dat', 'w')

    outfile_t.write('%d\n' % N)

    for i in range(N):
        t = singlerun(D_factor, ti, T)
        #print i, t

        if t != -1:
            outfile_t.write('%g\n' % t)
            outfile_t.flush()
 

        #[outfile_r.flush() for outfile_r in outfile_r_list]


    outfile_t.close()
    #[outfile_r.close() for outfile_r in outfile_r_list]



def singlerun(D_factor, ti, T):

    V = 1e-15
    #V = 1e-16
    D_ratio = 1
    
    if ti == 0:
        ki = float('inf')
    else:
        ki = math.log(2) / ti


    D_ref = 1e-12

    D_move = D_ref * D_factor

    D_react = D_ref

    # V in liter, L in meter
    L = math.pow(V * 1e-3, 1.0 / 3.0)

    s = EGFRDSimulator()
    s.set_world_size(L)

    N = 180
    matrix_size = min(max(3, int((3 * N) ** (1.0/3.0))), 60)
    s.set_matrix_size(matrix_size)


    box1 = CuboidalRegion([0,0,0],[L,L,L])

    radius = 2.5e-9

    m = ParticleModel()

#     K = m.new_species_type('K', D_move, radius)
    KK = m.new_species_type('KK', D_move, radius)
#     P = m.new_species_type('P', D_move, radius)
    Kp = m.new_species_type('Kp', D_move, radius)
#     Kpp = m.new_species_type('Kpp', D_move, radius)
#     K_KK = m.new_species_type('K_KK', D_move, radius)
    Kp_KK = m.new_species_type('Kp_KK', D_move, radius)
#     Kpp_P = m.new_species_type('Kpp_P', D_move, radius)
#     Kp_P = m.new_species_type('Kp_P', D_move, radius)

    # inactive forms
    KKi = m.new_species_type('KKi', D_move, radius)
#     Pi = m.new_species_type('Pi', D_move, radius)

    s.set_model(m)

    #  1 2   K + KK   <-> K_KK
    #  3     K_KK       -> Kp + KKi
    #  4 5   Kp + KK  <-> Kp_KK
    #  6     Kp_KK      -> Kpp + KKi 
    #  7 8   Kpp + P <-> Kpp_P
    #  9     Kpp_P     -> Kp + Pi
    # 10 11  Kp + P  <-> Kp_P
    # 12     Kp_P      -> K + Pi
    # 13     KKi     -> KK
    # 14     Pi      -> P


    sigma = radius * 2
    kD = k_D(D_react * 2, sigma)

    N_K = C2N(200e-9, V) 
    N_KK = C2N(50e-9, V)
    N_P = C2N(50e-9, V)

    #print N_KK
    #s.throw_in_particles(K, N_K, box1)
    #s.throw_in_particles(KK, N_KK, box1)
    #s.throw_in_particles(P, N_P, box1)
    
    s.place_particle(Kp, [0,0,0])
    s.place_particle(KKi, [0,0,sigma+1e-20])

    s.throw_in_particles(KK, N_KK-1, box1)

    # print kD
    # print k_a(Mtom3(0.02e9), kD)
    # print k_a(Mtom3(0.032e9), kD)
    # sys.exit(0)

#     end_time = 0
#     while 1:
#         s.step()
#         next_time = s.scheduler.getTopTime()
#         if next_time > end_time:
#             s.stop(end_time)
#             break

#     s.reset()
#     k1 = k_a(Mtom3(0.02e9), kD)
#     k2 = k_d(1.0, Mtom3(0.02e9), kD)
#     k3 = 1.5
    k4 = k_a(Mtom3(0.032e9), kD)
#     k5 = k_d(1.0, Mtom3(0.032e9), kD)
#     k6 = 15.0

#     r1 = create_binding_reaction_rule(K, KK, K_KK, k1)
#     m.network_rules.add_reaction_rule(r1)
#     r2 = create_unbinding_reaction_rule(K_KK, K, KK, k2)
#     m.network_rules.add_reaction_rule(r2)
#     r3 = create_unbinding_reaction_rule(K_KK, Kp, KKi, k3)
#     m.network_rules.add_reaction_rule(r3)

    r4 = create_binding_reaction_rule(Kp, KK, Kp_KK, k4)
    m.network_rules.add_reaction_rule(r4)
#     r5 = create_unbinding_reaction_rule(Kp_KK, Kp, KK, k5)
#     m.network_rules.add_reaction_rule(r5)
#     r6 = create_unbinding_reaction_rule(Kp_KK, Kpp, KKi, k6)
#     m.network_rules.add_reaction_rule(r6)


#     r7 = create_binding_reaction_rule(Kpp, P, Kpp_P, k1)
#     m.network_rules.add_reaction_rule(r7)
#     r8 = create_unbinding_reaction_rule(Kpp_P, Kpp, P, k2)
#     m.network_rules.add_reaction_rule(r8)
#     r9 = create_unbinding_reaction_rule(Kpp_P, Kp, Pi, k3)
#     m.network_rules.add_reaction_rule(r9)
    
#     r10 = create_binding_reaction_rule(Kp, P, Kp_P, k4)
#     m.network_rules.add_reaction_rule(r10)
#     r11 = create_unbinding_reaction_rule(Kp_P, Kp, P, k5)
#     m.network_rules.add_reaction_rule(r11)
#     r12 = create_unbinding_reaction_rule(Kp_P, K, Pi, k6)
#     m.network_rules.add_reaction_rule(r12)


    r13 = create_unimolecular_reaction_rule(KKi, KK, ki)
    m.network_rules.add_reaction_rule(r13)
#     r14 = create_unimolecular_reaction_rule(Pi, P, ki)
#     m.network_rules.add_reaction_rule(r14)

    s.set_model(m)


#     logname = model + '_' + '_'.join(sys.argv[1:6]) + '_' +\
#         os.environ['SGE_TASK_ID']

#     outfile = open('data/' + logname + '_t.dat', 'w')


    while s.t < T:
        s.step()

        if s.last_reaction:
            r = s.last_reaction
            for p in r.products:
#                if p.species == Kpp:
                if p.species == Kp_KK:
                    if s.t <= T:
                        return s.t
                    else:
                        return -1
        if s.get_next_time() > T:
            return -1

    return -1

    
if __name__ == '__main__':

    import os

    outfilename = 'data/model3-smallt_' + '_'.join(sys.argv[1:3]) +\
        '_' + os.environ['SGE_TASK_ID']

    def runmain():
        run(outfilename, float(sys.argv[1]), 
            float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]))



    runmain()
#     try:
#         import cProfile as profile
#     except:
#         import profile
#     profile.run('runmain()', 'fooprof')
        

#     import pstats
#     pstats.Stats('fooprof').sort_stats('time').print_stats(40)

