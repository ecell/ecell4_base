#!/usr/bin/env python

'''
# D_factor ti T N

LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py 1 1 100 10

'''


from egfrd import *
from bd import *

def run( outfilename, D_factor, ti, T, N ):
    print outfilename

    outfile_t = open( outfilename + '_t.dat', 'w' )

    outfile_t.write( '%d\n' % N )

    for i in range( N ):
        t = singlerun( D_factor, ti, T )
        #print i, t

        if t != -1:
            outfile_t.write( '%g\n' % t )
            outfile_t.flush()
 

        #[ outfile_r.flush() for outfile_r in outfile_r_list ]


    outfile_t.close()
    #[ outfile_r.close() for outfile_r in outfile_r_list ]



def singlerun( D_factor, ti, T ):

    V = 1e-15
    #V = 1e-16
    D_ratio = 1
    
    if ti == 0:
        ki = float( 'inf' )
    else:
        ki = math.log( 2 ) / ti


    D_ref = 1e-12

    D_move = D_ref * D_factor

    D_react = D_ref

    # V in liter, L in meter
    L = math.pow( V * 1e-3, 1.0 / 3.0 )

    s = EGFRDSimulator()
    s.setWorldSize( L )

    N = 180
    matrixSize = min( max( 3, int( (3 * N) ** (1.0/3.0) ) ), 60 )
    s.setMatrixSize( matrixSize )


    box1 = CuboidalSurface( [0,0,0],[L,L,L] )

    radius = 2.5e-9

#     K = Species( 'K', D_move, radius )
#     s.addSpecies( K )
    KK = Species( 'KK', D_move, radius )
    s.addSpecies( KK )
#     P = Species( 'P', D_move, radius )
#     s.addSpecies( P )
    Kp = Species( 'Kp', D_move, radius )
    s.addSpecies( Kp )
#     Kpp = Species( 'Kpp', D_move, radius )
#     s.addSpecies( Kpp )
#     K_KK = Species( 'K_KK', D_move, radius )
#     s.addSpecies( K_KK )
    Kp_KK = Species( 'Kp_KK', D_move, radius )
    s.addSpecies( Kp_KK )
#     Kpp_KK = Species( 'Kpp_KK', D_move, radius )
#     s.addSpecies( Kpp_KK )
#     Kpp_P = Species( 'Kpp_P', D_move, radius )
#     s.addSpecies( Kpp_P )
#     Kp_P = Species( 'Kp_P', D_move, radius )
#     s.addSpecies( Kp_P )

    # inactive forms
    KKi = Species( 'KKi', D_move, radius )
    s.addSpecies( KKi )
#     Pi = Species( 'Pi', D_move, radius )
#     s.addSpecies( Pi )



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
    kD = k_D( D_react * 2, sigma )

    N_K = C2N( 200e-9, V ) 
    N_KK = C2N( 50e-9, V )
    N_P = C2N( 50e-9, V )

    #print N_KK
    #s.throwInParticles( K, N_K, box1 )
    #s.throwInParticles( KK, N_KK, box1 )
    #s.throwInParticles( P, N_P, box1 )
    
    s.placeParticle( Kp, [0,0,0] )
    s.placeParticle( KKi, [0,0,sigma+1e-23] )

#    s.throwInParticles( KK, N_KK-1, box1 )

    # print kD
    # print k_a( Mtom3( 0.02e9 ), kD )
    # print k_a( Mtom3( 0.032e9 ), kD )
    # sys.exit(0)

#     endTime = 0
#     while 1:
#         s.step()
#         nextTime = s.scheduler.getTopTime()
#         if nextTime > endTime:
#             s.stop( endTime )
#             break

#     s.reset()
#     k1 = k_a( Mtom3( 0.02e9 ), kD )
#     k2 = k_d( 1.0, Mtom3( 0.02e9 ), kD )
#     k3 = 1.5
    k4 = k_a( Mtom3( 0.032e9 ), kD )
#     k5 = k_d( 1.0, Mtom3( 0.032e9 ), kD )
#     k6 = 15.0

#     r1 = BindingReactionType( K, KK, K_KK, k1 )
#     s.addReactionType( r1 )
#     r2 = UnbindingReactionType( K_KK, K, KK, k2 )
#     s.addReactionType( r2 )
#     r3 = UnbindingReactionType( K_KK, Kp, KKi, k3 )
#     s.addReactionType( r3 )

    r4 = BindingReactionType( Kp, KK, Kp_KK, k4 )
    s.addReactionType( r4 )
#     r5 = UnbindingReactionType( Kp_KK, Kp, KK, k5 )
#     s.addReactionType( r5 )
#     r6 = UnbindingReactionType( Kp_KK, Kpp, KKi, k6 )
#     s.addReactionType( r6 )


#     r7 = BindingReactionType( Kpp, P, Kpp_P, k1 )
#     s.addReactionType( r7 )
#     r8 = UnbindingReactionType( Kpp_P, Kpp, P, k2 )
#     s.addReactionType( r8 )
#     r9 = UnbindingReactionType( Kpp_P, Kp, Pi, k3 )
#     s.addReactionType( r9 )
    
#     r10 = BindingReactionType( Kp, P, Kp_P, k4 )
#     s.addReactionType( r10 )
#     r11 = UnbindingReactionType( Kp_P, Kp, P, k5 )
#     s.addReactionType( r11 )
#     r12 = UnbindingReactionType( Kp_P, K, Pi, k6 )
#     s.addReactionType( r12 )


    r13 = UnimolecularReactionType( KKi, KK, ki )
    s.addReactionType( r13 )
#     r14 = UnimolecularReactionType( Pi, P, ki )
#     s.addReactionType( r14 )


#     logname = model + '_' + '_'.join( sys.argv[1:6] ) + '_' +\
#         os.environ[ 'SGE_TASK_ID' ]

#     outfile = open( 'data/' + logname + '_t.dat', 'w' )


    while s.t < T:
        s.step()

        if s.lastReaction:
            r = s.lastReaction
            for p in r.products:
#                if p.species == Kpp:
                if p.species == Kp_KK:
                    if s.t <= T:
                        return s.t
                    else:
                        return -1
        if s.getNextTime() > T:
            return -1

    return -1

    
if __name__ == '__main__':

    import os

    outfilename = 'data/model3-smallt_' + '_'.join( sys.argv[1:3] ) +\
        '_' + os.environ['SGE_TASK_ID']

    def runmain():
        run( outfilename, float( sys.argv[1] ), 
             float(sys.argv[2]), float( sys.argv[3] ), int( sys.argv[4] ) )



    runmain()
#     try:
#         import cProfile as profile
#     except:
#         import profile
#     profile.run('runmain()', 'fooprof')
        

#     import pstats
#     pstats.Stats('fooprof').sort_stats('time').print_stats(40)

