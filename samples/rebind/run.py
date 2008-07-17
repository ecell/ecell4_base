#!/usr/bin/env python

'''
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py a 1e-5 1 10000 100

'''


from egfrd import *
from bd import *

def run( outfilename, T, DX_factor, N_X, seq, N ):
    print outfilename

    outfile_r = open( outfilename + '_r.dat', 'w' )
    outfile_t = open( outfilename + '_t.dat', 'w' )

    for i in range( N ):
        d, t, t_a = singlerun( T, DX_factor, N_X )
        outfile_r.write( '%g\n' % d )
        outfile_r.flush()
        if t_a != 0.0:
            outfile_t.write( '%g\n' % t_a )
            outfile_t.flush()

        print i, d, t
        assert d == 0 or t == T

    outfile_t.close()
    outfile_r.close()



def singlerun( T, DX_factor, N_X ):

    s = EGFRDSimulator()
    #s.setUserMaxShellSize( 1e-6 )
    #s = BDSimulator()


    # 100 nM = 100e-9 * N_A * 100 / m^3 = 6.02e19
    # V = 1 / 6.02e19 = 1.66e-20 m^3
    # L = 2.55e-7 m

    # 1 uM = 6.02e20 / m^3
    # V = 1.66e-21 m^3
    # L = 1.18e-7

    V = 1.66e-21 # m^3
    L = V ** (1.0/3.0) 

    s.setWorldSize( L )

    matrixSize = min( max( 3, int( (3 * N_X) ** (1.0/3.0) ) ), 60 )
    print 'matrixSize=', matrixSize
    s.setMatrixSize( matrixSize )

    box1 = CuboidalSurface( [0,0,0],[L,L,L] )

    radius = 2.5e-9
    sigma = radius * 2
    r0 = sigma
    D = 1e-12
    D_tot = D * 2

    tau = sigma**2 / D
    print 'tau=', tau

    kf = 10 * sigma * D_tot

    A = Species( 'A', D, radius )
    s.addSpecies( A )
    B = Species( 'B', D, radius )
    s.addSpecies( B )
    C = Species( 'C', D, radius )
    s.addSpecies( C )

    DX = D * DX_factor

    X = Species( 'X', DX, radius )
    s.addSpecies( X )

    if N_X != 0:
        s.throwInParticles( X, N_X, box1 )

    endTime = tau * 1
    while 1:
        s.step()
        nextTime = s.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break


    s.reset()

    r1 = BindingReactionType( A, B, C, kf )
    s.addReactionType( r1 )

    r2 = UnbindingReactionType( C, A, B, 1e3 )
    s.addReactionType( r2 )

    A_pos = [0,0,0]
    B_pos = [(A.radius + B.radius)+1e-23,0,0]

    while 1:
        pp = s.getParticlesWithinRadius( A_pos, A.radius )
        if not pp:
            break
        for p in pp:
            s.removeParticle( p )
        s.throwInParticles( X, len( pp ), box1 )

    s.placeParticle( A, A_pos )

    while 1:
        pp = s.getParticlesWithinRadius( B_pos, B.radius )
        if not pp:
            break
        for p in pp:
            s.removeParticle( p )
        s.throwInParticles( X, len( pp ), box1 )    

    s.placeParticle( B, B_pos )


    endTime = T
    s.step()

    t_a = 0.0
    while 1:
        if s.populationChanged and t_a == 0.0:
            t_a = s.t
        #    print 'reaction'
        #    return 0.0, s.t
        nextTime = s.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()

    if C.pool.size != 0:
        return 0, s.t, t_a

    distance = s.distance( A.pool.positions[0], B.pool.positions[0] )

    return distance, s.t, t_a
    
if __name__ == '__main__':

    outfilename = 'data/rebind_' + '_'.join( sys.argv[1:5] )
    run( outfilename, float( sys.argv[1] ), float( sys.argv[2] ), 
         int( sys.argv[3] ), int( sys.argv[4] ), int( sys.argv[5] )  )
