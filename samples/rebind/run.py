#!/usr/bin/env python

'''
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py a 1e-5 1 10000 100

'''


from egfrd import *
from bd import *

def run( outfilename, T, DX_factor, N_X, N ):
    print outfilename

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun( T, DX_factor, N_X )
        outfile.write( '%g\n' % d )
        outfile.flush()
        print i, d, t
        assert d == 0 or t == T

    outfile.close()



def singlerun( T, DX_factor, N_X ):

    s = EGFRDSimulator()
    #s.setUserMaxShellSize( 1e-6 )
    #s = BDSimulator()


    # 100 nM = 100e-9 * N_A * 100 / m^3 = 6.02e19
    # V = 1 / 6.02e19 = 1.66e-20 m^3
    V = 1.66e-20 # m^3
    L = V ** (1.0/3.0) # 2.55e-7 m

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

    s.throwInParticles( X, N_X, box1 )

    endTime = tau * 10
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

    while 1:
        nextTime = s.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()

    if C.pool.size != 0:
        return 0, s.t

    distance = s.distance( A.pool.positions[0], B.pool.positions[0] )

    return distance, s.t
    
if __name__ == '__main__':
    run( sys.argv[1], float( sys.argv[2] ), float( sys.argv[3] ), 
         int( sys.argv[4] ), int( sys.argv[5] )  )
