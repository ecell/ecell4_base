#!/usr/bin/env python

# tau = 0.0001
#PYTHONPATH=../../.. python run.py rev.0.out 0.0001 1000000

#from egfrd import *
from bd import *

def run( outfilename, T, N ):
    print outfilename

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun( T )
        outfile.write( '%g\n' % d )
        outfile.flush()
        print d, t
        assert d == 0 or t == T

    outfile.close()



def singlerun( T ):

    #s = EGFRDSimulator()
    #s.setMaxShellSize( 1e-6 )

    s = BDSimulator()

    s.setCellSize( 1e-3 )

    sigma = 1e-8
    r0 = sigma
    D = 1e-12
    kf = 10 * sigma * D

    A = Species( 'A', D/2, sigma/2 )
    s.addSpecies( A )
    B = Species( 'B', D/2, sigma/2 )
    s.addSpecies( B )
    C = Species( 'C', D/2, sigma/2 )
    s.addSpecies( C )

    r1 = BindingReactionType( A, B, C, kf )
    s.addReactionType( r1 )

    r2 = UnbindingReactionType( C, A, B, 1e3 )
    s.addReactionType( r2 )

    s.placeParticle( A, [0,0,0] )
    s.placeParticle( B, [(A.radius + B.radius)+1e-23,0,0] )

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
    run( sys.argv[1], float( sys.argv[2] ), int( sys.argv[3] ) )
