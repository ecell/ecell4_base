#!/usr/bin/env python

from egfrd import *

def run( outfilename, T, N ):
    print outfilename

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun( T )
        outfile.write( '%g\n' % d )

        print d, t
        assert d == 0 or t == T

    outfile.close()



def singlerun( T ):

    s = EGFRDSimulator()
    s.setCellSize( 1e-3 )

    s.setMaxShellSize( 1e-6 )

    A = Species( 'A', 5e-12, 5e-8 )
    s.addSpecies( A )
    B = Species( 'B', 5e-12, 5e-8 )
    s.addSpecies( B )
    C = Species( 'C', 5e-12, 5e-8 )
    s.addSpecies( C )
    
    r1 = BindingReactionType( A, B, C, 1e6 / N_A )
    s.addReactionType( r1 )

    r2 = UnbindingReactionType( C, A, B, 1e3 )
    s.addReactionType( r2 )

    particleA = s.placeParticle( A, [0,0,0] )
    particleB = s.placeParticle( B, [(A.radius + B.radius)+1e-23,0,0] )

    endTime = T
    s.step()

    while 1:
        nextTime = s.scheduler.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()

    if C.pool.size != 0:
        return 0, s.t

    distance = s.distance( particleB.getPos(), particleA.getPos() )

    return distance, s.t
    
if __name__ == '__main__':
    run( sys.argv[1], float( sys.argv[2] ), int( sys.argv[3] ) )
