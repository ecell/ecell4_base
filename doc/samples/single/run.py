#!/usr/bin/env python

from egfrd import *

def run( outfilename, T, S, N ):
    print outfilename

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun( T, S )
        outfile.write( '%g\n' % d )

        print d, t
        assert d == 0 or t == T

    outfile.close()



def singlerun( T, S ):

    s = EGFRDSimulator()
    s.setCellSize( 1e-3 )

    s.setMaxShellSize( S )

    A = Species( 'A', 1e-11, 5e-8 )
    s.addSpecies( A )
    
    particleA = s.placeParticle( A, [0,0,0] )

    endTime = T
    s.step()

    while 1:
        nextTime = s.scheduler.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()

    distance = s.distance( [0,0,0], particleA.getPos() )

    return distance, s.t
    
if __name__ == '__main__':
    run( sys.argv[1], float( sys.argv[2] ), float( sys.argv[3] ),
         int( sys.argv[4] ) )
