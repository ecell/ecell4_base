#!/usr/bin/env python

'''
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.3.out 0.05 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.2.out 0.005 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.1.out 0.0005 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.0.out 5e-5 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.-1.out 5e-6 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.-2.out 5e-7 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py rev.-3.out 5e-8 1000000 &
'''


from egfrd import *
from bd import *

def run( outfilename, T, N ):
    print outfilename

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun( T )
        outfile.write( '%g\n' % d )
        outfile.flush()
        #print d, t
        assert d == 0 or t == T

    outfile.close()



def singlerun( T ):

    s = EGFRDSimulator()
    #s.setUserMaxShellSize( 1e-6 )
    #s = BDSimulator()

    s.setWorldSize( 1e-3 )

    sigma = 5e-9
    r0 = sigma
    D = 1e-12
    D_tot = D * 2
    kf = 10 * sigma * D_tot

    A = Species( 'A', D, sigma/2 )
    s.addSpecies( A )
    B = Species( 'B', D, sigma/2 )
    s.addSpecies( B )
    C = Species( 'C', D, sigma/2 )
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
