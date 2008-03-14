#!/usr/bin/env python

from bd import *

def run( outfilename, T, N ):

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun( T )
        outfile.write( '%g\n' % d )

        print d, t
        #assert d == 0 or t == T

    outfile.close()


def singlerun( T ):

    s = BDSimulator()
    s.setCellSize( 1e-3 )

    sigma = 1e-8
    r0 = sigma
    D = 1e-12
    kf = 10 * sigma * D

    A = Species( 'A', D, sigma/2 )
    s.addSpecies( A )
    B = Species( 'B', D, sigma/2 )
    s.addSpecies( B )
    C = Species( 'C', D, sigma/2 )
    s.addSpecies( C )
    
    r1 = BindingReactionType( A, B, C, kf )
    s.addReactionType( r1 )
    
    particleA = s.placeParticle( A, [0,0,0] )
    particleB = s.placeParticle( B, [(A.radius + B.radius)+1e-23,0,0] )

    endTime = T
    #s.initialize()

    while 1:
        nextTime = s.t + s.dt
        if nextTime > endTime:
            break
        s.step()
        #print s.populationChanged()
        if s.populationChanged():
            print 'reaction'
            return 0.0, s.t

    distance = s.distance( particleB.getPos(), particleA.getPos() )

    return distance, s.t

    

if __name__ == '__main__':
    run( sys.argv[1], float( sys.argv[2] ), int( sys.argv[3] ) )
