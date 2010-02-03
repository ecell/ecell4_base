#!/usr/bin/env python

from bd import *
import sys

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
    s.setWorldSize( 1e-3 )

    sigma = 5e-9
    r0 = sigma
    D = 1e-12
    kf = 10 * sigma * D

    m = ParticleModel()

    A = m.new_species_type( 'A', D/2, sigma/2 )
    B = m.new_species_type( 'B', D/2, sigma/2 )
    C = m.new_species_type( 'C', D/2, sigma/2 )
    
    r1 = createBindingReactionRule( A, B, C, kf )
    s.newwork_rules.add_reaction_rule( r1 )

    s.setModel( m )
    
    particleA = s.placeParticle( A, [0,0,0] )
    particleB = s.placeParticle( B, [(A.radius + B.radius)+1e-23,0,0] )

    endTime = T
    #s.initialize()

    while 1:
        nextTime = s.t + s.dt
        if nextTime > endTime:
            break
        s.step()

        if s.lastReaction:
            print 'reaction'
            return 0.0, s.t

    distance = s.distance( particleB.getPos(), particleA.getPos() )

    return distance, s.t

    

if __name__ == '__main__':
    run( sys.argv[1], float( sys.argv[2] ), int( sys.argv[3] ) )
