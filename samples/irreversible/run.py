#!/usr/bin/env python

'''
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.3.out 5e-2 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.2.out 5e-3 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.1.out 5e-4 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.0.out 5e-5 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-1.out 5e-6 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-2.out 5e-7 1000000 &
LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py irr.-3.out 5e-8 1000000 &
'''


from egfrd import *
#from bd import *

def run( outfilename, T, N ):

    outfile = open( outfilename, 'w' )

    for i in range( N ):
        d, t = singlerun2( T )
        outfile.write( '%g\n' % d )
        outfile.flush()
        print i
        #print d, t
        assert d == 0 or t == T


    outfile.close()



def singlerun1( T ):

    s = BDSimulator()
    s.setWorldSize( 1e-3 )

    #s.setMaxShellSize( 1e-6 )


    sigma = 1e-8
    r0 = sigma
    D = 2e-12
    kf = 10 * sigma * D

    A = Species( 'A', 0.0, sigma/2 )
    s.addSpecies( A )
    B = Species( 'B', D, sigma/2 )
    s.addSpecies( B )
    C = Species( 'C', 0.0, sigma/2 )
    s.addSpecies( C )

    r1 = BindingReactionType( A, B, C, kf )
    s.addReactionType( r1 )
    
    particleA = s.placeParticle( A, [0,0,0] )
    particleB = s.placeParticle( B, [(A.radius + B.radius)+1e-23,0,0] )

    endTime = T
    s.step()

    while 1:
        nextTime = s.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()
        if s.populationChanged():
            print 'reaction'
            return 0.0, s.t

    distance = s.distance( particleB.getPos(), particleA.getPos() )

    return distance, s.t


def singlerun2( T ):

    s = EGFRDSimulator()
    s.setWorldSize( 1e-3 )

    #s.setUserMaxShellSize( 1e-6 )
    s.setUserMaxShellSize( 1e-3 )

    sigma = 1e-8
    r0 = sigma
    D = 1e-12
    kf = 10 * sigma * D
    tau = sigma*sigma/D 

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

    while 1:
        s.step()
        if s.populationChanged:
            #print 'reaction'
            return 0.0, s.t

        nextTime = s.getNextTime()
        if nextTime > endTime:
            s.stop( endTime )
            break

    distance = s.distance( particleB.getPos(), particleA.getPos() )

    return distance, s.t

    

if __name__ == '__main__':
    run( sys.argv[1], float( sys.argv[2] ), int( sys.argv[3] ) )
