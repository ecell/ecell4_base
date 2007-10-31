#!/usr/bin/env python

from egfrd import *

import gc

def run( outfilename, N ):
    print outfilename

    outfile = open( outfilename, 'w' )

    T = .1

    for i in range( N ):
        d, t = singlerun( T )
        outfile.write( '%g\n' % d )

        print d, t
        assert d == 0 or t == T

        #gc_classlist()

    outfile.close()



def singlerun( T ):

    s = EGFRDSimulator()
    s.setCellSize( 1e-3 )

    s.setMaxShellSize( 1e-6 )

    A = Species( 'A', 0.0, 5e-8 )
    s.addSpecies( A )
    B = Species( 'B', 1e-11, 5e-8 )
    s.addSpecies( B )
    C = Species( 'C', 0.0, 5e-8 )
    s.addSpecies( C )
    
    r1 = BindingReactionType( A, B, C, 1e6 / N_A )
    s.addReactionType( r1 )
    
    particleA = s.placeParticle( A, [0,0,0] )
    particleB = s.placeParticle( B, [1e-7 * (1.0 + 1e-8),0,0] )

    endTime = T
    s.step()

    while 1:
        nextTime = s.scheduler.getTopEvent().getTime()
        if nextTime > endTime:
            s.stop( endTime )
            break
        s.step()
        if s.populationChanged():
            print 'reaction'
            t = s.t
            del s
            return 0.0, t

    print particleA.getPos()
    assert s.distance( particleA.getPos(), [0,0,0] ) == 0.0

    distance = s.distance( particleB.getPos(), [0,0,0] )
    t = s.t

    #gc.collect()
    #print 'GC referrers', gc.get_referrers( s)
    #for i in  gc.get_referrers( s):
    #print 'GC i ref', i, gc.get_referrers( i )
    #print 'GC referrents', gc.get_referents( s)


    del s
    return distance, t
    

import gc

def gc_classlist():
    gc.collect()
    objCount = {}
    objList = gc.get_objects()

    tuplist = []
    for obj in objList:
        if getattr(obj, "__class__", None):
            name = obj.__class__.__name__
            if objCount.has_key(name):
                objCount[name] += 1
            else:
                objCount[name] = 1
            if name == 'tuple':
                tuplist.append( str( obj ) )

                
    tuplist.sort()
    #print 'GC', tuplist

    objs = objCount.items()
    objs.sort(key=lambda x: x[1],reverse=True) 
    print 'GC', objs
    #print 'GC', objs[0]
    del objList, objs, objCount

if __name__ == '__main__':
    run( sys.argv[1], int( sys.argv[2] ) )
