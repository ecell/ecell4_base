#!/usr/env python


import math
import random
import sys
import operator

import numpy
import scipy
import scipy.optimize


from utils import *
#from surface import *
from _gfrd import *


N_A = 6.0221367e23


def p_free( r, t, D ):
    Dt4 = D * t * 4.0
    Pi4Dt = numpy.pi * Dt4
    rsq = r * r
    
    p = math.exp( - rsq / Dt4 ) / math.sqrt( Pi4Dt * Pi4Dt * Pi4Dt )

    jacobian = 4.0 * numpy.pi * rsq

    return p * jacobian
    

def drawR_free( t, D ):
    ro = math.sqrt( 2.0 * D * t )

    displacement = numpy.random.normal( 0.0, ro, 3 )

    return displacement

class OverlapError( Exception ):
    pass


class Species( object ):
    
    def __init__( self, id, D, radius ):
        self.id = id
        self.D = D
        self.radius = radius
        self.pool = ParticlePool()

    def newParticle( self, position ):

        serial = self.pool.newParticle( position )
        return serial


    def removeParticleByIndex( self, index ):
        self.pool.removeByIndex( index )

    def removeParticleBySerial( self, serial ):
        self.pool.removeBySerial( serial )


class ReactionType( object ):

    def __init__( self, reactants=[], products=[], k=0.0 ):
        self.reactants = reactants
        self.products = products
        self.k = k

    def order( self ):
        return len( self.reactants )

    def str( self ):
        s = ''
        for i in self.reactants:
            s += i.id
            s += ' '

        s += ' -> '

        for i in self.products:
            s += i.id
            s += ' '

        return s

class UnimolecularReactionType( ReactionType ):

    def __init__( self, s1, p1, k ):
        ReactionType.__init__( self, [ s1, ], [ p1, ], k )

class DecayReactionType( ReactionType ):

    def __init__( self, s1, k ):
        ReactionType.__init__( self, [ s1, ], [], k )


class BindingReactionType( ReactionType ):

    def __init__( self, s1, s2, p1, k ):
        ReactionType.__init__( self, [ s1, s2 ], [ p1, ], k )
        D = s1.D + s2.D
        sigma = s1.radius + s2.radius
        self.pairGreensFunction = BasicPairGreensFunction( D, k, sigma )

class RepulsionReactionType( ReactionType ):

    def __init__( self, s1, s2 ):
        ReactionType.__init__( self, [ s1, s2 ], [], 0.0 )

        D = s1.D + s2.D
        sigma = s1.radius + s2.radius
        self.pairGreensFunction = BasicPairGreensFunction( D, self.k,\
                                                                sigma )

class UnbindingReactionType( ReactionType ):

    def __init__( self, s1, p1, p2, k ):
        ReactionType.__init__( self, [ s1, ], [ p1, p2 ], k )


class Particle( object ):

    def __init__( self, species, serial=None, index=None ):

        self.species = species

        if serial != None and index == None:
            self.serial = serial
        elif serial == None and index != None:
            self.serial = species.pool.getSerialByIndex( index )
        else:
            raise ValueError, 'give either serial or index.'

        self.pool = self.species.pool

    def __str__( self ):
        return str( ( self.species.id, self.serial ) )

    def __cmp__( self, other ):
        if self.species == other.species:
            return self.serial - other.serial
        elif self.species < other.species:
            return -1
        else:
            return 1

    def getPos( self ):
        pool = self.pool
        return pool.positions[ pool.indexMap[ self.serial ] ]

    def setPos( self, newpos ):
        pool = self.pool
        pool.positions[ pool.indexMap[ self.serial ] ] = newpos

    pos = property( getPos, setPos )

    def getIndex( self ):
        return self.pool.indexMap[ self.serial ]


class ParticlePool( object ):

    def __init__( self ):

        self.indexMap = {}
        self.serialCounter = 0

        self.serials = numpy.array( [0,], numpy.integer )

        self.positions = numpy.array( [], numpy.floating )
        self.positions.shape = ( 0, 3 )

        self.size = 0

    def newParticle( self, position ):

        newindex = self.size
        newserial = self.serialCounter
        
        self.size += 1
        self.serialCounter += 1

        self.__resizeArrays( self.size )

        self.indexMap[ newserial ] = newindex
        self.serials[ newindex ] = newserial
        self.positions[ newindex ] = position

        return newserial
    

    def __resizeArrays( self, newsize ):

        self.serials.resize( newsize )
        self.positions.resize( ( newsize, 3 ) )


    def removeBySerial( self, serial ):

        index = self.getIndex( serial )
        self.removeByIndex( index )


    def removeByIndex( self, index ):

        self.size -= 1

        serialOfRemovedItem = self.serials[ index ]
        serialOfMovedItem = self.serials[ self.size ]

        # swap items at index and the end, and discard the item
        # to be removed, which is now at the end, by shrinking
        # the arrays by one.
        self.serials[ index ] = self.serials[ self.size ]
        self.positions[ index ] = self.positions[ self.size ]
        self.__resizeArrays( self.size )

        # book keeping
        del self.indexMap[ serialOfRemovedItem ]
        self.indexMap[ serialOfMovedItem ] = index

    def getSerialByIndex( self, index ):
        
        return self.serials[index]

    def getIndex( self, serial ):

        return self.indexMap[ serial ]



class GFRDSimulatorBase( object ):
    
    def __init__( self ):
        self.speciesList = {}
        self.reactionTypeMap1 = {}
        self.reactionTypeMap2 = {}

        self.surfaceList = []

        self.dt = 1e-7
        self.t = 0.0

        self.H = 3.0
        
        self.dtLimit = 1e-3
        self.dtMax = self.dtLimit

        self.nextReaction = None


        self.setCellSize( INF )


        # counters
        self.rejectedMoves = 0
        self.reactionEvents = 0


    def getTime( self ):
        return self.t


    def setCellSize( self, size ):

        self.fsize = size

        if self.fsize == INF:
            self._distanceSq = distanceSq_Simple
            self._distanceSqArray = distanceSqArray_Simple
        else:
            self._distanceSq = distanceSq_Cyclic
            self._distanceSqArray = distanceSqArray_Cyclic

    def getCellSize( self ):
        return self.fsize

    def applyBoundary( self, pos ):
        if self.fsize != INF:
            pos %= self.fsize

        return pos

    def getReactionType1( self, species ):
        return self.reactionTypeMap1.get( species, None )

    def getReactionType2( self, species1, species2 ):
        return self.reactionTypeMap2.get( ( species1, species2 ), None )

    def getSpeciesByIndex( self, i ):
        return self.speciesList.values()[i]

    def simpleDiffusion( self, speciesIndex, particleIndex ):

        species = self.speciesList.values()[speciesIndex]
        if species.D == 0.0:
            return

        limitSq = self.H * self.H * ( 6.0 * species.D * self.dt )

        while True:
            displacement = drawR_free( self.dt, species.D )
            
            distSq = ( displacement * displacement ).sum()

            if distSq <= limitSq:
                break

            self.rejectedMoves += 1

        pos = species.pool.positions[particleIndex]
        pos += displacement

        #FIXME: SURFACE
        pos %= self.fsize


    def simpleDiffusionWithSurface( self, speciesIndex, particleIndex,\
                                    surface ):
        species = self.speciesList.values()[speciesIndex]

        limitSq = self.H * self.H * ( 6.0 * species.D * self.dt )

        while True:
            displacement = drawR_free( self.dt, species.D )
            
            distSq = ( displacement * displacement ).sum()

            if distSq <= limitSq:
                break

            self.rejectedMoves += 1

        pos = species.pool.positions[particleIndex]
        pos += displacement

        #SURFACE
        print surface
        


    def distanceSq( self, position1, position2 ):
        return self._distanceSq( position1, position2, self.fsize )

    def distance( self, position1, position2 ):
        return numpy.sqrt( self.distanceSq( position1, position2 ) )
        
    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.fsize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1,\
                                                 positions ) )

    def addSurface( self, surface ):
        self.surfaceList.append( surface )


    def addSpecies( self, species ):
        self.speciesList[ species.id ] = species

    def addReactionType( self, r ):

        numReactants = len( r.reactants )

        if numReactants == 1:
            species1 = r.reactants[0]
            if self.reactionTypeMap1.has_key( species1 ):
                self.reactionTypeMap1[species1].append( r )
            else:
                self.reactionTypeMap1[species1] = [ r, ]

        elif numReactants == 2:
            species1 = r.reactants[0]
            species2 = r.reactants[1]
            self.reactionTypeMap2[ (species1,species2) ] = r
            if species1 != species2:
                self.reactionTypeMap2[ (species2,species1) ] = r

        else:
            raise RuntimeError, 'Invalid ReactionType.'

    def setAllRepulsive( self ):
        for species1 in self.speciesList.values():
            for species2 in self.speciesList.values():
                try:
                    _ = self.reactionTypeMap2[ ( species1, species2 ) ]
                except:
                    self.reactionTypeMap2[ ( species1, species2 ) ] =\
                                            RepulsionReactionType( species1,\
                                                                   species2 )
        

    def throwInParticles( self, species, n, surface=[] ):
        print 'throwing in %s %s particles' % ( n, species.id )

        for _ in range( int( n ) ):

            while True:

                #position= numpy.random.uniform( 0, self.fsize, 3 )
                position = surface.randomPosition()
                if self.checkOverlap( position, species.radius ):
                    break
                else:
                    print _
            
            species.newParticle( position )


    def placeParticle( self, species, pos ):

        radius = species.radius

        if not self.checkOverlap( pos, radius ):
            raise OverlapError, 'overlap check failed'
            
        particle = self.createParticle( species, pos )
        return particle


    def createParticle( self, species, pos ):
        newserial = species.newParticle( pos )
        return Particle( species, serial=newserial )

    def removeParticle( self, particle ):
        particle.species.pool.removeBySerial( particle.serial )

        
    def checkOverlap( self, position, radius ):
        
        for species2 in self.speciesList.values():

            if species2.pool.size == 0:
                continue

            positions2 = species2.pool.positions

            radius12 = radius + species2.radius
            radius12sq = radius12 * radius12

            closestdistsq = self.distanceSqArray( position, positions2 ).min()

            if closestdistsq <= radius12sq:
                print 'reject: closest distance = ', math.sqrt(closestdistsq),\
                    ', which must be at least ', radius12
                return False

        return True
    

    def clear( self ):

        self.dtMax = self.dtLimit
        self.dt = self.dtLimit

        self.nextReaction = None

        self.pairs = []
        self.singles = []



    def nextReactionTime1( self, rt, pool ):

        if pool.size == 0:
            return None, None

        dts = numpy.array( [ random.random() for i in range( pool.size ) ] )
        dts = - numpy.log( dts ) / rt.k

        i = numpy.argmin( dts )

        return dts[i], i 
        


    def nextReactionTime2( self, pair ):

        ( dt, ndt, s1, i1, s2, i2, rt ) = pair

        if rt == None:
            return INF

        species1 = self.speciesList.values()[s1]
        species2 = self.speciesList.values()[s2]
        pos1 = species1.pool.positions[i1]
        pos2 = species2.pool.positions[i2]

        radius = species1.radius + species2.radius
        r0 = self.distance( pos1, pos2 )

        if radius > r0:
            print 'CRITICAL: radius > r0', str(radius), str(r0)
            #            return scipy.Inf
            print pair
            print rt.str()
            #sys.exit(-1)

        u = random.random()


        res = rt.pairGreensFunction.drawTime( u, r0, self.dtMax )

        return res


        
    '''
    Get closest n Particles.

    When the optional argument speciesList is given, only Particles of
    species in the list are considered.  When speciesList is not given
    or is None, all species in the simulator are considered.
    
    This method returns a tuple ( neighbors, distances ), where neighbors
    is a list of Particle objects.
    '''


    def getNeighborParticles( self, pos, n=2, speciesList=None ):

        neighbors = []
        distances = []

        if not speciesList:
            speciesList = self.speciesList.values()

        for species in speciesList:

            # empty
            if species.pool.size == 0:
                continue

            dist = self.distanceSqArray( pos, species.pool.positions )
        
            indices = dist.argsort()[:n]
            dist = numpy.sqrt( dist.take( indices ) )
            dist -= species.radius

            distances.extend( dist )
            neighbors.extend( [ ( species, i ) for i in indices ] )

        topargs = numpy.argsort( distances )[:n]
        distances = numpy.take( distances, topargs )
        neighbors = [ neighbors[arg] for arg in topargs ]
        neighbors = [ Particle( arg[0], index=arg[1] ) for arg in neighbors ]

        return neighbors, distances


    def getClosestParticle( self, pos, ignore=[] ):

        neighbors, distances = self.getNeighborParticles( pos, 
                                                          len( ignore ) + 1 )

        for i in range( len( neighbors ) ): 
            if neighbors[i] not in ignore:
                closest, distance = neighbors[i], distances[i]

                assert not closest in ignore
                return closest, distance

        # default case: none left.
        return None, numpy.inf

    def checkSurfaces( self, speciesIndex1, particleIndex ):

        speciesList = self.speciesList.values()

        species = speciesList[ speciesIndex1 ]
        pos = species.pool.positions[ particleIndex ].copy()

        dist = [ surface.distance( pos ) for surface in self.surfaceList ]

        if len( dist ) == 0:
            return -1, 0.0

        idx = numpy.argmin( dist )
        dist = dist[idx]

        dt = ( dist - species.radius ) ** 2 / \
                ( self.H * self.H * 6.0 * species.D ) 
        
        return dt, idx
        
    def dumpPopulation( self ):
        buf = ''
        for species in self.speciesList.values():
            buf += species.id + ':' + str( species.pool.size ) + '\t'

        print buf
