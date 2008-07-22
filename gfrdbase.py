 #!/usr/env python


import math
import random
import sys

import numpy
import scipy


from utils import *
#from surface import *
from _gfrd import *

from cObjectMatrix import *

import os
import logging
import logging.handlers


def setupLogging():
    global log 

    log = logging.getLogger( 'pdpp' )

    if 'LOGFILE' in os.environ:
        if 'LOGSIZE' in os.environ and int( os.environ[ 'LOGSIZE' ] ) != 0:
            handler = logging.handlers.\
                RotatingFileHandler( os.environ[ 'LOGFILE' ], mode='w',
                                     maxBytes=int( os.environ[ 'LOGSIZE' ] ) )
        else:
            handler = logging.FileHandler( os.environ[ 'LOGFILE' ], 'w', )
            
    else:
        handler = logging.StreamHandler( sys.stdout )

    formatter = logging.Formatter( '%(message)s' )
    handler.setFormatter( formatter )
    log.addHandler( handler )
        
    LOGLEVELS = { 'CRITICAL': logging.CRITICAL,
                  'ERROR': logging.ERROR,
                  'WARNING': logging.WARNING,
                  'INFO': logging.INFO, 
                  'DEBUG': logging.DEBUG, 
                  'NOTSET': logging.NOTSET }

    if 'LOGLEVEL' in os.environ:
        log.setLevel( LOGLEVELS[ os.environ[ 'LOGLEVEL' ] ] )
    else:
        log.setLevel( logging.INFO )


setupLogging()


def p_free( r, t, D ):
    Dt4 = D * t * 4.0
    Pi4Dt = numpy.pi * Dt4
    rsq = r * r
    
    p = math.exp( - rsq / Dt4 ) / math.sqrt( Pi4Dt * Pi4Dt * Pi4Dt )

    jacobian = 4.0 * numpy.pi * rsq

    return p * jacobian
    

def drawR_free( t, D ):
    ro = math.sqrt( 2.0 * D * t )
    return numpy.random.normal( 0.0, ro, 3 )


class NoSpace( Exception ):
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
        self.check()


    def check( self ):
        if len( self.reactants ) == 1 and len( self.products ) <= 2:
            totalProductRadii = 0.0
            for product in self.products:
                totalProductRadii += product.radius
            if totalProductRadii > self.reactants[0].radius * 2:
                raise ValueError, \
                    'total product radii must be smaller than ' \
                    + 'reactant radius * 2'
        if self.k < 0:
            raise ValueError, 'k < 0'

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


class RepulsionReactionType( ReactionType ):

    def __init__( self, s1, s2 ):
        ReactionType.__init__( self, [ s1, s2 ], [], 0.0 )

        D = s1.D + s2.D
        sigma = s1.radius + s2.radius


class UnbindingReactionType( ReactionType ):

    def __init__( self, s1, p1, p2, k ):
        ReactionType.__init__( self, [ s1, ], [ p1, p2 ], k )


class Particle( object ):

    def __init__( self, species, serial=None, index=None ):

        self.species = species

        if not serial is None:
            self.serial = serial
        elif not index is None:
            self.serial = species.pool.getSerialByIndex( index )
        else:
            raise ValueError, 'give either serial or index.'

        self.hash = hash( self.species ) ^ self.serial


    def __str__( self ):

        return '( ' + self.species.id + ', ' + str( self.serial ) + ' )'


    def __eq__( self, other ):

        return self.species == other.species and self.serial == other.serial

    def __ne__( self, other ):

        return self.species != other.species or self.serial != other.serial

    def __cmp__( self, other ):

        if self.species == other.species:
            return self.serial - other.serial
        elif self.species < other.species:
            return -1
        else:
            return 1


    def __hash__( self ):

        return self.hash


    def getPos( self ):
        pool = self.species.pool
        return pool.positions[ pool.indexMap[ self.serial ] ]

    def setPos( self, newpos ):
        pool = self.species.pool
        pool.positions[ pool.indexMap[ self.serial ] ] = newpos

    pos = property( getPos, setPos )

    def getRadius( self ):
        return self.species.radius

    radius = property( getRadius )


    def getIndex( self ):
        return self.species.pool.indexMap[ self.serial ]


class DummyParticle( object ):
    def __init__( self ):
        self.species = None
        self.serial = -1



class ParticlePool( object ):

    def __init__( self ):

        self.indexMap = {}
        self.serialCounter = 0

        self.serials = numpy.array( [], numpy.integer )

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

        if self.size == 0:
            self.indexMap = {}
            self.__resizeArrays( self.size )
            return

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



class ParticleSimulatorBase( object ):
    
    def __init__( self ):

        self.speciesList = {}
        self.reactionTypeMap1 = {}
        self.reactionTypeMap2 = {}

        self.surfaceList = []

        #self.dt = 1e-7
        #self.t = 0.0

        self.H = 3.0
        
        self.dtLimit = 1e-3
        self.dtMax = self.dtLimit

        self.nextReaction = None

        # counters
        self.rejectedMoves = 0
        self.reactionEvents = 0

        self.particleMatrix = ObjectMatrix()


        self.setWorldSize( INF )

        self.populationChanged = False


    def initialize( self ):
        pass


    def reconstructParticleMatrix( self ):

        self.particleMatrix.clear()
        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                self.particleMatrix.add( particle,
                                         particle.pos, species.radius )



    def setWorldSize( self, size ):

        if isinstance( size, list ) or isinstance( size, tuple ):
            size = numpy.array( size )

        self.worldSize = size

        self.particleMatrix.setWorldSize( size )

        if isinstance( size, float ) and size == INF:
            self._distance = distance_Simple
            #self._distanceArray = distanceArray_Simple
            self._distanceSq = distanceSq_Simple
            self._distanceSqArray = distanceSqArray_Simple
        else:
            self._distance = distance_Cyclic
            #self._distanceArray = distanceSqArray_Cyclic
            self._distanceSq = distanceSq_Cyclic
            self._distanceSqArray = distanceSqArray_Cyclic

    def getWorldSize( self ):
        return self.worldSize

    def setMatrixSize( self, size ):
        self.particleMatrix.setMatrixSize( size )

    def applyBoundary( self, pos ):

        pos %= self.worldSize

    def getReactionType1( self, species ):
        return self.reactionTypeMap1.get( species, None )

    def getReactionType2( self, species1, species2 ):
        return self.reactionTypeMap2.get( ( species1, species2 ), None )

    def getSpeciesByIndex( self, i ):
        return self.speciesList.values()[i]


    def distanceSq( self, position1, position2 ):
        return self._distanceSq( position1, position2, self.worldSize )

    def distance( self, position1, position2 ):
        return self._distance( position1, position2, self.worldSize )
        
    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1,\
                                                 positions ) )

    def addSurface( self, surface ):
        self.surfaceList.append( surface )


    def addSpecies( self, species ):
        self.speciesList[ species.id ] = species

    def addReactionType( self, rt ):

        numReactants = len( rt.reactants )

        if numReactants == 1:
            species1 = rt.reactants[0]

            if len( rt.products ) == 1:
                if species1.radius * 2 < rt.products[0].radius:
                    raise RuntimeError,\
                        'radius of product must be smaller ' \
                        + 'than radius of reactant.'
            elif len( rt.products ) == 2:
                if species1.radius < rt.products[0].radius or\
                        species1.radius < rt.products[1].radius:
                    raise RuntimeError,\
                        'radii of both products must be smaller than ' \
                        + 'reactant.radius.'

            if self.reactionTypeMap1.has_key( species1 ):
                self.reactionTypeMap1[species1].append( rt )
            else:
                self.reactionTypeMap1[species1] = [ rt, ]

        elif numReactants == 2:
            species1 = rt.reactants[0]
            species2 = rt.reactants[1]
            self.reactionTypeMap2[ (species1,species2) ] = rt
            if species1 != species2:
                self.reactionTypeMap2[ (species2,species1) ] = rt

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
        log.info( 'throwing in %s %s particles' % ( n, species.id ) )

        for i in range( int( n ) ):

            while 1:

                #position= numpy.random.uniform( 0, self.worldSize, 3 )
                position = surface.randomPosition()
                if self.checkOverlap( position, species.radius ):
                    break
                else:
                    log.info( '%d-th particle rejected.' %i )
            
            self.createParticle( species, position )


    def placeParticle( self, species, pos ):

        pos = numpy.array( pos )
        radius = species.radius

        if not self.checkOverlap( pos, radius ):
            raise NoSpace, 'overlap check failed'
            
        particle = self.createParticle( species, pos )
        return particle


    def createParticle( self, species, pos ):
        newserial = species.newParticle( pos )
        newparticle = Particle( species, serial=newserial )
        self.addToParticleMatrix( newparticle, pos )
        return newparticle

    def removeParticle( self, particle ):
        particle.species.pool.removeBySerial( particle.serial )
        self.removeFromParticleMatrix( particle )

    def moveParticle( self, particle, newpos ):
        particle.pos = newpos
        self.updateOnParticleMatrix( particle, newpos )


    def addToParticleMatrix( self, particle, pos ):
        self.particleMatrix.add( particle,
                                 pos, particle.species.radius )

    def removeFromParticleMatrix( self, particle ):
        self.particleMatrix.remove( particle )

    def updateOnParticleMatrix( self, particle, pos ):
        self.particleMatrix.update( particle,
                                    pos, particle.species.radius )


    def checkOverlap( self, pos, radius, ignore=[] ):
        
        particles, _ = \
            self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
        if len( particles ) == 0:
            return True

        if [ p for p in particles if p not in ignore ]:
            return False
        else:
            return True


        
    def getParticlesWithinRadius( self, pos, radius, ignore=[] ):
        particles, _ =\
            self.particleMatrix.getNeighborsWithinRadius( pos, radius )

        return [ p for p in particles if p not in ignore ]


    def getParticlesWithinRadiusNoSort( self, pos, radius, ignore=[] ): 
        particles, _ =\
            self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        return [ p for p in particles if p not in ignore ]


    def clear( self ):

        self.dtMax = self.dtLimit
        self.dt = self.dtLimit

        self.nextReaction = None

        
    '''
    Get closest n Particles.

    When the optional argument speciesList is given, only Particles of
    species in the list are considered.  When speciesList is not given
    or is None, all species in the simulator are considered.
    
    This method returns a tuple ( neighbors, distances ), where neighbors
    is a list of Particle objects.
    '''


    def getNeighborParticles( self, pos, n=None ):
        n, d = self.particleMatrix.getNeighbors( pos, n )
        neighbors = [ Particle( i[0], i[1] ) for i in n ]
        return neighbors, d

    def getNeighborParticlesNoSort( self, pos, n=None ):
        n, d = self.particleMatrix.getNeighborsNoSort( pos, n )
        neighbors = [ Particle( i[0], i[1] ) for i in n ]
        return neighbors, d


    def getClosestParticle( self, pos, ignore=[] ):

        neighbors, distances =\
            self.getNeighborParticles( pos, len( ignore ) + 1 )

        for i in range( len( neighbors ) ): 
            if neighbors[i] not in ignore:
                closest, distance = neighbors[i], distances[i]

                #assert not closest in ignore
                return closest, distance

        # default case: none left.
        return None, INF
        #return DummyParticle(), INF


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
        

    def checkParticleMatrix( self ):

        if self.worldSize != self.particleMatrix.worldSize:
            raise RuntimeError,\
                'self.worldSize != self.particleMatrix.worldSize'


        total = numpy.array( [ species.pool.size for species 
                               in self.speciesList.values() ] ).sum()


        if total != self.particleMatrix.size:
            raise RuntimeError,\
                'total number of particles %d != self.particleMatrix.size %d'\
                % ( total, self.particleMatrix.size )

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                pos, radius = self.particleMatrix.get( particle )

                if ( particle.pos - pos ).sum() != 0:
                    raise RuntimeError,\
                        'particleMatrix positions consistency broken'
                if particle.species.radius != radius:
                    raise RuntimeError,\
                        'particleMatrix radii consistency broken'

    def check( self ):

        self.checkParticleMatrix()


    def dumpPopulation( self ):
        buf = ''
        for species in self.speciesList.values():
            buf += species.id + ':' + str( species.pool.size ) + '\t'

        return buf
