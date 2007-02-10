#!/usr/env python


import math
import bisect
import random
import sys
import operator

import numpy
import scipy
import scipy.optimize


from utils import *
from surface import *
import gfrdfunctions
import _gfrd


class Species:
    
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


class ReactionType:

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
        self.reactants = ( s1, )
        self.products = ( p1, )
        self.k  = k

class BindingReactionType( ReactionType ):

    def __init__( self, s1, s2, p1, k ):
        self.reactants = ( s1, s2 )
        self.products = ( p1, )
        self.k  = k

        D = s1.D + s2.D
        sigma = s1.radius + s2.radius
        self.pairGreensFunction = _gfrd.PlainPairGreensFunction( D, k, sigma )

class RepulsionReactionType( ReactionType ):

    def __init__( self, s1, s2 ):
        self.reactants = ( s1, s2 )
        self.products = ( )
        self.k  = 0.0

        D = s1.D + s2.D
        sigma = s1.radius + s2.radius
        self.pairGreensFunction = _gfrd.PlainPairGreensFunction( D, self.k,\
                                                                sigma )

class UnbindingReactionType( ReactionType ):

    def __init__( self, s1, p1, p2, k ):
        self.reactants = ( s1, )
        self.products = ( p1, p2 )
        self.k  = k
    
class Pair:

    def __init__( self, single1, single2, dt, ndt, rt=None, rdt=-1.0 ):
        self.single1 = single1
        self.single2 = single2
        self.dt = dt
        self.ndt = ndt
        self.rt = rt
        self.rdt = rdt

    def __str__( self ):
        return str( (self.rdt, self.dt, self.ndt,\
                    single1, single2, self.rt) )


class Particle:

    def __init__( self, si, i ):

        self.si = si
        self.i = i

    def __str__( self ):
        return str( ( self.si, self.i ) )





class ParticlePool:

    def __init__( self ):

        self.indexMap = {}
        self.serialCounter = 0

        self.serials = numpy.array( [0,], numpy.integer )

        self.positions = numpy.array( [], numpy.floating )
        self.positions.shape = ( 0, 3 )

        self.drs = numpy.array( [], numpy.floating )

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

        self.drs[ newindex ] = 0.0

        return newserial
    

    def __resizeArrays( self, newsize ):

        self.serials = numpy.resize( self.serials, newsize )
        self.positions = numpy.resize( self.positions, ( newsize, 3 ) )
        self.drs = numpy.resize( self.drs, newsize )


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
        self.drs[ index ] = self.drs[ self.size ]
        self.__resizeArrays( self.size )

        # book keeping
        del self.indexMap[ serialOfRemovedItem ]
        self.indexMap[ serialOfMovedItem ] = index

    def getSerialByIndex( self, index ):
        
        return self.serials[index]

    def getIndex( self, serial ):

        index = self.indexMap[ serial ]

        return index



'''
class Particle:
    
    def __init__( self, pos, species ):
        self.species = species

        pool = species.pool
        
        self.serial = pool.newLot()

        index = pool.getIndexBySerial( self.serial )


        print 'idx', index
        print pool.positions
        pool.positions[ index ] = pos

        self.partner = None

        # temporary area, should be gotten rid of
        self.dt = array( (scipy.Inf,scipy.Inf) )

    def getPos( self ):
        pool = self.species.pool
        return pool.positions[ pool.getIndexBySerial( self.serial ) ]

    def setPos( self ):
        pool = self.species.pool
        return pool.positions[ pool.getIndexBySerial( self.serial ) ]
'''        



class GFRDSimulatorBase:
    
    def __init__( self ):
        self.speciesList = {}
        self.reactionTypeList1 = {}
        self.reactionTypeList2 = {}

        self.surfaceList = []

        self.dt = 1e-7
        self.t = 0.0

        self.H = 3.0
        
        self.dtLimit = 1e-3
        self.dtMax = self.dtLimit

        self.nextReaction = None

        self.fsize = 0.0


        # counters
        self.rejectedMoves = 0
        self.reactionEvents = 0


        # internal variables
        #self._distanceSq = distanceSq_Simple
        #self._distanceSqArray = distanceSqArray_Simple
        self._distanceSq = distanceSq_Cyclic
        self._distanceSqArray = distanceSqArray_Cyclic


    def getReactionType2( self, species1, species2 ):
        return self.reactionTypeList2.get( ( species1, species1 ), None )

    def getSpeciesByIndex( self, i ):
        return self.speciesList.values()[i]

    def simpleDiffusion( self, speciesIndex, particleIndex ):

        species = self.speciesList.values()[speciesIndex]
        if species.D == 0.0:
            return

        limitSq = self.H * self.H * ( 6.0 * species.D * self.dt )

        while True:
            displacement = gfrdfunctions.p1( species.D, self.dt )
            
            distanceSq = ( displacement * displacement ).sum()

            if distanceSq <= limitSq:
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
            displacement = gfrdfunctions.p1( species.D, self.dt )
            
            distanceSq = ( displacement * displacement ).sum()

            if distanceSq <= limitSq:
                break

            self.rejectedMoves += 1

        pos = species.pool.positions[particleIndex]
        pos += displacement

        #SURFACE

        


    def distanceSq( self, position1, position2 ):
        return self._distanceSq( position1, position2, self.fsize )

    def distance( self, position1, position2 ):
        return math.sqrt( self.distanceSq( position1, position2 ) )
        
    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.fsize )

    def setSize( self, size ):
        self.fsize = size

    def addSurface( self, surface ):
        self.surfaceList.append( surface )


    def addSpecies( self, species ):
        self.speciesList[ species.id ] = species

    def addReactionType( self, r ):

        numReactants = len( r.reactants )


        if numReactants == 1:
            species1 = r.reactants[0]
            self.reactionTypeList1[species1] = r
        elif numReactants == 2:
            species1 = r.reactants[0]
            species2 = r.reactants[1]
            self.reactionTypeList2[ (species1,species2) ] = r
            if species1 != species2:
                self.reactionTypeList2[ (species2,species1) ] = r
        else:
            raise 'unexpected'

    def setAllRepulsive( self ):
        for species1 in self.speciesList.values():
            for species2 in self.speciesList.values():
                try:
                    rt = self.reactionTypeList2[ (species1,species2) ]
                except:
                    self.reactionTypeList2[ (species1,species2) ] =\
                                            RepulsionReactionType( species1,\
                                                                   species2 )
        

    def throwInParticles( self, id, n, surface ):
        print 'throwing in %s %s particles' % ( n, id )

        species = self.speciesList[ id ]
        
        for i in range( n ):

            while True:

                #position= numpy.random.uniform( 0, self.fsize, 3 )
                position = surface.randomPosition()
                if self.checkOverlap( position, species.radius ):
                    break
            
            species.newParticle( position )



    def placeParticle( self, id, position ):

        position = numpy.array( position )
        species = self.speciesList[ id ]

        if not self.checkOverlap( position, species.radius ):
            raise 'placeParticle: overlap check failed'
            
        species.newParticle( position )
        
        
    def checkOverlap( self, position, radius ):
        
        for species2 in self.speciesList.values():

            if species2.pool.size == 0:
                continue

            positions2 = species2.pool.positions

            radius12 = radius + species2.radius
            radius12sq = radius12 * radius12

            closestdistsq = self.distanceSqArray( position, positions2 ).min()

            if closestdistsq <= radius12sq:
                print 'reject:', math.sqrt(closestdistsq), math.sqrt( radius12sq )
                return False

        return True
    

    def clear( self ):

        self.dtMax = self.dtLimit
        self.dt = self.dtLimit

        self.nextReaction = None

        self.pairs = []
        self.singles = []



    def isPopulationChanged( self ):
        return self.nextReaction != None


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

    def formPairs( self ):

        # 1. form pairs in self.pairs
        # 2. list singles in self.singles

        speciesList = self.speciesList.values()

        # list up pair candidates

        # partner -> nearest particle
        # neighbor -> second nearest particle

        # dtCache[ speciesIndex ][ particleIndex ][ 0 .. 1 ]
        dtCache = []
        neighborCache = []
        checklist = []

        for speciesIndex in range( len( speciesList ) ):
            size = speciesList[speciesIndex].pool.size

            dtCache.append( numpy.zeros( ( size, 2 ), numpy.floating ) )
            neighborCache.append( [[[ -1, -1 ],[-1,-1]]] * size )

            checklist.append( numpy.ones( speciesList[speciesIndex].pool.size ) )
            for particleIndex in range( size ):

                dt, neighbor = self.checkPairs( speciesIndex, particleIndex )

                dtCache[ speciesIndex ][ particleIndex ] = dt
                neighborCache[ speciesIndex ][ particleIndex ] = neighbor

        self.pairs = []
        for speciesIndex1 in range( len( speciesList ) ):

            species1 = speciesList[speciesIndex1]

            for particleIndex1 in range( species1.pool.size ):

                #bdt, surface = self.checkSurfaces( speciesIndex,\
                #particleIndex )

                # skip if this particle has already taken in a pair.
                if checklist[speciesIndex1][particleIndex1] == 0:
                    #print 'skip', speciesIndex1, particleIndex1
                    continue

                # A partner: the other of the pair.
                # A neighbor of a pair: closer of the second closest of
                #                       the particles in the pair.
                #                       This is different from neighbors of
                #                       a particle.
                
                # (1) Find the closest particle (partner).
                partner = neighborCache[ speciesIndex1 ][ particleIndex1 ][0]

                ( speciesIndex2, particleIndex2 ) = partner

                if speciesIndex2 == -1:
                    continue

                dts = dtCache[ speciesIndex1 ][ particleIndex1 ]

                partnersPartner = neighborCache\
                                  [ speciesIndex2 ][ particleIndex2 ][0]
                partnerDts = dtCache[ speciesIndex2 ][ particleIndex2 ]

                # (2) The partner's partner has to be this, otherwise
                #     this combination isn't a pair.
                # (3) 'Neighbor' of this pair is the closer of
                #     this and the partner's second closest.
                #     We take the particle that has this neighbor.
                if partnersPartner != ( speciesIndex1, particleIndex1 ) or \
                       partnerDts[1] < dts[1]:
                    continue
                
                # (4) Now we have a candidate pair.
                species2 = speciesList[speciesIndex2]
                rt = self.reactionTypeList2.get( ( species1, species2 ) )

                #pair = ( dts[0], dts[1], speciesIndex1, particleIndex1,\
                #speciesIndex2, particleIndex2, rt )
                pair = Pair( dts[0], dts[1], speciesIndex1, particleIndex1,\
                             speciesIndex2, particleIndex2, rt )
                self.pairs.append( pair )

                # (5) dtMax = the minimum neighbor dt of all pairs.
                self.dtMax = min( self.dtMax, dts[1] )

                # (6) book keeping
                checklist[speciesIndex1][particleIndex1] = 0
                checklist[speciesIndex2][particleIndex2] = 0



        # screening pairs
        self.pairs.sort(key=operator.attrgetter('dt'))

        checklist = []
        for i in range( len( speciesList ) ):
            checklist.append( numpy.ones( speciesList[i].pool.size ) )

        for i in range( len( self.pairs ) ):
            #( dt, ndt, si1, i1, si2, i2, rt ) = self.pairs[i]
            pair = self.pairs[i]


            # Don't take pairs with partner dt greater than dtMax.
            if pair.dt > self.dtMax:
                self.pairs = self.pairs[:i]
                break   # pairs are sorted by dt.  break here.

            if checklist[pair.si1][pair.i1] == 0 or \
                   checklist[pair.si2][pair.i2] == 0:
                print self.pairs[:i+1]
                print dtCache[pair.si1][pair.i1], dtCache[pair.si2][pair.i2]
                print neighborCache[pair.si1][pair.i1], \
                      neighborCache[pair.si2][pair.i2]
                print 'pairs not mutually exclusive.'
                self.pairs = self.pairs[:i]
                break

            checklist[pair.si1][pair.i1] = 0
            checklist[pair.si2][pair.i2] = 0


        # now we have the final list of pairs

        # next, make the list of singles.
        # a single is a particle that doesn't appear in the list
        # of pairs.
        self.singles = []

        for i in range( len( checklist ) ):
            singleIndices = numpy.nonzero( checklist[i] )[0]
            for j in singleIndices:
                self.singles.append( Single( INF, i, j ) )
        #    singleIndices = numpy.nonzero( checklist[i] )[0] #== flatnonzero()
        #    self.singles.append( singleIndices )

        #debug
        numSingles = len( self.singles )

        print '# pairs = ', len(self.pairs),\
              ', # singles = ', numSingles



    def checkPairs( self, speciesIndex1, particleIndex ):

        speciesList = self.speciesList.values()

        neighbordts = [ INF, INF ]
        #           [ ( index, reaction type ), (...,...), .. ]
        neighbors = [ ( -1, -1 ), ( -1, -1 ) ]
        drSqs = []

        species1 = speciesList[ speciesIndex1 ]
        positions = species1.pool.positions
        position1 = positions[ particleIndex ].copy()

        if self.reactionTypeList2.get( ( species1, species1 ), None ) != None \
           and len( position1 ) >= 2 and species1.D != 0.0:

            # temporarily displace the particle
            positions[particleIndex] = NOWHERE
            
            topDts, topIndices = self.checkDistance( position1, positions,
                                                     species1, species1 )

            # restore the particle.
            positions[particleIndex] = position1

            neighbordts.extend( topDts )

            if len( topIndices ) == 2:
                neighbors.extend( ( ( speciesIndex1, topIndices[0] ),\
                                    ( speciesIndex1, topIndices[1] ) ) )
            else:
                neighbors.extend( ( ( speciesIndex1, topIndices[0] ), ) )

        #for speciesIndex2 in range( speciesIndex1 + 1, len( speciesList ) ):
        for speciesIndex2 in range( speciesIndex1 )\
                + range( speciesIndex1 + 1, len( speciesList ) ):
            species2 = speciesList[speciesIndex2]

            # non reactive
            if self.reactionTypeList2.get( ( species1, species2 ), None )\
                   == None:
                continue
            
            if species2.pool.size == 0:
                continue

            if species1.D + species2.D == 0.0:
                continue
                    
            positions = species2.pool.positions

            if species2.pool.size == 1:  # insert a dummy
                positions = numpy.concatenate( ( positions, [NOWHERE,] ) )
                
            topDts, topIndices = self.checkDistance( position1, positions,
                                                     species1, species2 )
            neighbordts.extend( topDts )
            neighbors.extend( ( ( speciesIndex2, topIndices[0] ),\
                                ( speciesIndex2, topIndices[1] ) ) )

        topargs = numpy.argsort( neighbordts )[:2]
        topNeighborDts = numpy.take( neighbordts, topargs )
        topNeighbors = ( neighbors[topargs[0]], neighbors[topargs[1]] )

        return topNeighborDts, topNeighbors


    def checkDistance( self, position1, positions2, species1, species2 ):

        #positions2 = species2.pool.positions

        distanceSq = self.distanceSqArray( position1, positions2 )
        sortedindices = distanceSq.argsort()

        #print sortedindices
        #debug
        radius12 = species1.radius + species2.radius
        radius12sq = radius12 * radius12

        # check if particles overlap
        if distanceSq[ sortedindices[0] ] < radius12sq - 1e-20 and \
               distanceSq[ sortedindices[0] ] != 0.0:
            print position1, positions2[ sortedindices[0] ]
            print 'dr<radius', math.sqrt(distanceSq[sortedindices[0]]), radius12
            print species1.id, species2.id, sortedindices[0]
            raise "critical"

        factor = ( math.sqrt( species1.D ) + math.sqrt( species2.D ) ) * self.H
        factor *= factor
        factor *= 6.0

        distanceSqSorted = distanceSq.take( sortedindices[:2] )

        # instead of just
        #dts = distanceSqSorted / factor
        # below takes into account of particle radii.
        dts = numpy.sqrt( distanceSqSorted ) - radius12
        dts *= dts
        dts /= factor

        #if min( dts ) < 0.0:
        #print 'negative dts occured, clipping.', dts
        #dts = numpy.clip( dts, 0.0, INF )
        # raise 'stop'

        indices = sortedindices[:2]

        return dts, indices # , distanceSqSorted[0]


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
        
