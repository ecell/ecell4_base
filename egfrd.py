#!/usr/env python


import math

import numpy
#import scipy
#import scipy.optimize


from utils import *
from surface import *

from gfrdbase import *




"""
class DistanceMatrix:

    def __init__( self, numSpecies ):
        row = [ numpy.array([], numpy.floating ), ] * numSpecies
        self.matrix = [ row, ] * numSpecies

    def __getitem__( self, si1, i1, si2, i2 ):
        return self.matrix[si1][si2][
"""



class Single:

    def __init__( self, sim, particle ):

        self.particle = particle
        self.sim = sim
        self.lastTime = 0.0
        self.dt = 0.0
        self.closest = None
        self.eventID = None

        self.gf = FirstPassageGreensFunction( particle.species.D )


    '''
    Initialize this Single.

    The protective sphere size, self.lastTime, self.dt, and self.closest
    are updated.  The protective sphere size is determined without
    taken into account of other shells.

    '''

    def initialize( self ):

        neighbors, drs = self.sim.getNeighbors( self.particle.getPos() )
        closest = neighbors[1]
        dr = drs[1] - self.particle.species.radius
        dr *= .5
        closestParticle = Particle( closest[0], index=closest[1] )
        closestSingle = self.sim.findSingle( closestParticle )
        
        self.closest = closestSingle
        
        #print 'dr', dr
        #FIXME: take different D into account
        self.setDr( dr )
        
        self.lastTime = self.sim.t
        self.dt = self.calculateFirstPassageTime()
        


    def fire( self ):

        rnd = numpy.random.uniform( size=2 )

        displacementS = [ self.getDr(), rnd[0] * Pi, rnd[1] * 2 * Pi ]
        displacement = sphericalToCartesian( displacementS )

        #self.checkShellForAll()
        
        pos = self.particle.getPos()
        pos += displacement

        # BOUNDARY
        self.sim.applyBoundary( pos )
        
        #print 'displacement', length(displacement), self.getDr()

        neighbors, distances = self.sim.getNeighborShells( pos )
        newdr = distances[1] - self.particle.species.radius
        
        self.closest = self.sim.findSingle( Particle( neighbors[1][0],
                                                    neighbors[1][1] ) )

        #print 'newdr', newdr
        if newdr <= 0:
            print newdr, self.closest
            raise RuntimeError, 'Fatal newdr <= 0'

        self.setDr( newdr )

        #print self, self.closest
        #debug
        #self.checkShellForAll()

        self.dt = self.calculateFirstPassageTime()

        return self.dt


    '''
    Update the position of the particle and the protective sphere.
    '''
    
    def update( self, t ):

        assert t > self.lastTime
        assert self.getDr() > 0.0

        rnd = numpy.random.uniform( size=3 )

        dt = t - self.lastTime
        r = self.gf.drawR( rnd[0], dt, self.getDr() )
        theta = rnd[1] * Pi
        phi = rnd[2] * 2 * Pi
        displacement = sphericalToCartesian( [ r, theta, phi ] )
        newPos = self.particle.getPos() + displacement

        self.particle.setPos( newPos )

        self.updateShell()
        self.updateDt()

        self.lastTime = t

    '''
    Update the protective sphere.
    self.closest is updated too.
    '''

    def updateShell( self ):

        neighbors, drs = self.sim.getNeighborShells( self.particle.getPos() )
        closest = neighbors[1]
        dr = drs[1] - self.particle.species.radius
        print 'dr dr1', dr, drs[1]
        closestParticle = Particle( closest[0], index=closest[1] )
        closestSingle = self.sim.findSingle( closestParticle )
        
        self.closest = closestSingle
        
        #print 'dr', dr
        #FIXME: take different D into account
        self.setDr( dr )

    def updateDt( self ):
        self.dt = self.calculateFirstPassageTime()        
        

    def isDependentOn( self, event ):
        #print event
        return False

    def setPosition( self, pos ):
        self.particle.species.pool.positions[self.particle.i] = pos

    def getPosition( self ):
        return self.particle.species.pool.positions[self.particle.i]

    def setDr( self, dr ):
        pool = self.particle.species.pool
        pool.drs[ pool.getIndex( self.particle.serial ) ] = dr

    def getDr( self ):
        pool = self.particle.species.pool
        return pool.drs[ pool.getIndex( self.particle.serial ) ]

    def calculateFirstPassageTime( self ):
        
        species = self.particle.species
        rnd = numpy.random.uniform()
        dr = self.getDr()
        if dr <= 0.0:
            raise RuntimeError, 'dr <= 0.0: %s' % str(dr)
        dt = self.gf.drawTime( rnd, dr )
        print dt
        if dt <= 0.0:
            raise RuntimeError, 'dt <= 0.0: %s' % str(dt)
        return dt


    def __str__( self ):
        return str( self.particle )



class Pair:
    
    def __init__( self, sim, single1, single2, rt ):

        self.single1 = single1
        self.single2 = single2

        self.rt = rt
        
        self.sim = sim
        self.lastTime = 0.0
        self.dt = 0.0
        self.closest = None

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        D1, D2 = particle1.species.D, particle2.species.D
        D12 = D1 + D2
        self.sqrtD1D2 = math.sqrt( D1 / D2 )
        
        self.sigma = particle1.species.radius + particle2.species.radius

        self.sgf = FirstPassageGreensFunction( D12 / 4.0 )
        self.pgf = FirstPassagePairGreensFunction( D12, rt.k, self.sigma )

        self.eventID = None



    def __del__( self ):
        print 'del', str( self )



    '''
    Calculate and return the "Center of Mass" (== CoM) of this pair.
    '''

    def getCoM( self ):
        particle1 = self.single1.particle
        particle2 = self.single2.particle
        
        pos1 = particle1.getPos()
        pos2 = particle2.getPos()
        
        #com = sqrtD2D1 * pos1 + self.sqrtD1D2 * pos2
        com = ( pos1 + self.sqrtD1D2 * pos2 ) * .5

        return com

    '''
    Calculate new positions of the pair particles using
    a new center-of-mass, a new inter-particle vector, and
    an old inter-particle vector.

    '''

    def newPositions( self, newCoM, newInterParticle, oldInterParticle ):

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotationAxis = crossproduct( [ 0,0,1 ], interParticle )
        
        rotationAxis = crossproductAgainstZAxis( oldInterParticle )
        rotationAxis = normalize( rotationAxis )
        
        angle = vectorAngleAgainstZAxis( oldInterParticle )
        
        newInterParticle = rotateVector( newInterParticle,
                                         rotationAxis,
                                         angle )
        
        
        newpos1 = ( 2 * newCoM - self.sqrtD1D2 * newInterParticle ) \
                  / ( 1 + self.sqrtD1D2 )
        newpos2 = newpos1 + newInterParticle

        return newpos1, newpos2
        

    def nextEvent( self, dr ):

        self.dr = dr
        
        rnd = numpy.random.uniform( size=3 )

        pos1 = self.single1.particle.getPos()
        pos2 = self.single2.particle.getPos()

        self.r0 = self.sim.distance( pos1, pos2 )

        self.a_r = ( self.dr + self.r0 ) * .5
        self.a_R = self.a_r - self.r0

        #print 'dr', dr, 'r0', r0, 'a_r', a_r, 'a_R', a_R, dr - a_r - a_R

        self.pgf.seta( self.a_r )

        self.t_R = self.sgf.drawTime( rnd[0], self.a_R )
        self.t_r = self.pgf.drawTime( rnd[1], self.r0 )

        if self.t_R < self.t_r:
            self.dt = self.t_R
            self.eventType = 2
        else:
            self.dt = self.t_r
            self.eventType = self.pgf.drawEventType( rnd[2], self.r0,
                                                     self.t_r )

        return self.dt, self.eventType


    def fire( self ):

        print 'fire:', self


        particle1 = self.single1.particle
        particle2 = self.single2.particle
        species1 = particle1.species
        species2 = particle2.species

        pos1 = particle1.getPos()
        pos2 = particle2.getPos()

        D1 = species1.D
        D2 = species2.D

        oldInterParticle = pos2 - pos1

        # Three cases:
        #  1. Reaction
        #  2.1 Escaping through a_r.
        #  2.2 Escaping through a_R.

        # 1. Reaction
        if self.eventType == EventType.REACTION:

            print 'reaction'

            if len( self.rt.products ) == 1:
                
                species3 = self.rt.products[0]

                if D1 == 0.0:
                    newR = pos1
                elif D2 == 0.0:
                    newR = pos2
                else:
                    R0 = self.getCoM()
                    dR = p_free( ( D1 + D2 ) / 4, self.dt )
                    newR = R0 + dR
                
                
                #FIXME: SURFACE
                newPos = self.sim.applyBoundary( newR )

                self.sim.removeParticle( particle1 )
                self.sim.removeParticle( particle2 )

#                species1.removeParticleBySerial( particle1.serial )
#                species2.removeParticleBySerial( particle2.serial )

                particle = self.sim.createParticle( species3, newPos )
                self.sim.insertParticle( particle )
                return -1

            else:
                raise NotImplementedError,\
                      'num products >= 2 not supported yet.'

        # 2.1 Escaping through a_r.
        elif self.eventType == EventType.ESCAPE:

            print 'escape r'

            rnd = numpy.random.uniform( size=5 )

            # calculate new R
            
            r_R = self.sgf.drawR( rnd[0], self.dt, self.a_R )
            
            displacement_R_S = [ r_R,
                                 rnd[1] * Pi,
                                 rnd[2] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = self.getCoM() + displacement_R


            # calculate new r
            print ( rnd[3], self.a_r, self.r0, self.dt )
            theta_r = self.pgf.drawTheta( rnd[3], self.a_r*.09, self.r0, self.dt )
            phi_r = rnd[4] * 2 * Pi
            newInterParticleS = numpy.array( [ self.a_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = self.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )


        # 2.2 escaping through a_R.
        elif self.eventType == 2:

            print 'escape R'

            rnd = numpy.random.uniform( size = 5 )

            # calculate new r
            r = self.pgf.drawR( rnd[0], self.a_r, self.dt )
            print ( rnd[1], r, self.r0, self.dt )
            theta_r = self.pgf.drawTheta( rnd[1], r, self.r0, self.dt )
            phi_r = rnd[2] * 2*Pi
            newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            # calculate new R
            displacement_R_S = [ self.a_R,
                                 rnd[3] * Pi,
                                 rnd[4] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = self.getCoM() + displacement_R


            newpos1, newpos2 = self.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
                
        else:
            raise SystemError, 'Bug: invalid eventType.'

        particle1.setPos( newpos1 )
        particle2.setPos( newpos2 )

        # check if the new positions are valid:
        # FIXME: check if positions are inside the shells too.
        newParticleDistance = distance( newpos1, newpos2 )
        if newParticleDistance <= species1.radius + species2.radius:
            print 'rejected move: ', 'radii, interp',\
                  species1.radius + species2.radius, newParticleDistance
            print 'DEBUG: r0, dt, pos1, pos2, newpos1, newpos2',\
                  self.r0, self.dt, pos1, pos2, newpos1, newpos2
            raise RuntimeError

        # here decide whether this pair still continues or breaks up

        # temporarily, it alwayw breaks up to singles.
        single1 = self.sim.findSingle( particle1 )
        single2 = self.sim.findSingle( particle2 )

        # protect the singles with shells.
        single1.setDr( single1.particle.species.radius )
        single2.setDr( single2.particle.species.radius )
        single1.updateShell()
        single2.updateShell()

#         if single1.closest == single2:
#             single1.setDr( newParticleDistance * .49 )
#         if single2.closest == single1:
#             single2.setDr( newParticleDistance * .49 )
        print 'aa', newParticleDistance
        single1.updateDt()
        single2.updateDt()

        self.sim.checkShell( single1 )
        self.sim.checkShell( single2 )

        self.sim.createSingleEvent( single1 )
        self.sim.createSingleEvent( single2 )

        return -1


    def update( self, t ):
        print 'update ', t

    def isDependentOn( self, event ):
        #print event
        return False

    def __str__( self ):
        return str(self.single1.particle) + ' ' + str(self.single2.particle)


class EGFRDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        GFRDSimulatorBase.__init__( self )

        self.isDirty = True

        self.scheduler = EventScheduler()

        self.t = 0.0
        self.dtMax = INF
        self.dt = INF

        #self.pairList = []
        self.singleMap = {}

        self.lastEvent = None
        self.hoggerCounter = 0

    def initialize( self ):

        self.scheduler.clear()

        self.initializeSingleMap()
        self.initializeSingles()

        self.formPairs()

        #debug
        #self.checkShellForAll()

        #self.scheduler.updateAllEventDependency()

        self.isDirty = False



    def step( self ):

        self.checkInvariants()
        if self.isDirty:
            self.initialize()

        self.t, self.lastEvent = self.scheduler.getTopEvent()

        self.scheduler.step()

        nextTime, nextEvent = self.scheduler.getTopEvent()
        self.dt = nextTime - self.t
        
        assert self.scheduler.getSize() != 0


        # if the same event stepped in the last n steps,
        # reinitialize everything.
        # FIXME: don't need to initialize everything.
        #        (1) recalculate dr to the closest with
        #            its dr shrunken.
        #        (2) new dr is
        #            min( dr to the second closest,
        #                 dr to the closest, with its dr
        #                 shrunken )

        if self.lastEvent is nextEvent:
            self.hoggerCounter += 1
        else:
            self.hoggerCounter = 0

        if self.hoggerCounter >= 10: # or self.dt < 1e-15:
                print 'reinitialize'
                self.hoggerCounter = 0
                self.initialize()
                #nextEvent = self.scheduler.getTopEvent()[1]
                self.dt = self.scheduler.getTime() - self.t

        #if self.dt == 0.0:
        #    raise 'dt=0'


        print 'dt', self.dt,\
              'reactions', self.reactionEvents,\
              'rejected moves', self.rejectedMoves,\
              'hogger counter', self.hoggerCounter
        print ''
        


    def findSingle( self, particle ):
        return self.singleMap.get( ( particle.species, particle.serial ) )

    def createSingle( self, particle ):
        single = Single( self, particle )
        self.singleMap[ ( particle.species, particle.serial ) ] = single
        return single

    def removeSingle( self, single ):
        particle = single.particle
        del self.singleMap[ ( particle.species, particle.serial ) ]
        return single


    def removeParticle( self, particle ):
        single = self.findSingle( particle )
        self.removeSingle( single )
        particle.species.removeParticleBySerial( particle.serial )

    def initializeSingleMap( self ):

        self.singleMap = {}

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
                single.setDr( 0.0 )

    def initializeSingles( self ):

        for single in self.singleMap.values():
            single.initialize()
            nextt = single.lastTime + single.dt
            self.scheduler.addEvent( nextt, single )

            #self.createSingleEvent( single )
        
        #FIXME: here perhaps a second pass to get drs maximally large.


    def createSingleEvent( self, single ):
        print 'create', self.t+single.dt, self.t, single.dt, single
        self.scheduler.addEvent( self.t + single.dt, single )

    def insertParticle( self, particle ):

        single = self.createSingle( particle )
        single.updateShell()
        single.updateDt()
        self.createSingleEvent( single )


    def formPairs( self ):

        for single in self.singleMap.values():

            #neighbors, drs = self.getNeighbors( single.particle.getPos() )
            #closest = Particle( neighbors[1][0], neighbors[1][1] )
            closest = single.closest

            if single.particle < closest.particle:
                continue

            partnerNeighbors, partnerDrs = self.getNeighbors( closest.particle.getPos())
            partnerClosest = Particle( partnerNeighbors[1][0],
                                       index=partnerNeighbors[1][1] )

            if single.particle == partnerClosest:
                
                pos1 = single.particle.getPos()
                pos2 = closest.particle.getPos()
                r0 = self.distance( pos1, pos2 )

                pair = self.createPair( single, closest )
                com = pair.getCoM()

                neighbors, drs = self.getNeighbors( com, 3 )
                pairDr = drs[2]

                if pairDr < r0:   # this happens with a small probability
                    print 'pairDr < r0', pairDr, r0, pairDr - r0
                    #raise ''
                    break
                
                if pairDr > r0 * 5:
                    pairDr = r0 * 5

                nextEvent = pair.nextEvent( pairDr )
                dt = nextEvent[0]
                #eventType = nextEvent[1]

                if pair.single1.eventID != None:
                    self.scheduler.removeEvent( pair.single1.eventID )
                if pair.single2.eventID != None:
                    self.scheduler.removeEvent( pair.single2.eventID )

                pair.eventID = self.scheduler.addEvent( self.t + dt, pair ) 
                

    def createPair( self, single1, single2 ):

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.reactionTypeMap2.get( ( species1, species2 ) )

        # If either one of D1 or D2 are zero, it must be D1.
        # D1 and D2 must not be zero at the same time.
        if single2.particle.species.D == 0.0:
            if single1.particle.species.D == 0.0:
                raise RuntimeError, 'createPair: D1 == D2 == 0.'
            else:
                s1, s2 = single2, single1 # D1 must be nonzero.
        else:
            s1, s2 = single1, single2

        
        return Pair( self, s1, s2, rt )

            
        
    def checkShell( self, single ):
        neighbors, drs = self.getNeighborShells( single.particle.getPos() )
        closest, distance = neighbors[1], drs[1]
        distance -= single.particle.species.radius
        if single.getDr() - distance >= 1e-18:
            dr = single.getDr()
            print single.particle, closest, dr, distance, dr - distance
            raise RuntimeError, 'Fatal: shells overlap.'

    def checkShellForAll( self ):
        for single in self.singleMap.values():
            self.checkShell( single )



    def getNeighbors( self, pos, n=2, speciesList=None ):

        topNeighbors = []
        topDistances = []

        if speciesList == None:
            speciesList = self.speciesList.values()

        for species in speciesList:

            # empty
            if species.pool.size == 0:
                continue

            positions = species.pool.positions

            indices, distances = self.getNeighborsInSpecies( pos,
                                                             positions,
                                                             n )
            distances -= species.radius

            topDistances.extend( distances )
            topNeighbors.extend( [ ( species, i ) for i in indices ] )

        topargs = numpy.argsort( topDistances )[:n]
        topDistances = numpy.take( topDistances, topargs )
        topNeighbors = [ topNeighbors[arg] for arg in topargs ]

        return topNeighbors, topDistances


    def getNeighborsInSpecies( self, position1, positions, n=2 ):

        distances = self.distanceSqArray( position1, positions )
        
        indices = distances.argsort()[:n]
        distances = distances.take( indices )
        distances = numpy.sqrt( distances )
        
        return indices, distances


    def getNeighborShells( self, pos, n=2, speciesList=None ):

        topNeighbors = []
        topDistances = []

        if speciesList == None:
            speciesList = self.speciesList.values()

        for species in speciesList:

            # empty
            if species.pool.size == 0:
                continue

            positions = species.pool.positions
            indices, distances = self.getNeighborShellsInSpecies( pos,
                                                                  species,
                                                                  positions,
                                                                  n )
            topDistances.extend( distances )
            topNeighbors.extend( [ ( species, i ) for i in indices ] )

        topargs = numpy.argsort( topDistances )[:n]
        topDistances = numpy.take( topDistances, topargs )
        topNeighbors = [ topNeighbors[arg] for arg in topargs ]

        return topNeighbors, topDistances


    def getNeighborShellsInSpecies( self, position1, species, positions, n=2 ):

        distances = self.distanceSqArray( position1, positions )

        distances = numpy.sqrt( distances )
        distances -= species.pool.drs

        indices = distances.argsort()[:n]
        distances = distances.take( indices )
        return indices, distances

    def checkInvariants( self ):

        assert self.t >= 0.0
        assert self.dt >= 0.0
        
        self.checkShellForAll()
        
        
