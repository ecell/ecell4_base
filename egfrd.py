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

    def setPos( self, pos ):
        self.particle.setPos( pos )

    def getPos( self ):
        return self.particle.getPos()

    def setShellSize( self, size ):

        if size < self.getRadius():
            raise RuntimeError, 'shell size < radius; %g %g %g' % \
                  ( size, self.getRadius(), size - self.getRadius() )

        pool = self.particle.species.pool
        pool.drs[ pool.getIndex( self.particle.serial ) ] = size


    '''
    A shell size of a particle is the distance from the current position
    of the particle to the farthest point in space that it can occupy
    when it made the maximum displacement defined by the mobility radius
    of the particle.
    '''

    def getShellSize( self ):
        pool = self.particle.species.pool
        return pool.drs[ pool.getIndex( self.particle.serial ) ]

    def getRadius( self ):
        return self.particle.species.radius

    '''
    A mobility radius indicates the maximum displacement this single
    particle can make.

    Mobility radius of a particle is calculated as follows;

    mobility radius = shell size - radius.

    '''
    
    def getMobilityRadius( self ):
        return self.getShellSize() - self.getRadius()


    '''
    Initialize this Single.

    The shell size is shrunken to the particle radius.
    self.lastTime is reset to the current time, and self.dt
    is set to zero.

    '''

    def initialize( self ):

        neighbors, distances = self.sim.getNeighbors( self.particle.getPos() )
        closest = neighbors[1]
        closestParticle = Particle( closest[0], index=closest[1] )
        closestSingle = self.sim.findSingle( closestParticle )
        
        self.closest = closestSingle

        self.setShellSize( self.getRadius() )
        self.lastTime = self.sim.t
        self.dt = 0.0


    def reinitialize( self ):
        
        self.update( self.sim.t )


    def fire( self ):

        #debug
        self.sim.checkShellForAll()

        # (1) propagate

        rnd = numpy.random.uniform( size=2 )

        radius = self.getRadius()

        r = self.getMobilityRadius()
        displacementS = [ r, rnd[0] * Pi, rnd[1] * 2 * Pi ]
        displacement = sphericalToCartesian( displacementS )

        pos = self.particle.getPos()
        self.particle.setPos( pos + displacement )

        # BOUNDARY
        self.sim.applyBoundary( pos )

        self.lastTime = self.sim.t


        # (2) pair check
        neighborShells, distances = self.sim.getNeighbors( pos )

        # (3) determine new shell size and dt.

        neighborShells, distances = self.sim.getNeighborShells( pos )
        self.closest = neighborShells[1]
        distanceToClosestShell = distances[1]
        #print neighborShells, distances

        shellSize = self.getShellSize()

        ShellSizeDisparityFactor = 2

        closestMobilityRadius = self.closest.getMobilityRadius()
        shellSize = min( closestMobilityRadius * ShellSizeDisparityFactor
                         + ( distanceToClosestShell - radius ) * 0.5 + radius,
                         distanceToClosestShell )

        shellSize = shellSize * ( 1.0 - 1e-8 ) # safety
        shellSize = max( shellSize, radius )

        assert shellSize <= distanceToClosestShell

        self.setShellSize( shellSize )

        self.updateDt()

        return self.dt


    '''
    Update the position of the particle.

    Shell size shrunken to the radius.   self.lastTime is reset.
    self.dt is set to 0.0.
    '''
    
    def update( self, t ):

        assert t >= self.lastTime
        assert self.getShellSize() >= self.getRadius()

        if t == self.lastTime or self.getMobilityRadius() == 0.0:
            return

        rnd = numpy.random.uniform( size=3 )

        dt = t - self.lastTime
        #print rnd[0], dt, self.getMobilityRadius()
        r = self.gf.drawR( rnd[0], dt, self.getMobilityRadius() )
        theta = rnd[1] * Pi
        phi = rnd[2] * 2 * Pi
        displacement = sphericalToCartesian( [ r, theta, phi ] )
        newPos = self.particle.getPos() + displacement

        self.particle.setPos( newPos )

        self.lastTime = t

        self.setShellSize( self.getRadius() )
        self.dt = 0.0


    def getNeighborShells( self, n=2 ):
        return self.sim.getNeighborShells( self.particle.getPos(), n )

#     '''
#     Update the protective sphere.
#     self.closest is updated too.
#     '''

#     def updateShell( self ):

#         neighbors, distancess = self.getNeighborShells()
#         self.closest = neighbors[1]
#         distance = distances[1]
#         print 'distance', distance
        
#         #FIXME: take different D into account
#         self.setShellSize( distance )

    def updateDt( self ):
        self.dt = self.calculateFirstPassageTime()        
        

    def isDependentOn( self, event ):
        #print event
        return False

    def calculateFirstPassageTime( self ):
        
        rnd = numpy.random.uniform()
        r = self.getMobilityRadius()
        dt = self.gf.drawTime( rnd, r )
        if dt < 0.0:
            raise RuntimeError, 'dt < 0.0: %s' % str(dt)
        return dt


    def __str__( self ):
        return 'Single: ' + str( self.particle )



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

        self.radius = max( self.single1.particle.species.radius * 2,
                           self.single2.particle.species.radius )



    def __del__( self ):
        print 'del', str( self )

    def getPos( self ):
        return self.getCoM()

    def getShellSize( self ):
        return self.a_r + self.a_R

    def getRadius( self ):
        return self.radius

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

        oldCoM = self.getCoM()


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
                    dR = p_free( ( D1 + D2 ) / 4, self.dt )
                    newR = oldCoM + dR
                
                
                #FIXME: SURFACE
                newPos = self.sim.applyBoundary( newR )

                self.sim.removeParticle( particle1 )
                self.sim.removeParticle( particle2 )

                particle = self.sim.createParticle( species3, newPos )
                newsingle = self.sim.insertParticle( particle )
                
                #debug
                self.sim.checkShell( newsingle )
                
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
            newCoM = oldCoM + displacement_R


            # calculate new r
            print ( rnd[3], self.a_r, self.r0, self.dt )
            theta_r = self.pgf.drawTheta( rnd[3], self.a_r, self.r0, self.dt )
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
            print 'r0 = ', self.r0, 'dt = ', self.dt, self.pgf.dump()
            r = self.pgf.drawR( rnd[0], self.r0, self.dt )
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
            newCoM = oldCoM + displacement_R


            newpos1, newpos2 = self.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
                
        else:
            raise SystemError, 'Bug: invalid eventType.'

        newpos1 = self.sim.applyBoundary( newpos1 )
        newpos2 = self.sim.applyBoundary( newpos2 )

        particle1.setPos( newpos1 )
        particle2.setPos( newpos2 )

        # check if the new positions are valid:
        # FIXME: check if positions are inside the shells, too.
        newDistance = distance( newpos1, newpos2 )
        radius12 = species1.radius + species2.radius
        # debug
        if newDistance <= radius12:
            print 'rejected move: ', 'radii, interp',\
                  species1.radius + species2.radius, newDistance
            print 'DEBUG: r0, dt, pos1, pos2, newpos1, newpos2',\
                  self.r0, self.dt, pos1, pos2, newpos1, newpos2
            raise RuntimeError, 'New particles overlap'

        # debug
        if self.sim.distance( oldCoM, newpos1 ) > self.a_r + self.a_R or \
               self.sim.distance( oldCoM, newpos2 ) > self.a_r + self.a_R:
            raise RuntimeError, 'New particle(s) out of protective sphere.'
            

        pairFormingFactor = 10

        # here decide whether this pair still continues or breaks up
        if newDistance <= radius12 * pairFormingFactor:
            dt, type = self.nextEvent()
            return dt

        else: # breaking up to singles

            single1 = self.single1
            single2 = self.single2
            
            # protect the singles with shells.
            neighbors1, distances1 = single1.getNeighborShells( n=3 )
            closest1, distance1 = neighbors1[1], drs1[1]
            if closest1 == self:  # avoid this pair. ugly
                closest1, distance1 = neighbors1[2], distances1[2]
            neighbors2, distancess2 = single2.getNeighborShells()
            closest2, distance2 = neighbors2[1], distancess2[1]
            if closest2 == self:  # avoid this pair. ugly
                closest2, distance2 = neighbors2[2], distancess2[2]
            
            print closest1, single2, closest2, single1
            if closest1 == single2 and closest2 == single1:
                single1.setShellSize( newDistance * .4999 )
                single2.setShellSize( newDistance * .4999 )
                
            else:
                print 'd1, d2', distance1, distance2
                single1.setShellSize( distance1 )
                single2.setShellSize( distance2 )


            single1.updateDt()
            single2.updateDt()
            
            self.sim.checkShell( single1 )
            self.sim.checkShell( single2 )
            
            self.sim.addEvent( self.sim.t + single1.dt, single1 )
            self.sim.addEvent( self.sim.t + single2.dt, single2 )

        return -1


    def update( self, t ):
        print 'update ', t

    def isDependentOn( self, event ):
        #print event
        return False

    def __str__( self ):
        return 'Pair of ' + str(self.single1.particle) +\
               ' and ' + str(self.single2.particle)


class EGFRDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        GFRDSimulatorBase.__init__( self )

        self.isDirty = True

        self.scheduler = EventScheduler()

        self.t = 0.0
        self.dtMax = INF
        self.dt = INF

        self.pairList = []
        self.singleMap = {}

        self.lastEvent = None
        self.hoggerCounter = 0


    def initialize( self ):

        self.scheduler.clear()

        self.initializeSingleMap()
        self.initializeSingles()

        for single in self.singleMap.values():
            nextt = single.lastTime + single.dt
            self.addEvent( nextt, single )


        #self.formPairs()

        #debug
        self.checkShellForAll()

        #self.scheduler.updateAllEventDependency()

        self.isDirty = False


    def reinitialize( self ):

        for single in self.singleMap.values():
            single.reinitialize()
            self.updateEvent( single )

        self.dt = self.scheduler.getTime() - self.t
        assert self.dt >= 0.0

        #debug
        self.checkShellForAll()


    def step( self ):

        #self.checkInvariants()

        if self.isDirty:
            self.initialize()

        self.t, self.lastEvent = self.scheduler.getTopEvent()

        print self.lastEvent
        
        self.scheduler.step()

        nextTime, nextEvent = self.scheduler.getTopEvent()
        self.dt = nextTime - self.t
        
        assert self.scheduler.getSize() != 0


        # if the same event stepped in the last n steps,
        # reinitialize everything.
        # FIXME: don't need to initialize everything.
        #        (1) recalculate shell size to the closest with
        #            its shell size shrunken.
        #        (2) new shell size is
        #            min( shell size to the second closest,
        #                 shell size to the closest, with its shell size
        #                 shrunken )

        if self.lastEvent is nextEvent:
            self.hoggerCounter += 1
        else:
            self.hoggerCounter = 0

        if self.hoggerCounter >= 10: # or self.dt < 1e-15:
                print 'reinitialize'
                self.hoggerCounter = 0
                self.reinitialize()

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

    def insertParticle( self, particle ):
        single = self.createSingle( particle )
        single.updateShell()
        single.updateDt()
        self.addEvent( self.t + single.dt, single )
        return single

    def removeParticle( self, particle ):
        single = self.findSingle( particle )
        self.removeSingle( single )
        particle.species.removeParticleBySerial( particle.serial )


    def addEvent( self, t, event ):
        event.eventID = self.scheduler.addEvent( t, event )

    def removeEvent( self, event ):
        self.scheduler.removeEvent( event.eventID )

    def updateEvent( self, event ):
        print event.eventID
        self.scheduler.updateEvent( event.eventID )


    def initializeSingleMap( self ):

        self.singleMap = {}

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
                single.setShellSize( single.getRadius() )

    def initializeSingles( self ):

        for single in self.singleMap.values():
            single.initialize()


    def formPairs( self ):

        for single in self.singleMap.values():

            closest = single.closest

            if single.particle < closest.particle:
                continue

            partnerNeighbors, partnerDistancess =\
                              self.getNeighbors( closest.particle.getPos())
            partnerClosest = Particle( partnerNeighbors[1][0],
                                       index=partnerNeighbors[1][1] )

            if single.particle == partnerClosest:
                
                pos1 = single.particle.getPos()
                pos2 = closest.particle.getPos()
                r0 = self.distance( pos1, pos2 )

                pair = self.createPair( single, closest )
                com = pair.getCoM()

                neighbors, distances = self.getNeighbors( com, 3 )
                pairDistances = distances[2]

                if pairDistance < r0:   # this happens with a small probability
                    print 'pairDistance < r0', pairDistance, r0, pairDistance - r0
                    #raise ''
                    break
                
                if pairDistance > r0 * 5:
                    pairDistance = r0 * 5

                nextEvent = pair.nextEvent( pairDistance )
                dt = nextEvent[0]
                #eventType = nextEvent[1]

                if pair.single1.eventID != None:
                    self.removeEvent( pair.single1 )
                if pair.single2.eventID != None:
                    self.removeEvent( pair.single2 )

                self.addEvent( self.t + dt, pair ) 
                self.pairList.append( pair )
                

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

            
        
    def checkShell( self, obj ):
        neighbors, distances = self.getNeighborShells( obj.getPos() )
        closest, distance = neighbors[1], distances[1]
        shellSize = obj.getShellSize()
        if distance - shellSize <= 0.0:
            print obj, closest, shellSize, distance, shellSize - distance
            raise RuntimeError, 'Fatal: shells overlap.'

    def checkShellForAll( self ):
        scheduler = self.scheduler

        size = scheduler.getSize()

        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i)[1]
            self.checkShell( obj )


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

            distances = self.distanceSqArray( pos, positions )
        
            indices = distances.argsort()[:n]
            distances = distances.take( indices )
            distances = numpy.sqrt( distances )

            distances -= species.radius

            topDistances.extend( distances )
            topNeighbors.extend( [ ( species, i ) for i in indices ] )

        topargs = numpy.argsort( topDistances )[:n]
        topDistances = numpy.take( topDistances, topargs )
        topNeighbors = [ topNeighbors[arg] for arg in topargs ]

        return topNeighbors, topDistances



    def getNeighborShells( self, pos, n=2 ):

        scheduler = self.scheduler

        size = scheduler.getSize()
        topNeighbors = [None,] * size
        topDistances = numpy.zeros( size )


        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i)[1]
            #print obj.getPos(), obj.getShellSize()
            topNeighbors[i] = obj
            topDistances[i] = self.distance( obj.getPos(), pos )\
                              - obj.getShellSize()
            
        topargs = numpy.argsort( topDistances )[:n]
        topDistances = numpy.take( topDistances, topargs )
        topNeighbors = [ topNeighbors[arg] for arg in topargs ]

        return topNeighbors, topDistances



    def getNeighborSingleShells( self, pos, n=2, speciesList=None ):

        topNeighbors = []
        topDistances = []

        if speciesList == None:
            speciesList = self.speciesList.values()

        for species in speciesList:

            # empty
            if species.pool.size == 0:
                continue

            positions = species.pool.positions
            distances = self.distanceSqArray( pos, positions )
            distances = numpy.sqrt( distances )
            distances -= species.pool.distances
            
            indices = distances.argsort()[:n]
            distances = distances.take( indices )

            topDistances.extend( distances )
            topNeighbors.extend( [ ( species, i ) for i in indices ] )

        topargs = numpy.argsort( topDistances )[:n]
        topDistances = numpy.take( topDistances, topargs )
        topNeighbors = [ topNeighbors[arg] for arg in topargs ]

        return topNeighbors, topDistances


    def checkInvariants( self ):

        assert self.t >= 0.0
        assert self.dt >= 0.0
        
        self.checkShellForAll()
        
        
