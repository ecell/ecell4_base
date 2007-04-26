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

        self.partner = None

        self.gf = FirstPassageGreensFunction( particle.species.D )

    def __del__( self ):
        print 'del', str( self )

    def setPos( self, pos ):
        self.particle.setPos( pos )

    def getPos( self ):
        return self.particle.getPos()

    def getD( self ):
        return self.particle.species.D

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


    def propagate( self, r, t ):

        rnd = numpy.random.uniform( size=2 )

        r = self.getMobilityRadius()
        displacementS = [ r, rnd[0] * Pi, rnd[1] * 2 * Pi ]
        displacement = sphericalToCartesian( displacementS )

        pos = self.particle.getPos()
        self.particle.setPos( pos + displacement )

        # BOUNDARY
        self.sim.applyBoundary( pos )

        self.lastTime = t


    '''
    Burst the protective shell.

    Shell size is shrunken to the actual radius of the particle.
    self.dt is reset to 0.0.  Do not forget to reschedule this Single
    after calling this method.
    '''

    def burstShell( self ):

        self.setShellSize( self.getRadius() )
        self.dt = 0.0
        

    '''
    Initialize this Single.

    The shell size is shrunken to the particle radius.
    self.lastTime is reset to the current time, and self.dt
    is set to zero.

    '''

    def initialize( self ):

        #neighbors, distances = self.sim.getNeighborShells( self.getPos() )
        #closestParticle = neighbors[1]
        #closestSingle = self.sim.findSingle( closestParticle )
        #self.closest = closestSingle

        self.burstShell()
        self.lastTime = self.sim.t


    def fire( self ):

        #debug
        #self.sim.checkShellForAll()

        radius = self.getRadius()
        t = self.sim.t

        # (1) propagate

        self.propagate( self.getMobilityRadius(), t )
        pos = self.getPos()

        # (2) pair check

        neighbors, distances = self.sim.getNeighborParticles( self.getPos() )
        closest = neighbors[1]
        closestDistance = distances[1]
        radius12 = radius + closest.species.radius
        closestSingle = self.sim.findSingle( closest )


        if closestSingle.partner == None:
            print t, closestSingle.lastTime, closestSingle.dt
            closestSingle.update( t ) # update the partner candidate.
            self.sim.updateEvent( closestSingle )

            pair = self.sim.createPair( self, closestSingle )
        
            if self.sim.checkPair( pair ) != None:
                # if pair was formed, destroy this single.
                self.sim.removeEvent( closestSingle )
                dt = pair.nextEvent()[0]
                self.sim.addEvent( t + dt, pair )
                return -1
            else:
                # if not, this single continues. clean up the pair.
                pair.releaseSingles()
                del pair
            
        
        # (3) neighbor check
        neighborShells, shellDistances = self.sim.getNeighborShells( pos )
        self.closest = neighborShells[1]
        distanceToClosestShell = shellDistances[1]

        # (4) determine new shell size and dt.

        shellSize = self.getShellSize()

        ShellSizeDisparityFactor = 2

        #closestMobilityRadius = self.closest.getMobilityRadius()
        #closestMeanArrivalTime = closestMobilityRadius ** 2 / \
        #( 6.0 * self.closest.getD() )
        
        shellSize = min( self.closest.getShellSize() * ShellSizeDisparityFactor
                         + ( distanceToClosestShell - radius ) * 0.5 + radius,
                         distanceToClosestShell )
        shellSize = shellSize * ( 1.0 - 1e-8 ) # safety
        shellSize = max( shellSize, radius ) # cannot be smaller than radius

        assert shellSize <= distanceToClosestShell

        self.setShellSize( shellSize )

        meanArrivalTime = self.getMobilityRadius() ** 2 / \
                          ( 6.0 * self.particle.species.D )

        self.updateDt()

        #FIXME: use of closest.dt here is a temporary workaround.
        if meanArrivalTime <= 0.0 or \
           self.closest.dt / meanArrivalTime\
           >= ShellSizeDisparityFactor * 5:
            print 'burst'
            self.closest.update( t )
            self.sim.updateEvent( self.closest )

        return self.dt


    '''
    Update the position of the particle at time t.

    t must be after the last time this Single was propagated
    (self.lastTime) but before the next scheduled time
    (self.lastTime + self.dt).

    Shell size shrunken to the radius.   self.lastTime is reset.
    self.dt is set to 0.0.

    Update the scheduler after calling this method.
    '''
    
    def update( self, t ):

        assert t >= self.lastTime
        assert t <= self.lastTime + self.dt
        assert self.getShellSize() >= self.getRadius()

        if t == self.lastTime or self.getMobilityRadius() == 0.0:
            return

        dt = t - self.lastTime

        rnd = numpy.random.uniform()
        r = self.gf.drawR( rnd , dt, self.getMobilityRadius() )
        self.propagate( r, t )

        self.burstShell()



    def updateDt( self ):
        self.dt = self.calculateFirstPassageTime()        
        

    def isDependentOn( self, event ):
        #print event
        return False

    def calculateFirstPassageTime( self ):
        
        rnd = numpy.random.uniform()
        r = self.getMobilityRadius()
        dt = self.gf.drawTime( rnd, r )
        return dt


    def __str__( self ):
        return 'Single' + str( self.particle )



class Pair:
    
    def __init__( self, sim, single1, single2, rt ):

        # Order single1 and single2 so that D1 < D2.
        if single1.particle.species.D <= single1.particle.species.D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.single1.partner = self.single2
        self.single2.partner = self.single1

        self.rt = rt
        
        self.sim = sim
        self.lastTime = 0.0
        self.dt = 0.0
        self.closest = None

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        D1, D2 = particle1.species.D, particle2.species.D
        self.D = D1 + D2
        self.sqrtD1D2 = math.sqrt( D1 / D2 )
        self.sqrtD2D1 = math.sqrt( D2 / D1 )
        
        self.sigma = particle1.species.radius + particle2.species.radius

        #self.sgf = FirstPassageGreensFunction( self.D / 4.0 )
        self.sgf = FirstPassageGreensFunction( self.D )
        self.pgf = FirstPassagePairGreensFunction( self.D, rt.k, self.sigma )

        self.eventID = None

        self.radius = max( self.single1.particle.species.radius,
                           self.single2.particle.species.radius )

        self.shellSize = self.radius



    def __del__( self ):
        print 'del', str( self )

    def getPos( self ):
        return self.getCoM()

    def getD( self ):
        return self.D

    def setShellSize( self, shellSize ):
        assert shellSize >= self.radius
        self.shellSize = shellSize

    def getShellSize( self ):
        return self.shellSize


    '''
    Calculate and return the "Center of Mass" (== CoM) of this pair.
    '''

    def getCoM( self ):

        #FIXME: what if there are boundaries?
        
        particle1 = self.single1.particle
        particle2 = self.single2.particle
        
        pos1 = particle1.getPos()
        pos2 = particle2.getPos()
        
        com = ( self.sqrtD2D1 * pos1 + self.sqrtD1D2 * pos2 ) / \
              ( self.sqrtD2D1 + self.sqrtD1D2 )
        
        return com


    def releaseSingles( self ):
        self.single1.partner = None
        self.single2.partner = None

    '''
    Calculate new positions of the pair particles using
    a new center-of-mass, a new inter-particle vector, and
    an old inter-particle vector.

    '''

    def newPositions( self, CoM, newInterParticle, oldInterParticle ):

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
        
        newpos1 = CoM - ( self.sqrtD1D2 * newInterParticle / \
                          ( self.sqrtD2D1 + self.sqrtD1D2 ) )
        newpos2 = newpos1 + newInterParticle

        return newpos1, newpos2
        

    def nextEvent( self ):

        rnd = numpy.random.uniform( size=3 )

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        species1 = particle1.species
        species2 = particle2.species
        radius1 = species1.radius
        radius2 = species2.radius
        D1 = species1.D
        D2 = species2.D

        pos1 = particle1.getPos()
        pos2 = particle2.getPos()

        self.r0 = self.sim.distance( pos1, pos2 )
        
        # determine which of particle 1 or 2 defines the r_shell
        rmin_1 = self.r0 * D1 / self.D
        rmin_2 = self.r0 * D2 / self.D
        if rmin_1 + radius1 > rmin_2 + radius2:
            rmin_factor = D1 / self.D
            rmin = rmin_1
            r_sigma = radius1
        else:
            rmin_factor = D2 / self.D
            rmin = rmin_2
            r_sigma = radius2
            
        # mobilityRadius = a_R + a_r 
        # shellSize      = mobilityRadius + r_sigma
        #            a_r > rmin + r_sigma
        
        margin = self.getShellSize() - rmin - r_sigma
        
        self.a_r = ( margin * .5 + rmin ) / rmin_factor
        self.a_R = margin * .5

        assert margin > 0.0
        assert self.a_r > self.r0
        #print 'ne', self.getShellSize(), margin, self.r0, self.a_r, self.a_R

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
        radius1 = species1.radius
        radius2 = species2.radius
        
        pos1 = particle1.getPos()
        pos2 = particle2.getPos()

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

                rnd = numpy.random.uniform( size=5 )

                # calculate new R
            
                r_R = self.sgf.drawR( rnd[0], self.dt, self.a_R )
            
                displacement_R_S = [ r_R,
                                     rnd[1] * Pi,
                                     rnd[2] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R \
                         / ( self.sqrtD2D1 + self.sqrtD1D2 )
                
                #FIXME: SURFACE
                newPos = self.sim.applyBoundary( newCoM )

                self.releaseSingles()
                self.sim.removeParticle( particle1 )
                self.sim.removeParticle( particle2 )

                particle = self.sim.createParticle( species3, newPos )
                newsingle = self.sim.createSingle( particle )
                newsingle.initialize()
                self.sim.addSingle( newsingle )

                self.sim.reactionEvents += 1
                
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
            newCoM = oldCoM + displacement_R \
                     / ( self.sqrtD2D1 + self.sqrtD1D2 )

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
            newCoM = oldCoM + displacement_R \
                     / ( self.sqrtD2D1 + self.sqrtD1D2 )

            newpos1, newpos2 = self.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
                
        else:
            raise SystemError, 'Bug: invalid eventType.'

        newpos1 = self.sim.applyBoundary( newpos1 )
        newpos2 = self.sim.applyBoundary( newpos2 )

        particle1.setPos( newpos1 )
        particle2.setPos( newpos2 )

        # debug: check if the new positions are valid:
        newDistance = distance( newpos1, newpos2 )
        radius12 = species1.radius + species2.radius
        # check 1: particles don't overlap.
        if newDistance <= radius12:
            print 'rejected move: ', 'radii, interp',\
                  species1.radius + species2.radius, newDistance
            print 'DEBUG: r0, dt, pos1, pos2, newpos1, newpos2',\
                  self.r0, self.dt, pos1, pos2, newpos1, newpos2
            raise RuntimeError, 'New particles overlap'
        # check 2: particles within mobility radius.
        if self.sim.distance( oldCoM, newpos1 ) + radius1 \
               > self.getShellSize() or \
               self.sim.distance( oldCoM, newpos2 ) + radius2 \
               > self.getShellSize():
            raise RuntimeError, 'New particle(s) out of protective sphere.'
            

        # here decide whether this pair still continues or breaks up

        if self.sim.checkPair( self ) != None:  # pair continues.

            pairClosest, pairDistance = self.findClosestShell()
            self.setShellSize( pairDistance )

            dt = self.nextEvent()[0]
            return dt

        else: # breaks up to singles

            self.releaseSingles()

            self.single1.initialize()
            self.single2.initialize()
            
            self.sim.addEvent( self.sim.t + self.single1.dt, self.single1 )
            self.sim.addEvent( self.sim.t + self.single2.dt, self.single2 )

        return -1


    def findClosestShell( self ):

        pairNeighbors, pairDistances = \
                       self.sim.getNeighborShells( self.getCoM(), n=4 )

        n = 0
        if pairNeighbors[0] in ( self, self.single1, self.single2 ):
            n = 1
            if pairNeighbors[1] in ( self, self.single1, self.single2 ):
                n = 2
                if pairNeighbors[2] in ( self, self.single1, self.single2 ):
                    n = 3
        
        pairClosest, pairDistance = pairNeighbors[n], pairDistances[n]
        assert not pairClosest in ( self.single1, self.single2 )

        return pairClosest, pairDistance


    def update( self, t ):
        print 'update ', t

    def isDependentOn( self, event ):
        #print event
        return False

    def __str__( self ):
        return 'Pair( ' + str(self.single1.particle) +\
               ', ' + str(self.single2.particle) + ' )'


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


    def initialize( self ):

        self.scheduler.clear()

        self.initializeSingleMap()

        for single in self.singleMap.values():
            single.initialize()


        for single in self.singleMap.values():
            nextt = single.lastTime + single.dt
            self.addEvent( nextt, single )

        #debug
        self.checkShellForAll()

        self.isDirty = False


    def reinitialize( self ):

        for single in self.singleMap.values():
            single.update( self.t )
            self.updateEvent( single )

        #self.dt = self.scheduler.getTime() - self.t
        assert self.dt >= 0.0

        #debug
        self.checkShellForAll()


    def step( self ):

        self.checkInvariants()

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
#         if self.lastEvent is nextEvent:
#             self.hoggerCounter += 1
#         else:
#             self.hoggerCounter = 0

#         if self.hoggerCounter >= 10: # or self.dt < 1e-15:
#             print 'reinitialize'
#             self.hoggerCounter = 0
#             self.reinitialize()
#             self.reinitializationCounter += 1

        #if self.dt == 0.0:
        #    raise 'dt=0'


        print 'dt', self.dt
        print 'reactions', self.reactionEvents,\
              'rejected moves', self.rejectedMoves
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

    def addSingle( self, single ):
        self.addEvent( self.t + single.dt, single )

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
        self.scheduler.updateEvent( event.eventID, self.t + event.dt, event )


    def initializeSingleMap( self ):

        self.singleMap = {}

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
                single.setShellSize( single.getRadius() )


    def createPair( self, single1, single2 ):

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.reactionTypeMap2.get( ( species1, species2 ) )

        return Pair( self, single1, single2, rt )


    '''
    Determine if given single forms a pair with another.

    If a pair is formed, this method returns the pair instance, with
    the partner of the given single removed from and the pair added
    to the scheduler.  It is the responsibility of the caller to remove
    the given single from the scheduler.

    If a pair is not formed, it returns None.

    '''

    def checkPair( self, pair ):


        PairMakingFactor = 5

        # pair making criteria:
        # 1. Distance between particles to form pair is closer than
        #    the total radii * PairMakingFactor, and
        # 2. Distance from the center-of-mass of the pair to the pair
        #    neighbor is larger than the distance between particles.

        single1 = pair.single1
        single2 = pair.single2
        species1 = single1.particle.species
        species2 = single2.particle.species
        D1 = species1.D
        D2 = species2.D
        radius1 = species1.radius
        radius2 = species2.radius
        radius12 = radius1 + radius2

        pairDistance = self.distance( single1.particle.getPos(),
                                      single2.particle.getPos() )

        # 1
        if pairDistance > radius12 * PairMakingFactor:
            return None
            

        pairClosest, pairClosestDistance = pair.findClosestShell()
        print pairClosest, pairClosestDistance
        
        # now, shellSize of this pair must be at minimum larger
        # than r0 * max( D1/(D1+D2)+raidus1, D2/(D1+D2)+radius2 )
        
        rmax = max(  D1 / pair.D + radius1, D2 / pair.D + radius2 ) *\
               pairDistance
        
        # 2
        if pairClosestDistance < rmax:
            print 'pairClosestDistance < rmax; %g, %g' % \
                  ( pairClosestDistance, rmax )
            return None

        pair.setShellSize( pairClosestDistance * (1.0 - 1e-8) )
        
        return pair
    
    
        
    def checkShell( self, obj ):
        neighbors, distances = self.getNeighborShells( obj.getPos() )
        closest, distance = neighbors[1], distances[1]
        shellSize = obj.getShellSize()
        if distance - shellSize < 0.0:
            raise RuntimeError,\
                  '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                  % ( str( obj ), str( closest ), shellSize, distance,\
                      shellSize - distance )

    def checkShellForAll( self ):
        scheduler = self.scheduler

        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i)[1]
            self.checkShell( obj )


    def getNeighborParticles( self, pos, n=2, speciesList=None ):

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
        topNeighbors = [ Particle( arg[0], index=arg[1] )\
                         for arg in topNeighbors ]

        return topNeighbors, topDistances


    '''
    Get neighbors simply by distance.

    This does not take into account of particle radius or shell size.
    '''
    def getNeighbors( self, pos, n=2 ):

        scheduler = self.scheduler

        size = scheduler.getSize()
        topNeighbors = [None,] * size
        topDistances = numpy.zeros( size )


        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i)[1]
            topNeighbors[i] = obj
            topDistances[i] = self.distance( obj.getPos(), pos )
            
        topargs = numpy.argsort( topDistances )[:n]
        topDistances = numpy.take( topDistances, topargs )
        topNeighbors = [ topNeighbors[arg] for arg in topargs ]

        return topNeighbors, topDistances

    '''
    Find closest shells.

    '''
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

        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i)[1]
            if isinstance( obj, Single ):
                assert obj.partner == None
