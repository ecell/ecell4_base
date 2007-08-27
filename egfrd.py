#!/usr/env python


import math

import numpy
#import scipy
#import scipy.optimize


from utils import *
from surface import *

from gfrdbase import *


class Single:

    def __init__( self, sim, particle ):

        self.particle = particle
        self.sim = sim
        self.lastTime = 0.0
        self.dt = 0.0
        self.setShellSize( self.getRadius() )
        self.eventID = None

        self.partner = None

        self.gf = FirstPassageGreensFunction( particle.species.D )

    def __del__( self ):
        #pass
        print 'del', str( self )


    def isPair( self ):
        return False
        
    def getD( self ):
        return self.particle.species.D

    def fire( self ):
        print 'fireSingle', self, self.dt
        self.sim.fireSingle( self )
        print 'single new t dt', self.sim.t + self.dt, self.dt
        return self.dt
        
    def setPos( self, pos ):
        self.particle.setPos( pos )

    def getPos( self ):
        return self.particle.getPos()

    def getD( self ):
        return self.particle.species.D

    def setShellSize( self, shellSize ):

        if shellSize < self.getRadius():
            raise RuntimeError, 'shell size < radius; %g %g %g' % \
                  ( size, self.getRadius(), shellSize - self.getRadius() )

        self.shellSize = min( shellSize, self.sim.getCellSize() )


    '''
    A shell size of a particle is the distance from the current position
    of the particle to the farthest point in space that it can occupy
    when it made the maximum displacement defined by the mobility radius
    of the particle.
    '''

    def getShellSize( self ):
        return self.shellSize

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

        displacementS = [ r, rnd[0] * Pi, rnd[1] * 2 * Pi ]
        displacement = sphericalToCartesian( displacementS )

        pos = self.particle.getPos()
        self.particle.setPos( pos + displacement )

        # BOUNDARY
        pos = self.sim.applyBoundary( pos )

        self.lastTime = t


    '''
    Burst the protective shell.

    Shell size is shrunken to the actual radius of the particle.
    self.dt is reset to 0.0.  Do not forget to reschedule this Single
    after calling this method.
    '''

    def resetShell( self ):

        self.setShellSize( self.getRadius() )
        self.dt = 0.0
        

    '''
    Initialize this Single.

    The shell size is shrunken to the particle radius.
    self.lastTime is reset to the current time, and self.dt
    is set to zero.

    '''

    def initialize( self ):

        self.resetShell()
        self.lastTime = self.sim.t


    '''
    Update the position of the particle at time t.

    t must be after the last time this Single was propagated
    (self.lastTime) but before the next scheduled time
    (self.lastTime + self.dt).

    Shell size shrunken to the radius.   self.lastTime is reset.
    self.dt is set to 0.0.

    This method updates the scheduler.
    '''
    
    def burst( self, t ):

        #print 'b', self, 't ', t, 'last ', self.lastTime, 'dt ', self.dt
        assert t >= self.lastTime
        assert t <= self.lastTime + self.dt
        assert self.getShellSize() >= self.getRadius()

        dt = t - self.lastTime

        if dt != 0.0:

            rnd = numpy.random.uniform()
            self.gf.seta( self.getMobilityRadius() )
            r = self.gf.drawR( rnd , dt )
            self.propagate( r, t )  # self.lastTime = t

        self.resetShell()
        assert self.dt == 0.0
        self.sim.updateEvent( t, self )


    def calculateFirstPassageTime( self ):
        
        rnd = numpy.random.uniform()
        self.gf.seta( self.getMobilityRadius() )
        dt = self.gf.drawTime( rnd )
        return dt


    def __str__( self ):
        return 'Single' + str( self.particle )



'''
Just a free func ver of Pair.getCoM().
'''

def getPairCoM( pos1, pos2, D1, D2, fsize ):

    #FIXME: what if there are boundaries?
    
    sqrtD1D2 = math.sqrt( D1 / D2 )
    sqrtD2D1 = math.sqrt( D2 / D1 )

    pos2t = cyclicTranspose( pos2, pos1, fsize )

    com = ( sqrtD2D1 * pos1 + sqrtD1D2 * pos2t ) / \
          ( sqrtD2D1 + sqrtD1D2 )
    
    return com


class Pair:
    
    # H is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    H = 5.6


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
        self.lastTime = self.sim.t
        self.dt = 0.0

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        D1, D2 = particle1.species.D, particle2.species.D
        self.D = D1 + D2
        self.sqrtD1D2 = math.sqrt( D1 / D2 )
        self.sqrtD2D1 = math.sqrt( D2 / D1 )
        
        self.sigma = particle1.species.radius + particle2.species.radius

        #self.sgf = FirstPassageGreensFunction( self.D / 4.0 )
        self.sgf = FirstPassageGreensFunction( self.D )
        self.sgf_free = FreeGreensFunction( self.D )
        self.pgf = FirstPassagePairGreensFunction( self.D, rt.k, self.sigma )
        self.pgf_free = FreePairGreensFunction( self.D )

        self.eventID = None

        self.radius = max( self.single1.particle.species.radius,
                           self.single2.particle.species.radius )

        self.shellSize = self.radius



    def __del__( self ):
        #pass
        print 'del', str( self )

    def isPair( self ):
        return True

    def getD( self ):
        return self.D

    def fire( self ):
        self.sim.firePair( self )
        return self.dt



    def getPos( self ):
        return self.getCoM()

    def getD( self ):
        return self.D

    def setShellSize( self, shellSize ):
        #assert shellSize >= self.radius
        self.shellSize = min( shellSize, self.sim.getCellSize() )

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

        pos2t = cyclicTranspose( pos2, pos1, self.sim.getCellSize() )
        
        com = ( self.sqrtD2D1 * pos1 + self.sqrtD1D2 * pos2t ) / \
              ( self.sqrtD2D1 + self.sqrtD1D2 )
        
        return self.sim.applyBoundary( com )


    def releaseSingles( self ):
        self.single1.partner = None
        self.single2.partner = None


    def chooseSingleGreensFunction( self, t ):

        shellSize = self.a_R
        thresholdDistance = Pair.H * math.sqrt( 6.0 * self.D * t );

        if shellSize < thresholdDistance:
            return self.sgf
        else:
            return self.sgf_free


    def choosePairGreensFunction( self, r0, t ):

        distanceFromSigma = r0 - self.sigma
        distanceFromShell = self.a_r - r0;

        thresholdDistance = Pair.H * math.sqrt( 6.0 * self.D * t );

        if distanceFromSigma < thresholdDistance:
        
            if distanceFromShell < thresholdDistance:
                # near both a and sigma;
                # use FirstPassagePairGreensFunction
                return self.pgf
            else:
                # near sigma; use PlainPairGreensFunction

                #FIXME:
                return self.pgf
        else:
            if distanceFromShell < thresholdDistance:
                # near a;

                #FIXME:
                return self.pgf
                
            else:
                # distant from both a and sigma; 
                return self.pgf_free


    def drawR_single( self, rnd, t, shellSize ):

        gf = self.chooseSingleGreensFunction( t )
        gf.seta( shellSize )
        r = gf.drawR( rnd, t )
        while r > self.a_R: # redraw; shouldn't happen often
            print 'drawR_single: redraw'
            self.sim.rejectedMoves += 1
            r = gf.drawR( rnd, t )

        return r


    '''
    Draw r for the pair inter-particle vector.
    '''
    def drawR_pair( self, rnd, r0, t, a ):

        gf = self.choosePairGreensFunction( r0, t )

        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        r = gf.drawR( rnd, r0, t )
        while r > self.a_r or r <= self.sigma: # redraw; shouldn't happen often
            print 'drawR_pair: redraw'
            self.sim.rejectedMoves += 1
            r = gf.drawR( rnd, r0, t )


        return r


    '''
    Draw theta for the pair inter-particle vector.
    '''
    def drawTheta_pair( self, rnd, r, r0, t ):

        gf = self.choosePairGreensFunction( r0, t )
        print gf
        theta = gf.drawTheta( rnd, r, r0, t )

        return theta

        
    '''
    Calculate new positions of the pair particles using
    a new center-of-mass, a new inter-particle vector, and
    an old inter-particle vector.

    '''

    def newPositions( self, CoM, newInterParticle, oldInterParticle ):

        #FIXME: needs better handling of angles near zero and pi.

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotationAxis = crossproduct( [ 0,0,1 ], interParticle )

        angle = vectorAngleAgainstZAxis( oldInterParticle )
        if angle % numpy.pi != 0.0:
            rotationAxis = crossproductAgainstZAxis( oldInterParticle )
            rotationAxis = normalize( rotationAxis )
            rotated = rotateVector( newInterParticle,
                                    rotationAxis,
                                    angle )
        elif angle == 0.0:
            rotated = newInterParticle
        else:
            rotated = numpy.array( [newInterParticle[0],newInterParticle[1],
                                    -newInterParticle[2]] )

        newpos1 = CoM - ( self.sqrtD1D2 * rotated / \
                          ( self.sqrtD2D1 + self.sqrtD1D2 ) )
        newpos2 = newpos1 + rotated

        return newpos1, newpos2
        

    def determineNextEvent( self ):

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

        #FIXME: not good
        if self.r0 < self.sigma:
            self.r0 = self.sigma


        factor_1 = D1 / self.D
        factor_2 = D2 / self.D
        r0_1 = self.r0 * D1 / self.D
        r0_2 = self.r0 * D2 / self.D

        shellSize = self.getShellSize()

        margin_1 = shellSize - r0_1 - radius1
        margin_2 = shellSize - r0_2 - radius2
        
        margin = min( margin_1, margin_2 )
        assert margin > 0.0

        # FIXME: equalize expected mean t_r and t_R
        self.a_r = self.r0 + margin * .5
        self.a_R = margin * .5

        print 'ar0', self.a_r, self.r0
        assert self.a_r > self.r0

        rnd = numpy.random.uniform( size=3 )

        self.sgf.seta( self.a_R )
        self.t_R = self.sgf.drawTime( rnd[0] )

        try:
            self.pgf.seta( self.a_r )
            self.t_r = self.pgf.drawTime( rnd[1], self.r0 )
        except:
            print self.r0, self.pgf.dump()
            raise

        if self.t_R < self.t_r:
            self.dt = self.t_R
            self.eventType = 2
        else:
            self.dt = self.t_r
            self.eventType = self.pgf.drawEventType( rnd[2],
                                                     self.r0, self.t_r )




    '''
    Burst the shell, update positions of the particles and
    release them as two Singles.

    '''

    def burst( self, t ):

        assert t >= self.lastTime

        dt = t - self.lastTime 

        if dt != 0.0:

            particle1 = self.single1.particle
            particle2 = self.single2.particle
            
            rnd = numpy.random.uniform( size = 6 )
            
            oldInterParticle = particle2.getPos() - particle1.getPos()
            oldCoM = self.getCoM()
            
            # calculate new CoM
            r_R = self.drawR_single( rnd[0], dt, self.a_R )
            
            displacement_R_S = [ r_R,
                                 rnd[1] * Pi,
                                 rnd[2] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R#\
            #/ ( self.sqrtD2D1 + self.sqrtD1D2 )
            
            # calculate new interparticle
            #print ( rnd[3], self.r0, dt )
            r_r = self.drawR_pair( rnd[3], self.r0, dt, self.a_r )
            theta_r = self.drawTheta_pair( rnd[4], r_r, self.r0, dt )
            phi_r = rnd[5] * 2 * Pi
            newInterParticleS = numpy.array( [ r_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = self.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )

            newpos1 = self.sim.applyBoundary( newpos1 )
            newpos2 = self.sim.applyBoundary( newpos2 )

            dist = self.sim.distance( newpos1, newpos2 )

            self.checkNewpos( newpos1, newpos2 )

            particle1.setPos( newpos1 )
            particle2.setPos( newpos2 )



        self.releaseSingles()

        self.single1.initialize()
        self.single2.initialize()
            
        self.sim.removeEvent( self )

        self.sim.addEvent( t + self.single1.dt, self.single1 )
        self.sim.addEvent( t + self.single2.dt, self.single2 )



    def isDependentOn( self, event ):
        #print event
        return False


    def checkNewpos( self, pos1, pos2 ):

        species1 = self.single1.particle.species
        species2 = self.single2.particle.species

        oldCoM = self.getCoM()
        
        # debug: check if the new positions are valid:
        newDistance = distance( pos1, pos2 )
        radius12 = species1.radius + species2.radius

        # check 1: particles don't overlap.
        if newDistance <= radius12:
            print 'rejected move: ', 'radii, interp',\
                  species1.radius + species2.radius, newDistance
            print 'DEBUG: r0, dt, pos1, pos2, pos1, pos2',\
                  self.r0, self.dt, pos1, pos2, pos1, pos2
            raise RuntimeError, 'New particles overlap'

        # check 2: particles within mobility radius.
        if self.sim.distance( oldCoM, pos1 ) + species1.radius \
               > self.getShellSize() or \
               self.sim.distance( oldCoM, pos2 ) + species2.radius \
               > self.getShellSize():
            raise RuntimeError, 'New particle(s) out of protective sphere.'



    def __str__( self ):
        return 'Pair( ' + str(self.single1.particle) +\
               ', ' + str(self.single2.particle) + ' )'

class SqueezingException:
    def __init__( self, s1, s2, s3 ):
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3


class EGFRDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        GFRDSimulatorBase.__init__( self )

        self.isDirty = True

        self.scheduler = EventScheduler()

        self.t = 0.0
        self.dtMax = INF
        self.dt = INF

        self.lastEvent = None

        self.clearPopulationChanged()

    def initialize( self ):

        self.setAllRepulsive()

        self.scheduler.clear()

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
                single.initialize()
                self.addEvent( self.t, single )


        #debug
        self.checkShellForAll()

        self.isDirty = False

    def stop( self, t ):

        self.t = t
        
        scheduler = self.scheduler
        
        pairList = []

        # first burst all Singles.
        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getObj()
            if obj.isPair():
                pairList.append( obj )
            else:
                obj.burst( t )

        # then burst all Pairs.
        for obj in pairList:
            obj.burst( t )

        self.dt = 0.0
#         event = self.scheduler.getTopEvent()
#         nextTime, nextEvent = event.getTime(), event.getObj()
#         self.dt = nextTime - self.t
#         assert self.dt == 0.0


    def step( self ):

        #self.checkInvariants()

        self.clearPopulationChanged()

        if self.isDirty:
            self.initialize()

        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getObj()

        print 't = ', self.t, ': event = ', self.lastEvent
        
        self.scheduler.step()

        event = self.scheduler.getTopEvent()
        nextTime, nextEvent = event.getTime(), event.getObj()
        self.dt = nextTime - self.t

        assert self.scheduler.getSize() != 0

        print 'next dt = ', self.dt, 'reactions', self.reactionEvents,\
              'rejected moves', self.rejectedMoves
        assert self.scheduler.check()
        print ''


    def populationChanged( self ):
        return self.isPopulationChanged

    def clearPopulationChanged( self ):
        self.isPopulationChanged = False

    def setPopulationChanged( self ):
        self.isPopulationChanged = True
        

    def createSingle( self, particle ):
        single = Single( self, particle )
        return single

    def addSingle( self, single ):
        self.addEvent( self.t + single.dt, single )

#     def findSingleByParticle( particle ):
#         scheduler = self.scheduler
#         for i in range( scheduler.getSize() ):
#             obj = scheduler.getEventByIndex(i)[1]
#             if obj.hasattr( 'particle' ):
#                 if obj.particle == particle:
#                     return obj

#     def removeParticle( self, particle ):
#         single = self.findSingle( particle )
#         self.removeSingle( single )
#         particle.species.removeParticleBySerial( particle.serial )


    def addEvent( self, t, event ):
        assert self.scheduler.check()
        event.eventID = self.scheduler.addEvent( t, event )
        assert self.scheduler.check()

    def removeEvent( self, event ):
        assert self.scheduler.check()
        self.scheduler.removeEvent( event.eventID )
        assert self.scheduler.check()

    def updateEvent( self, t, event ):
        self.scheduler.updateEvent( event.eventID, t, event )


    def fireSingle( self, single ):

        #debug
        #self.checkShellForAll()

        radius1 = single.getRadius()

        # (1) propagate
        #
        # Propagate this particle to the exit point on the surface.
        
        single.propagate( single.getMobilityRadius(), self.t )


        # (2) pair check
        # First, find the closest particle.
        neighbors, distances = self.getNeighbors( single.getPos() )
        closest = neighbors[1]
        closestDistance = distances[1]

        print 'closest:', closest, 'dt:', closest.dt, 'distance:',    closestDistance, 'shell', closest.getShellSize()

        # Check if this and the closest particle can form a Pair.

        # Try forming a Pair if the closest is a Single.
        if not closest.isPair():

            try:
                pair = self.formPair( single, closest )
                if pair:
                    # if a Pair was formed, destroy the pair singles including
                    # self (by rescheduling to the past).
                    pair.determineNextEvent()
                    print 'pair dt', pair.dt
                    self.addEvent( self.t + pair.dt, pair )
                    self.removeEvent( closest )
                    single.dt = -1
                    return
            except SqueezingException, e:
                print 'squeezed'
                self.recoverSqueezed( e.s1, e.s2, e.s3 )
            
        
        # (3) determine new shell size and dt.

        closest, distanceToClosestShell =\
                 self.getClosestShell( single.getPos(), ignore = [ single, ] )

        pairDistance = self.distance( single.getPos(), closest.getPos() )

        ShellSizeDisparityFactor = 2

        radius1 = single.particle.species.radius

        if not closest.isPair():  # closest is a Single
            D1, D2 = single.getD(), closest.getD()
            sqrtD1, sqrtD2 = math.sqrt( D1 ), math.sqrt( D2 )
            radius2 = closest.particle.species.radius
        
            shellSize = min( sqrtD1 / ( sqrtD1 + sqrtD2 ) *
                             ( pairDistance - radius1 - radius2 ) + radius1,
                             distanceToClosestShell )
        else: 
            # the closest is a Pair.
            print 'closest is pair'
            #shellSize = max( distanceToClosestShell*.5, radius1 )
            shellSize = distanceToClosestShell


        shellSize = shellSize * ( 1.0 - 1e-8 ) # safety
        shellSize = max( shellSize, radius1 ) # cannot be smaller than radius

        assert shellSize <= distanceToClosestShell

        single.setShellSize( shellSize )
        single.dt = single.calculateFirstPassageTime()

        print 'single shell', single.getShellSize(), 'dt', single.dt

        # (4) Burst the closest, either Single or Pair, if its dt is
        #     too much larger than this.
        if closest.dt != 0.0 and single.dt * 10 < closest.dt:
            print 'burst', closest, single.dt, closest.dt
            closest.burst( self.t )

        
#         meanArrivalTime = single.getMobilityRadius() ** 2 / \
#                           ( 6.0 * single.particle.species.D )
#         if meanArrivalTime == 0.0 or \
#            closest.dt / meanArrivalTime\
#            >= ShellSizeDisparityFactor * 5:
#             print 'burst', closest, 'distance= ', distanceToClosestShell
#             closest.burst( self.t )



    def firePair( self, pair ):


        print 'fire:', pair

        particle1 = pair.single1.particle
        particle2 = pair.single2.particle
        species1 = particle1.species
        species2 = particle2.species
        radius1 = species1.radius
        radius2 = species2.radius
        D1 = species1.D
        D2 = species2.D
        
        pos1 = particle1.getPos()
        pos2 = particle2.getPos()

        oldInterParticle = pos2 - pos1

        oldCoM = pair.getCoM()


        # Three cases:
        #  1. Reaction
        #  2.1 Escaping through a_r.
        #  2.2 Escaping through a_R.

        # 1. Reaction
        if pair.eventType == EventType.REACTION:

            print 'reaction'

            if len( pair.rt.products ) == 1:
                
                species3 = pair.rt.products[0]

                rnd = numpy.random.uniform( size=5 )

                # calculate new R
            
                r_R = pair.drawR_single( rnd[0], pair.dt, pair.a_R )
            
                displacement_R_S = [ r_R,
                                     rnd[1] * Pi,
                                     rnd[2] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R \
                         / ( pair.sqrtD2D1 + pair.sqrtD1D2 )
                
                #FIXME: SURFACE
                newPos = self.applyBoundary( newCoM )

                pair.releaseSingles()
                self.removeParticle( particle1 )
                self.removeParticle( particle2 )

                particle = self.createParticle( species3, newPos )
                newsingle = self.createSingle( particle )
                newsingle.initialize()
                self.addSingle( newsingle )

                self.reactionEvents += 1
                self.setPopulationChanged()
                
                pair.dt = -1
                return

            else:
                raise NotImplementedError,\
                      'num products >= 2 not supported yet.'

        # 2.1 Escaping through a_r.
        elif pair.eventType == EventType.ESCAPE:

            print 'escape r'

            print 'r0 = ', pair.r0, 'dt = ', pair.dt, pair.pgf.dump()
            
            rnd = numpy.random.uniform( size=5 )

            # calculate new R
            
            r_R = pair.drawR_single( rnd[0], pair.dt, pair.a_R )
            
            displacement_R_S = [ r_R,
                                 rnd[1] * Pi,
                                 rnd[2] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R \
                     / ( pair.sqrtD2D1 + pair.sqrtD1D2 )

            # calculate new r
            print ( rnd[3], pair.a_r, pair.r0, pair.dt )
            theta_r = pair.drawTheta_pair( rnd[3], pair.a_r, pair.r0, pair.dt )
            phi_r = rnd[4] * 2 * Pi
            newInterParticleS = numpy.array( [ pair.a_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )


        # 2.2 escaping through a_R.
        elif pair.eventType == 2:

            print 'escape R'

            rnd = numpy.random.uniform( size = 5 )

            # calculate new r
            print 'r0 = ', pair.r0, 'dt = ', pair.dt, pair.pgf.dump()
            r = pair.drawR_pair( rnd[0], pair.r0, pair.dt, pair.a_r )
            print 'new r = ', r
            #assert r >= pair.sigma

            theta_r = pair.drawTheta_pair( rnd[1], r, pair.r0, pair.dt )
            phi_r = rnd[2] * 2*Pi
            newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            # calculate new R
            displacement_R_S = [ pair.a_R,
                                 rnd[3] * Pi,
                                 rnd[4] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R \
                     / ( pair.sqrtD2D1 + pair.sqrtD1D2 )

            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
                
        else:
            raise SystemError, 'Bug: invalid eventType.'

        newpos1 = self.applyBoundary( newpos1 )
        newpos2 = self.applyBoundary( newpos2 )

        #assert self.distance( newpos1, newpos2 ) >= pair.sigma

        particle1.setPos( newpos1 )
        particle2.setPos( newpos2 )

        print newpos1, newpos2

        # here decide whether this pair still continues or breaks up

        pairClosest, pairClosestShellDistance =\
                     self.getClosestShell( pair.getCoM(),\
                                               ( pair, pair.single1,\
                                                     pair.single2 ) )

        if self.checkPairFormationCriteria( pair.single1, pair.single2,
                                            pairClosest,
                                            pairClosestShellDistance ): 

            pairDistance = self.distance( pair.single1.getPos(),
                                          pair.single2.getPos() )

            rmax = max( pairDistance * D1 / pair.D + radius1,
                        pairDistance * D2 / pair.D + radius2 )

            shellSize = rmax + (pairClosestShellDistance-rmax)*.5

            pair.setShellSize( shellSize )

            pair.lastTime = self.t

            pair.determineNextEvent()
            return

        else: # breaks up to singles

            pair.releaseSingles()

            single1, single2 = pair.single1, pair.single2

            single1.initialize()
            single2.initialize()
            
            # singles step immediately.
            # single dts are zero.
            self.addEvent( self.t, single1 )
            self.addEvent( self.t, single2 )

            print pair.eventID, single1.eventID, single2.eventID
            print self.scheduler.getEvent( single1.eventID ).getTime(),\
                self.scheduler.getEvent( single2.eventID ).getTime()

        pair.dt = -1
        return

    def recoverSqueezed( self, single1, single2, single3 ):
        # displace single1 a bit.
        origPos = single1.getPos().copy()
        displacement = math.sqrt( single1.particle.species.D * 1e-9 * 6.0 )
        raise 'not implemented yet'

            

    def createPair( self, single1, single2 ):

        print single1.dt, single2.dt
        assert single1.dt == 0.0
        assert single2.dt == 0.0
        assert single1.getMobilityRadius() == 0.0
        assert single2.getMobilityRadius() == 0.0

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.reactionTypeMap2.get( ( species1, species2 ) )

        return Pair( self, single1, single2, rt )


    def formPair( self, single1, single2 ):

        # First, don't form a pair if either one of the singles was
        # already a member of a pair
        if single1.partner != None or single2.partner != None:
            return None

        # Then, check if this pair of singles meets the pair formation
        # criteria defined in self.checkPairFormationCriteria().
        
        com = getPairCoM( single1.getPos(), single2.getPos(),\
                          single1.particle.species.D,\
                          single2.particle.species.D, self.getCellSize() )
        com = self.applyBoundary( com )
        pairClosest, pairClosestShellDistance =\
                     self.getClosestShell( com, ignore = ( single1, single2 ) )
        

        radius1 = single1.particle.species.radius
        radius2 = single2.particle.species.radius
        pairDistance = self.distance( single1.getPos(), single2.getPos() )
        if pairDistance <= (radius1 + radius2) * (1.0+1e-6):
            neighbors, neighborDistances = \
                      self.getNeighbors( com, n=3 )
            closest, distance = neighbors[2], neighborDistances[2]
            closestRadius = closest.particle.species.radius
            if distance - closestRadius <= \
                    (radius1 + radius2) * (1.0+1e-6):
                print single1.getPos(), single2.getPos(), closest.getPos()
                print self.distance( com, closest.getPos() )
                raise SqueezingException( single1, single2, closest )
        
        if not self.checkPairFormationCriteria( single1, single2,
                                                pairClosest,
                                                pairClosestShellDistance ):
            return None

        # burst shells of both Singles.  
        single1.burst( self.t )
        single2.burst( self.t )
            
        pair = self.createPair( single1, single2 )
            
        # find closest again; singles were propagated. can be faster?
        pairClosest, pairClosestShellDistance =\
           self.getClosestShell( pair.getCoM(),\
                                 ( pair, single1, single2 ) )

        pairDistance = self.distance( single1.getPos(), single2.getPos() )

        species1 = single1.particle.species
        species2 = single2.particle.species
        D1 = species1.D
        D2 = species2.D
        D12 = D1 + D2
        radius1 = species1.radius
        radius2 = species2.radius

        rmax = max( pairDistance * D1 / D12 + radius1,
                    pairDistance * D2 / D12 + radius2 )

        shellSize = rmax + (pairClosestShellDistance-rmax)*.5

        pair.setShellSize( shellSize * ( 1.0 - 1e-8 ) )

        print 'Pair formed: ', pair, 'shell size=', pair.getShellSize(),\
              ', closest = ', pairClosest,\
              ', distance to shell = ', pairClosestShellDistance

        return pair
            

    '''
    Determine if given couple singles meet the criteria of forming
    a Pair.

    '''

    def checkPairFormationCriteria( self, single1, single2,
                                    closest, closestShellDistance ):

        # FIXME: should be larger when the GF impl is improved.
        PairMakingFactor = 5

        # pair making criteria:
        # 1. Distance between particles to form a pair is closer than
        #    the total radii * PairMakingFactor, and
        # 2. Distance from the center-of-mass of the pair to the pair
        #    neighbor is larger than the distance between particles.

        species1 = single1.particle.species
        species2 = single2.particle.species
        D1 = species1.D
        D2 = species2.D
        D12 = D1 + D2

        radius1 = species1.radius
        radius2 = species2.radius
        radius12 = radius1 + radius2

        pos1, pos2 = single1.getPos(), single2.getPos()
        pairDistance = self.distance( pos1, pos2 )

        # 1
        if pairDistance > radius12 * PairMakingFactor:
            return False
            
        # 2
        
        # Shell size of this pair must be at least larger than
        # rmax = max( r0 * D1/(D1+D2)+raidus1, r0 * D2/(D1+D2)+radius2 )
        
        # here add single12's mobility radii when the singles are still
        # on the fly.

        rmax = max( pairDistance * D1 / D12 + radius1
                    + single1.getMobilityRadius(),
                    pairDistance * D2 / D12 + radius2
                    + single2.getMobilityRadius() )
        
        if closestShellDistance < rmax * 1.0:  # margin
            print 'closestShellDistance < rmax; %g, %g' % \
                  ( closestShellDistance, rmax )
            return False

        # 3 finally, check if Pair is better than two Singles.
        closestShellSize = closest.getShellSize()
        singleMobility = min( pairDistance - radius12,
                              min( self.distance( pos1, closest.getPos() )
                                   - closestShellSize,
                                   self.distance( pos2, closest.getPos() )
                                   - closestShellSize ) )
        pairMobility = closestShellDistance - rmax
        if singleMobility >= pairMobility:
            return False
        
        return True
    
        
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
            obj = scheduler.getEventByIndex(i).getObj()
            self.checkShell( obj )


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

        if speciesList == None:
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


    '''
    Get neighbors simply by distance.

    This method does not take into account of particle radius or shell size.
    This method pick top n neighbors simply by the distance between given
    position pos to positions of objects (either Singles or Pairs) around.

    This method returns a tuple ( neighbors, distances ).
    '''

    def getNeighbors( self, pos, n=2 ):

        scheduler = self.scheduler

        size = scheduler.getSize()
        neighbors = [None,] * size
        positions = numpy.zeros( ( size, 3 ) )
        distances = numpy.zeros( size )

        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getObj()
            neighbors[i] = obj
            positions[i] = obj.getPos()

        distances = self.distanceSqArray( positions, pos )
            
        topargs = numpy.argsort( distances )[:n]
        distances = numpy.take( distances, topargs )
        distances = numpy.sqrt( distances )
        neighbors = [ neighbors[arg] for arg in topargs ]

        return neighbors, distances



    '''
    Find closest n shells.

    This method returns a tuple ( neighbors, distances ).
    '''

    def getNeighborShells( self, pos, n=2 ):

        scheduler = self.scheduler

        size = scheduler.getSize()
        neighbors = [None,] * size
        distances = numpy.zeros( size )
        positions = numpy.zeros( ( size, 3 ) )
        shellSizes = numpy.zeros( size )

        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getObj()
            neighbors[i] = obj
            positions[i] = obj.getPos()
            shellSizes[i] = obj.getShellSize()
            
        distances = self.distanceArray( positions, pos ) - shellSizes
            
        topargs = numpy.argsort( distances )[:n]
        distances = numpy.take( distances, topargs )
        neighbors = [ neighbors[arg] for arg in topargs ]

        return neighbors, distances


    '''
    Find the closest shell from the position pos.

    When the sequence parameter ignore is given, this method finds the
    closest ignoring the objects in it.

    This method returns a tuple ( neighbors, distances ).
    '''

    def getClosestShell( self, pos, ignore=[] ):

        neighbors, distances = self.getNeighborShells( pos, len( ignore ) + 1 )

        for i in range( len( neighbors ) ): 
            if neighbors[i] not in ignore:
                closest, distance = neighbors[i], distances[i]

                assert not closest in ignore
                return closest, distance

        # default case: none left.
        return None, numpy.inf

    '''
    '''

    def getClosestNeighbor( self, pos, ignore=[] ):

        neighbors, distances = self.getNeighbors( pos, len( ignore ) + 1 )

        for i in range( len( neighbors ) ): 
            if neighbors[i] not in ignore:
                closest, distance = neighbors[i], distances[i]

                assert not closest in ignore
                return closest, distance

        # default case: none left.
        return None, numpy.inf



    def dumpScheduler( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex(i)
            print i, event.getTime(), event.getObj()

    def dump( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex(i)
            print i, event.getTime(), event.getObj(), event.getObj().getPos()

        
    def checkInvariants( self ):

        assert self.t >= 0.0
        assert self.dt >= 0.0
        
        self.checkShellForAll()

        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getObj()
            if hasattr( obj, 'partner' ):
                assert obj.partner == None



# not used
#     '''
#     Find closest n shells taking into account only Singles, ignoring
#     Pairs.

#     This method returns a tuple ( neighbors, distances ).
#     '''

#     def getNeighborSingleShells( self, pos, n=2, speciesList=None ):

#         neighbors = []
#         distances = []

#         if speciesList == None:
#             speciesList = self.speciesList.values()

#         for species in speciesList:

#             # empty
#             if species.pool.size == 0:
#                 continue

#             positions = species.pool.positions
#             distances = self.distanceArray( pos, positions )
#             distances -= species.pool.distances
            
#             indices = distances.argsort()[:n]
#             distances = distances.take( indices )

#             distances.extend( distances )
#             neighbors.extend( [ ( species, i ) for i in indices ] )

#         topargs = numpy.argsort( distances )[:n]
#         distances = numpy.take( distances, topargs )
#         neighbors = [ neighbors[arg] for arg in topargs ]

#         return neighbors, distances

