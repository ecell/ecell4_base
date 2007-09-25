#!/usr/env python


import math

import numpy
#import scipy
#import scipy.optimize


from utils import *
from surface import *

from gfrdbase import *


class Single:

    def __init__( self, particle, sim ):

        self.particle = particle

        self.sim = sim

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.setShellSize( self.getRadius() )
        self.eventID = None

        self.gf = FirstPassageGreensFunction( particle.species.D )

    def __del__( self ):
        pass
#        print 'del', str( self )


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

        self.shellSize = shellSize


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
    Initialize this Single.

    The shell size is shrunken to the particle radius.
    self.lastTime is reset to the current time, and self.dt
    is set to zero.

    '''

    def initialize( self, t ):

        self.resetShell()
        self.lastTime = t



    def calculateShellSize( self, closest, distance, shellDistance ):

        radius1 = self.getRadius()

        D1, D2 = self.getD(), closest.getD()
        radius2 = closest.getRadius()
        radius12 = radius1 + radius2
        sqrtD1 = math.sqrt( D1 )
            
        shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                         * ( distance - radius12 ) + radius1,
                         shellDistance )

        shellSize = shellSize * ( 1.0 - 1e-10 ) # safety
        shellSize = max( shellSize, radius1 ) # not smaller than the radius

        return shellSize

        


    '''
    A mobility radius indicates the maximum displacement this single
    particle can make.

    Mobility radius of a particle is calculated as follows;

    mobility radius = shell size - radius.

    '''
    
    def getMobilityRadius( self ):
        return self.getShellSize() - self.getRadius()


    def displace( self, r ):

        rnd = numpy.random.uniform( size=2 )

        displacementS = [ r, rnd[0] * Pi, rnd[1] * 2 * Pi ]
        displacement = sphericalToCartesian( displacementS )

        pos = self.particle.getPos()
        pos += displacement

        # BOUNDARY
        pos = self.sim.applyBoundary( pos )

        self.particle.setPos( pos )

    def propagate( self, r, t ):

        self.displace( r )
        self.lastTime = t
        self.resetShell()


    '''
    Reset the protective shell.

    Shell size is shrunken to the actual radius of the particle.
    self.dt is reset to 0.0.  Do not forget to reschedule this Single
    after calling this method.
    '''

    def resetShell( self ):

        self.setShellSize( self.getRadius() )
        self.dt = 0.0
        self.eventType = EventType.ESCAPE

    def isReset( self ):
        return self.getShellSize() == self.getRadius() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE
        
        


    '''
    Update the position of the particle at time t.

    t must be after the last time this Single was propagated
    (self.lastTime) but before the next scheduled time
    (self.lastTime + self.dt).

    Shell size shrunken to the radius.   self.lastTime is reset.
    self.dt is set to 0.0.

    This method updates the scheduler.
    '''
    
    def burstShell( self, t ):

        #print 'b', self, 't ', t, 'last ', self.lastTime, 'dt ', self.dt
        assert t >= self.lastTime
        assert t <= self.lastTime + self.dt
        assert self.getShellSize() >= self.getRadius()

        dt = t - self.lastTime

        rnd = numpy.random.uniform()
        self.gf.seta( self.getMobilityRadius() )
        r = self.gf.drawR( rnd , dt )
        self.propagate( r, t )

        return self

    def determineNextEvent( self ):
        firstPassageTime = self.calculateEscapeTime()
        reactionTime = self.calculateReactionTime()

        if firstPassageTime <= reactionTime:
            self.dt = firstPassageTime
            self.eventType = EventType.ESCAPE
        else:
            self.dt = reactionTime
            self.eventType = EventType.REACTION


    def calculateReactionTime( self ):

        reactionType = self.sim.getReactionType1( self.particle.species )
        if reactionType == None:
            return numpy.inf

        rnd = numpy.random.uniform()
        dt = ( 1.0 / reactionType.k ) * math.log( 1.0 / rnd )

        return dt

    def calculateEscapeTime( self ):
        
        rnd = numpy.random.uniform()
        self.gf.seta( self.getMobilityRadius() )
        dt = self.gf.drawTime( rnd )
        return dt


    def __str__( self ):
        return 'Single' + str( self.particle )



'''
Just a free func ver of Pair.getCoM().
'''

def calculatePairCoM( pos1, pos2, D1, D2, fsize ):

    #FIXME: what if there are boundaries?
    
    sqrtD1D2 = math.sqrt( D1 / D2 )
    sqrtD2D1 = math.sqrt( D2 / D1 )

    pos2t = cyclicTranspose( pos2, pos1, fsize )

    com = ( sqrtD2D1 * pos1 + sqrtD1D2 * pos2t ) / \
          ( sqrtD2D1 + sqrtD1D2 )
    
    return com


class Pair:
    
    # CUTOFF_FACTOR is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__( self, single1, single2, rt, sim ):

        # Order single1 and single2 so that D1 < D2.
        if single1.particle.species.D <= single1.particle.species.D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.rt = rt
        
        self.sim = sim

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        self.D1, self.D2 = particle1.species.D, particle2.species.D
        self.D = self.D1 + self.D2
        self.sqrtD1D2 = math.sqrt( self.D1 / self.D2 )
        self.sqrtD2D1 = math.sqrt( self.D2 / self.D1 )
        
        self.radius = max( particle1.species.radius,
                           particle2.species.radius )
        self.sigma = particle1.species.radius + particle2.species.radius

        self.sgf = FirstPassageGreensFunction( self.D )
        self.sgf_free = FreeGreensFunction( self.D )
        self.pgf = FirstPassagePairGreensFunction( self.D, rt.k, self.sigma )
        self.pgf_free = FreePairGreensFunction( self.D )

        self.eventID = None

        self.shellSize = self.radius
        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None


        self.squeezed = False


    def __del__( self ):
        pass
    #        print 'del', str( self )

    def initialize( self, t ):
        self.lastTime = t
        self.shellSize = self.radius
        self.dt = 0
        self.eventType = None

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
        self.shellSize = shellSize

    def getShellSize( self ):
        return self.shellSize

    '''
    This method returns the radius from its CoM that this Pair must reserve
    to remain mobile.
    '''

    def getRadius( self ):  #FIXME: should be renamed?
        pairDistance = self.sim.distance( self.single1.getPos(),
                                          self.single2.getPos() )
        radius = max( pairDistance * self.D1 /
                      self.D + self.single1.getRadius(),
                      pairDistance * self.D2 /
                      self.D + self.single2.getRadius() )
        return radius

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


    def chooseSingleGreensFunction( self, t ):

        shellSize = self.a_R
        thresholdDistance = Pair.CUTOFF_FACTOR * math.sqrt( 6.0 * self.D * t )

        if shellSize < thresholdDistance:
            return self.sgf
        else:
            return self.sgf_free


    def choosePairGreensFunction( self, r0, t ):

        distanceFromSigma = r0 - self.sigma
        distanceFromShell = self.a_r - r0;

        thresholdDistance = Pair.CUTOFF_FACTOR * math.sqrt( 6.0 * self.D * t )

        if distanceFromSigma < thresholdDistance:
        
            if distanceFromShell < thresholdDistance:
                # near both a and sigma;
                # use FirstPassagePairGreensFunction
                print 'normal'
                return self.pgf
            else:
                # near sigma; use PlainPairGreensFunction

                #FIXME:
                print 'near only sigma'
                return self.pgf
        else:
            if distanceFromShell < thresholdDistance:
                # near a;

                #FIXME:
                print 'near only a'
                return self.pgf
                
            else:
                # distant from both a and sigma; 
                print 'free'
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
            rotated = numpy.array( [ newInterParticle[0], newInterParticle[1],
                                     - newInterParticle[2] ] )

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

        print 'a r0', self.a_r, self.r0
        assert self.a_r > self.r0

        rnd = numpy.random.uniform( size=3 )

        self.sgf.seta( self.a_R )
        self.t_R = self.sgf.drawTime( rnd[0] )

        try:
            self.pgf.seta( self.a_r )
            print rnd[1]
            self.t_r = self.pgf.drawTime( rnd[1], self.r0 )
        except:
            print 'dump', self.pgf.dump()
            raise

        print 't_R', self.t_R, 't_r', self.t_r

        if self.t_R < self.t_r:
            self.dt = self.t_R
            self.eventType = 2
        else:
            self.dt = self.t_r
            self.eventType = self.pgf.drawEventType( rnd[2],
                                                     self.r0, self.t_r )




    '''
    Update positions of the particles and release them as a couple of Singles.
    '''

    def breakUp( self, t ):

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
            
            displacement_R_S = [ r_R, rnd[1] * Pi, rnd[2] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R
            
            # calculate new interparticle
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

        return ( self.single1, self.single2 )


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
        buf = 'Pair( ' + str(self.single1.particle) +\
              ', ' + str(self.single2.particle) + ' )'
        if self.squeezed:
            buf += '; squeezed.'

        return buf

class SqueezingException:
    def __init( self ):
        pass
#     def __init__( self, s1, s2, s3 ):
#         self.s1 = s1
#         self.s2 = s2
#         self.s3 = s3


class EGFRDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        GFRDSimulatorBase.__init__( self )

        self.isDirty = True

        self.scheduler = EventScheduler()

        self.t = 0.0
        self.dt = INF

        self.maxDt = INF
        self.minDt = 1e-12


        self.lastEvent = None

        self.clearPopulationChanged()

        self.squeezed = 0

    def initialize( self ):

        self.setAllRepulsive()

        self.scheduler.clear()

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
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
                self.burstSingle( obj )

        # then burst all Pairs.
        for obj in pairList:
            self.burstPair( obj )

        self.dt = 0.0
#         event = self.scheduler.getTopEvent()
#         nextTime, nextEvent = event.getTime(), event.getObj()
#         self.dt = nextTime - self.t
#         assert self.dt == 0.0


    def step( self ):

        self.clearPopulationChanged()

        if self.isDirty:
            self.initialize()

        #self.checkInvariants()


        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getObj()

        print 't = ', self.t, ': event = ', self.lastEvent
        
        self.scheduler.step()

        nextEvent = self.scheduler.getTopEvent()
        nextTime, nextEventObject = nextEvent.getTime(), nextEvent.getObj()
        self.dt = nextTime - self.t

        assert self.scheduler.getSize() != 0

        print 'next dt = ', self.dt, 'reactions', self.reactionEvents,\
              'rejected moves', self.rejectedMoves, 'squeezed', self.squeezed
        assert self.scheduler.check()
        print ''


    def populationChanged( self ):
        return self.isPopulationChanged

    def clearPopulationChanged( self ):
        self.isPopulationChanged = False

    def setPopulationChanged( self ):
        self.isPopulationChanged = True
        

    def createSingle( self, particle ):
        single = Single( particle, self )
        single.initialize( self.t )
        return single

    def createPair( self, single1, single2 ):

        print single1.dt, single2.dt
        assert single1.dt == 0.0
        assert single2.dt == 0.0
        assert single1.getMobilityRadius() == 0.0
        assert single2.getMobilityRadius() == 0.0

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.reactionTypeMap2.get( ( species1, species2 ) )

        pair = Pair( single1, single2, rt, self )
        pair.initialize( self.t )
        return pair


    def addSingle( self, single ):
        self.addEvent( self.t + single.dt, single )

    def addEvent( self, t, event ):
        event.eventID = self.scheduler.addEvent( t, event )

    def removeEvent( self, event ):
        self.scheduler.removeEvent( event.eventID )

    def updateEvent( self, t, event ):
        self.scheduler.updateEvent( event.eventID, t, event )

    def fireSingleReaction( self, single ):

        single.gf.seta( single.getMobilityRadius() )
        r = single.gf.drawR( numpy.random.uniform(), single.dt )
        
        single.propagate( r, self.t )
        self.updateEvent( self.t, single )
        
        rt = self.getReactionType1( single.particle.species )
        pos = single.getPos().copy()
        
        if len( rt.products ) == 0:
            
            self.removeParticle( single.particle )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]
            
            self.removeParticle( single.particle )
            
            newparticle = self.placeParticle( productSpecies, pos )
            newsingle = self.createSingle( newparticle )
            self.addSingle( newsingle )

            print 'product;', newsingle
            
        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]
            
            D1 = productSpecies1.D
            D2 = productSpecies2.D
            
            self.removeParticle( single.particle )
            
            unitVector = randomUnitVector()
            
            #print 'unit', self.distance( unitVector, numpy.array([0,0,0]) )
            radius1 = productSpecies1.radius
            radius2 = productSpecies2.radius
            distance = radius1 + radius2
            vector = unitVector * distance * (1.0 + 1e-10) # safety
            
            # place particles according to the ratio D1:D2
            # this way, species with D=0 doesn't move.
            # FIXME: what if D1 == D2 == 0?
            newpos1 = pos + vector * ( D1 / ( D1 + D2 ) )
            newpos2 = pos - vector * ( D2 / ( D1 + D2 ) )
            
            #FIXME: check surfaces here
            
            newpos1 = self.applyBoundary( newpos1 )
            newpos2 = self.applyBoundary( newpos2 )
            
            # debug
            d = self.distance( newpos1, newpos2 )
            if d < distance:
                raise "d = %s, %s" %( d, distance)

            particle1 = self.placeParticle( productSpecies1, newpos1 )
            particle2 = self.placeParticle( productSpecies2, newpos2 )
            newsingle1 = self.createSingle( particle1 )
            newsingle2 = self.createSingle( particle2 )
            
            self.addSingle( newsingle1 )
            self.addSingle( newsingle2 )

            print 'products;', newsingle1, newsingle2

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1
        self.setPopulationChanged()



    def fireSingle( self, single ):

        #debug
        #self.checkShellForAll()

        # Reaction.
        if single.eventType == EventType.REACTION:

            print 'single reaction', single
            self.fireSingleReaction( single )

            single.dt = -1  # remove this Single from the Scheduler
            return

        # If not reaction, propagate.


        # (1) propagate
        #
        # Propagate this particle to the exit point on the shell.
        
        single.propagate( single.getMobilityRadius(), self.t )
        self.updateEvent( self.t, single )

        # (2) Check shell size disparity.   Check if this Single needs
        #     to burst the closest Single or Pair.

        closest, distanceToClosestShell =\
                 self.getClosestShell( single.getPos(), ignore = [ single, ] )

        distanceToClosest = self.distance( single.getPos(), closest.getPos() )
        D0 = single.getD()
        sqrtD0 = math.sqrt( D0 ) 
        radius0 = single.getRadius()

        realDistance = distanceToClosest - radius0 - closest.getRadius()
            
        SHELLSIZE_DISPARITY_FACTOR = 0.8

        criticalPoint = SHELLSIZE_DISPARITY_FACTOR * math.sqrt( D0 ) / \
                   ( math.sqrt( D0 ) + math.sqrt( closest.getD() ) ) *\
                   realDistance + radius0

        print 'criticalPoint %g, closest shell %s, distance to shell %g' %\
              (criticalPoint,closest,distanceToClosestShell)

        if distanceToClosestShell < criticalPoint or \
               distanceToClosestShell < radius0 * 10:

            print 'pair making'
            # (2-1) Burst the closest and do pair check with that.

            # Determine the partnerCandidate, which is the closest single.
            if closest.isPair():
                # If the partner was a Pair, then this has bursted into two
                # singles.  Find the closer one.

                self.burstPair( closest )
                candidate1 = closest.single1
                candidate2 = closest.single2
                D1 = candidate1.getD()
                D2 = candidate2.getD()
                pos0 = single.getPos()
                pos1 = candidate1.getPos()
                pos2 = candidate2.getPos()

                dist01 = self.distance( pos0, pos1 ) # radius?
                dist02 = self.distance( pos0, pos2 )
                dist12 = self.distance( pos1, pos2 )

                radius1 = candidate1.getRadius()
                radius2 = candidate2.getRadius()
                
                # MTTC = Mean Time to Correlate
                MTTC01 = meanArrivalTime( dist01 - radius0 - radius1, D0 + D1 )
                MTTC02 = meanArrivalTime( dist02 - radius0 - radius2, D0 + D2 )
                MTTC12 = meanArrivalTime( dist12 - radius1 - radius2, D1 + D2 )

                if dist01 < dist02 and MTTC01 < MTTC12:
                    partnerCandidates = [ candidate1, candidate2 ]
                elif dist02 < dist01 and MTTC02 < MTTC12:
                    partnerCandidates = [ candidate2, candidate1 ]
                else:
                    partnerCandidates = []
                
            else:  # If the closest was a Single, that is the partnerCandidate.
                self.burstSingle( closest )
                partnerCandidates = [closest,]

            print 'partnerCandidates=',str(partnerCandidates)

            # try forming a Pair
            if len( partnerCandidates ) >= 1:

                pair = self.formPair( single, partnerCandidates[0] )
                print 'pair=',pair

                if pair:
                    pair.determineNextEvent()

                    print pair, 'dt=', pair.dt, 'type=', pair.eventType
                    
                    self.addEvent( self.t + pair.dt, pair )
                    self.removeEvent( partnerCandidates[0] )
                    
                    for remainingCandidate in partnerCandidates[1:]:
                        self.updateSingle( remainingCandidate )
                        self.updateEvent( self.t + remainingCandidate.dt,
                                          remainingCandidate )
                        
                    single.dt = -1 # remove by rescheduling to past.
                    return

                else:
                    for remainingCandidate in partnerCandidates:
                        self.updateSingle( remainingCandidate )
                        self.updateEvent( self.t + remainingCandidate.dt,
                                          remainingCandidate )



        # (3) If a new Pair was not formed, this Single continues.
        #     Determine a new shell size and dt.

        # recheck the closest and distance to it.
        self.updateSingle( single )

        print 'single shell', single.getShellSize(), 'dt', single.dt


    def updateSingle( self, single ):
        closest, distanceToClosestShell =\
                 self.getClosestShell( single.getPos(), ignore = [ single, ] )

        distanceToClosest = self.distance( single.getPos(), closest.getPos() )

        shellSize = single.calculateShellSize( closest, distanceToClosest,
                                               distanceToClosestShell )

        assert shellSize <= distanceToClosestShell,\
               '%g %g' % (shellSize, distanceToClosestShell)

        single.setShellSize( shellSize )
        single.determineNextEvent()



    def firePair( self, pair ):

        print 'fire:', pair, pair.eventType

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

        #
        # 1. Reaction
        #
        if pair.eventType == EventType.REACTION:

            print 'reaction'

            if len( pair.rt.products ) == 1:
                
                species3 = pair.rt.products[0]


                
                rnd = numpy.random.uniform( size=5 )

                # calculate new R
            
                r_R = pair.drawR_single( rnd[0], pair.dt, pair.a_R )
            
                displacement_R_S = [ r_R, rnd[1] * Pi, rnd[2] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R \
                         / ( pair.sqrtD2D1 + pair.sqrtD1D2 )
                
                #FIXME: SURFACE
                newPos = self.applyBoundary( newCoM )

                self.removeParticle( particle1 )
                self.removeParticle( particle2 )

                particle = self.createParticle( species3, newPos )
                newsingle = self.createSingle( particle )
                self.addSingle( newsingle )

                self.reactionEvents += 1
                self.setPopulationChanged()
                
        
            else:
                raise NotImplementedError,\
                      'num products >= 2 not supported yet.'

            pair.dt = -1
            return


        #
        # 2 Escape
        #


        # 2.1 Escaping through a_r.
        if pair.eventType == EventType.ESCAPE:

            print 'escape r'

            print 'r0 = ', pair.r0, 'dt = ', pair.dt, pair.pgf.dump()
            
            for i in range(1000):

                rnd = numpy.random.uniform( size=5 )

                # calculate new R
            
                r_R = pair.drawR_single( rnd[0], pair.dt, pair.a_R )
                
                displacement_R_S = [ r_R, rnd[1] * Pi, rnd[2] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R \
                         / ( pair.sqrtD2D1 + pair.sqrtD1D2 )

                # calculate new r
                print ( rnd[3], pair.a_r, pair.r0, pair.dt )
                theta_r = pair.drawTheta_pair( rnd[3], pair.a_r, pair.r0,
                                               pair.dt )
                phi_r = rnd[4] * 2 * Pi
                newInterParticleS = numpy.array( [ pair.a_r, theta_r, phi_r ] )
                newInterParticle = sphericalToCartesian( newInterParticleS )
                
                newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                      oldInterParticle )
                newpos1 = self.applyBoundary( newpos1 )
                newpos2 = self.applyBoundary( newpos2 )
                
                if not pair.squeezed or \
                   ( self.checkOverlap( newpos1, radius1 ) and \
                     self.checkOverlap( newpos2, radius2 ) ):
                    break
                else:
                    self.rejectedMoves += 1
                    print '%s:ESCAPE_r: rejected move. redrawing..' % pair
            else:
                print 'redrawing limit reached.  hanging up..'
                raise RuntimeError,\
                      'redrawing limit reached under squeezing in Pair.ESCAPE_r'


        # 2.2 escaping through a_R.
        elif pair.eventType == 2:

            print 'escape R'

            for i in range(1000):
                
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
                displacement_R_S = [ pair.a_R, rnd[3] * Pi, rnd[4] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R \
                         / ( pair.sqrtD2D1 + pair.sqrtD1D2 )
                
                newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                      oldInterParticle )
                newpos1 = self.applyBoundary( newpos1 )
                newpos2 = self.applyBoundary( newpos2 )

                if not pair.squeezed or \
                   ( self.checkOverlap( newpos1, radius1 ) and \
                     self.checkOverlap( newpos2, radius2 ) ):
                    break
                else:
                    self.rejectedMoves += 1
                    print '%s:ESCAPE_r: rejected move. redrawing..' % pair
            else:
                print 'redrawing limit reached.  hanging up..'
                raise RuntimeError,\
                      'redrawing limit reached under squeezing in Pair.ESCAPE_R'

                
        else:
            raise SystemError, 'Bug: invalid eventType.'

        pair.squeezed = False

        #assert self.distance( newpos1, newpos2 ) >= pair.sigma

        particle1.setPos( newpos1 )
        particle2.setPos( newpos2 )

        print newpos1, newpos2

        single1, single2 = pair.single1, pair.single2

        single1.initialize( self.t )
        single2.initialize( self.t )
            
        self.addSingle( single1 )
        self.addSingle( single2 )

        print pair.eventID, single1.eventID, single2.eventID
        print self.scheduler.getEvent( single1.eventID ).getTime(),\
              self.scheduler.getEvent( single2.eventID ).getTime()

        pair.dt = -1
        return


    def burstSingle( self, single ):
        single.burstShell( self.t )
        self.updateEvent( self.t, single )

    def burstPair( self, pair ):
        single1, single2 = pair.breakUp( self.t )
        single1.initialize( self.t )
        single2.initialize( self.t )
        
        self.removeEvent( pair )
        self.addEvent( self.t + single1.dt, single1 )
        self.addEvent( self.t + single2.dt, single2 )


    def formPair( self, single1, single2 ):

        assert single1.isReset()
        assert single2.isReset()

        # Then, check if this pair of singles meets the pair formation
        # criteria defined in self.checkPairFormationCriteria().
        
        com = calculatePairCoM( single1.getPos(), single2.getPos(),\
                                single1.getD(), single2.getD(),\
                                self.getCellSize() )
        com = self.applyBoundary( com )
        pairClosest, pairClosestShellDistance =\
                     self.getClosestShell( com, ignore = ( single1, single2 ) )
        
        radius1 = single1.getRadius()
        radius2 = single2.getRadius()
        radius12 = radius1 + radius2
        pairDistance = self.distance( single1.getPos(), single2.getPos() )

        shellSize = self.checkPairFormationCriteria( single1, single2,
                                                     pairClosest,
                                                     pairClosestShellDistance )
        print 'pair shell size', shellSize

        if shellSize <= 0.0:  # Pair not formed
            return None

        pair = self.createPair( single1, single2 )

        # Squeezed; Pair must be formed but shell size bigger than given space.
        if shellSize > pairClosestShellDistance:
            print 'squeezed', shellSize, pairClosestShellDistance
            pair.squeezed = True
            self.squeezed += 1
            
        pair.setShellSize( shellSize * ( 1.0 - 1e-8 ) )

        print 'Pair formed: ', pair, 'pair distance', pairDistance,\
              'shell size=', pair.getShellSize(),\
              ', closest = ', pairClosest,\
              ', distance to shell = ', pairClosestShellDistance

        return pair
            

    '''
    Determine if given couple singles meet the criteria of forming
    a Pair.

    '''

    def checkPairFormationCriteria( self, single1, single2,
                                    closest, closestShellDistance ):


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
        #PairMakingFactor = 10
        #if pairDistance > radius12 * PairMakingFactor:
        #    return -0.0
            
        # 2
        
        # Shell size of this pair must be at least larger than
        # minShellSize = max( r0 * D1/(D1+D2)+raidus1, r0 * D2/(D1+D2)+radius2 )
        
        minShellSize = max( pairDistance * D1 / D12 + radius1,
                            pairDistance * D2 / D12 + radius2 )

        shellSizeMargin = 1e-9 #minShellSize * 1.01 # margin; dummy
        minShellSizeWithMargin = minShellSize + shellSizeMargin

        # pairGap = real distance including radii
        pairGap = pairDistance - radius12

        
        if closestShellDistance <= minShellSizeWithMargin:
            print 'closestShellDistance < minShellSizeWithMargin; %g, %g' % \
                  ( closestShellDistance, minShellSizeWithMargin )
            if pairGap < shellSizeMargin:
                print 'pairGap < shellSizeMargin'
                return minShellSizeWithMargin
            else:
                return -0.0

        # 3 finally, check if a Pair is better than two Singles.
        closestShell = closest.getShellSize()
        closestPos = closest.getPos()
        singleMobility = min( pairDistance - radius12,
                              self.distance( pos1, closestPos )
                              - closestShell - radius1,
                              self.distance( pos2, closestPos )
                              - closestShell - radius2 )
        
        pairMobility = closestShellDistance - minShellSize
        if singleMobility >= pairMobility:
            print 'singleMobility %g >= pairMobility %g' %\
                  (singleMobility, pairMobility)
            return -0.0
        
        #FIXME: dummy?
        shellSize = minShellSize + ( closestShellDistance - minShellSize ) * .5

        return shellSize
    
        
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


    #
    # consistency checkers
    #
    
    def checkShell( self, obj ):
        closest, distance = self.getClosestShell( obj.getPos(), [obj,] )
        shellSize = obj.getShellSize()
        if distance - shellSize < 0.0:
            if ( obj.isPair() and obj.squeezed ) or \
                   ( closest.isPair() and closest.squeezed ):
                print '%s overlaps with %s.  ignoring because squeezed.' \
                          % ( str( obj ), str( closest ) )
            else:
                raise RuntimeError,\
                      '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                      % ( str( obj ), str( closest ), shellSize, distance,\
                          distance - shellSize )

    def checkShellForAll( self ):
        scheduler = self.scheduler

        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getObj()
            self.checkShell( obj )

    def checkEventStoichiometry( self ):

        population = 0
        for species in self.speciesList.values():
            population += species.pool.size
        
        eventPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getObj()
            if obj.isPair():
                eventPopulation += 2
            else:
                eventPopulation += 1

        if population != eventPopulation:
            raise RuntimeError, 'population %d != eventPopulation %d' %\
                  ( population, eventPopulation )
        
                
        
        
    def checkInvariants( self ):

        assert self.t >= 0.0
        assert self.dt >= 0.0
        
        self.checkShellForAll()

        self.checkEventStoichiometry()

    #
    # methods for debugging.
    #


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



