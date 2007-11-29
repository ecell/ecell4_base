#!/usr/env python

import weakref

import math

import numpy
#import scipy
#import scipy.optimize


from utils import *
from surface import *

from gfrdbase import *

class Delegate( object ):

    def __init__( self, obj, method ):
        self.obj = weakref.proxy( obj )
        self.method = method


    def __call__( self, arg ):
        return self.method( self.obj, arg )



class Single( object ):

    def __init__( self, particle, rt ):

        self.particle = particle
        self.rt = rt

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.setShellSize( self.getRadius() )
        self.eventID = None

        self.gf = FirstPassageGreensFunction( particle.species.D )


    def isPair( self ):

        return False

        
    def getD( self ):

        return self.particle.species.D

        
    def getPos( self ):

        return self.particle.pos


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

        pos = self.particle.pos
        pos += displacement

        self.particle.pos = pos


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
        if self.getD() == 0:
            self.dt = numpy.inf
            self.eventType = EventType.ESCAPE

        firstPassageTime = self.calculateEscapeTime()
        reactionTime = self.calculateReactionTime()

        if firstPassageTime <= reactionTime:
            self.dt = firstPassageTime
            self.eventType = EventType.ESCAPE
        else:
            self.dt = reactionTime
            self.eventType = EventType.REACTION


    def calculateReactionTime( self ):

        if not self.rt:
            return numpy.inf

        rnd = numpy.random.uniform()
        dt = ( 1.0 / self.rt.k ) * math.log( 1.0 / rnd )

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
    
    pos2t = cyclicTranspose( pos2, pos1, fsize )

    return ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 )


class Pair( object ):
    
    # CUTOFF_FACTOR is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__( self, single1, single2, rt, distFunc, cellSize ):

        # Order single1 and single2 so that D1 < D2.
        if single1.particle.species.D <= single2.particle.species.D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.rt = rt

        self.distance = distFunc
        self.cellSize = cellSize
        
        particle1 = self.single1.particle
        particle2 = self.single2.particle

        self.D1, self.D2 = particle1.species.D, particle2.species.D

        self.D_tot = self.D1 + self.D2
        self.D_geom = math.sqrt( self.D1 * self.D2 )  # geometric mean

        self.radius = max( particle1.species.radius,
                           particle2.species.radius )
        self.sigma = particle1.species.radius + particle2.species.radius

        self.sgf = FirstPassageGreensFunction( self.D_geom )
        self.sgf_free = FreeGreensFunction( self.D_geom )
        self.pgf = FirstPassagePairGreensFunction( self.D_tot, 
                                                   rt.k, self.sigma )
        self.pgf_free = FreePairGreensFunction( self.D_tot )
        self.pgf_nocol = FirstPassageNoCollisionPairGreensFunction( self.D_tot )

        self.eventID = None

        self.shellSize = self.radius
        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None


        self.squeezed = False


    #def __del__( self ):
    #pass
    #        print 'del', str( self )

    def initialize( self, t ):

        self.lastTime = t
        self.shellSize = self.radius
        self.dt = 0
        self.eventType = None

    def isPair( self ):

        return True

    def getPos( self ):

        return self.getCoM()

    def getD( self ):

        return self.D_tot #FIXME: is this correct?

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

        pairDistance = self.distance( self.single1.getPos(),
                                      self.single2.getPos() )
        radius = max( pairDistance * self.D1 /
                      self.D_tot + self.single1.getRadius(),
                      pairDistance * self.D2 /
                      self.D_tot + self.single2.getRadius() )
        return radius


    '''
    Calculate and return the "Center of Mass" (== CoM) of this pair.
    '''

    def getCoM( self ):

        particle1 = self.single1.particle
        particle2 = self.single2.particle
        
        pos1 = particle1.pos
        pos2 = particle2.pos

        pos2t = cyclicTranspose( pos2, pos1, self.cellSize ) #FIXME:
        
        com = ( pos1 * self.D2 + pos2t * self.D1 ) / self.D_tot
        
        return com


    def chooseSingleGreensFunction( self, t ):

        shellSize = self.a_R
        thresholdDistance = Pair.CUTOFF_FACTOR \
            * math.sqrt( 6.0 * self.D_geom * t )

        if shellSize < thresholdDistance:
            return self.sgf
        else:
            return self.sgf_free


    def choosePairGreensFunction( self, r0, t ):

        distanceFromSigma = r0 - self.sigma
        distanceFromShell = self.a_r - r0;

        thresholdDistance = Pair.CUTOFF_FACTOR * \
            math.sqrt( 6.0 * self.D_tot * t )

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
                print 'near only a'
                return self.pgf_nocol
                
            else:
                # distant from both a and sigma; 
                print 'free'
                return self.pgf_free


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

        newpos1 = CoM - rotated * ( self.D1 / self.D_tot )
        newpos2 = CoM + rotated * ( self.D2 / self.D_tot )

        return newpos1, newpos2
        

    def determineNextEvent( self ):

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        species1 = particle1.species
        species2 = particle2.species
        radius1 = species1.radius
        radius2 = species2.radius

        pos1 = particle1.pos
        pos2 = particle2.pos

        D1 = self.D1
        D2 = self.D2

        D1_factor = D1 / self.D_tot
        D2_factor = D2 / self.D_tot

        shellSize = self.getShellSize()

        sqrtD_tot = math.sqrt( self.D_tot )
        sqrtD_geom = math.sqrt( self.D_geom )

        r0 = self.distance( pos1, pos2 )

        assert r0 >= self.sigma, \
            'r0 %g, sigma %g' % ( r0, self.sigma )

        # equalize expected mean t_r and t_R.

        r0_1 = r0 * D1_factor
        r0_2 = r0 * D2_factor

        D_factor = sqrtD_tot + sqrtD_geom

        qrrtD1D25 = ( D1 * D2**5 ) ** 0.25
        qrrtD15D2 = ( D1**5 * D2 ) ** 0.25

        if qrrtD15D2 * r0 + ( qrrtD15D2 + qrrtD1D25 ) * radius1 \
                + D1 * ( sqrtD_tot * ( shellSize - radius2 ) 
                         - sqrtD_geom * radius2 )\
                - D2 * ( sqrtD_geom * r0 + sqrtD_tot * 
                         ( shellSize - radius1 ) )\
                         - qrrtD1D25 * radius2 >= 0:

            den1 = qrrtD1D25 + D1 * ( sqrtD_geom + sqrtD_tot )

            a_R_1 = sqrtD_geom * ( D2 * ( shellSize - radius1) + 
                                   D1 * ( shellSize - r0 - radius1 ) ) / den1

            a_r_1 = self.D_tot * ( sqrtD_geom * r0 + sqrtD_tot * 
                                   ( shellSize - radius1 ) ) / den1

            assert a_R_1 + a_r_1 * D1_factor + radius1 >= \
                a_R_1 + a_r_1 * D2_factor + radius2

            assert abs( a_R_1 + a_r_1 * D1_factor + radius1 - shellSize ) \
                < 1e-12 * shellSize

            self.a_r = a_r_1
            self.a_R = a_R_1
        else:
            den2 = qrrtD15D2 + D2 * ( sqrtD_geom + sqrtD_tot )

            a_R_2 = sqrtD_geom * ( D1 * ( shellSize - radius2 ) + 
                                   D2 * ( shellSize - r0 - radius2 ) ) / den2

            a_r_2 = self.D_tot * ( sqrtD_geom * r0 + sqrtD_tot * 
                                   ( shellSize - radius2 ) ) / den2

            assert a_R_2 + a_r_2 * D2_factor + radius2 >= \
                a_R_2 + a_r_2 * D1_factor + radius1

            assert abs( a_R_2 + a_r_2 * D2_factor + radius2 - shellSize ) \
                < 1e-12 * shellSize

            self.a_r = a_r_2
            self.a_R = a_R_2

        print 'r R', self.a_r, self.a_R
        print 'tr, tR', (( self.a_r - r0 ) / math.sqrt(6 * self.D_tot))**2,\
              (self.a_R / math.sqrt( 6*self.D_geom ))**2

        #print 'a_r a_R r0', self.a_r, self.a_R, r0
        assert self.a_r > 0
        assert self.a_R > 0
        #assert self.a_r > r0, '%g %g' % ( self.a_r, r0 )


        rnd = numpy.random.uniform( size=3 )

        self.sgf.seta( self.a_R )
        self.t_R = self.sgf.drawTime( rnd[0] )

        try:
            self.pgf.seta( self.a_r )
            self.t_r = self.pgf.drawTime( rnd[1], r0 )
        except:
            print 'dump', self.pgf.dump()
            raise

        #print 't_R', self.t_R, 't_r', self.t_r

        if self.t_R < self.t_r:
            self.dt = self.t_R
            self.eventType = 2
        else:
            self.dt = self.t_r
            self.eventType = self.pgf.drawEventType( rnd[2],
                                                     r0, self.t_r )

        #assert False


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
            
            oldInterParticle = particle2.pos - particle1.pos
            oldCoM = self.getCoM()

            r0 = self.distance( particle1.pos, particle2.pos )
            
            # calculate new CoM
            r_R = self.drawR_single( rnd[0], dt, self.a_R )
            
            displacement_R_S = [ r_R, rnd[1] * Pi, rnd[2] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R
            
            # calculate new interparticle
            r_r = self.drawR_pair( rnd[3], r0, dt, self.a_r )
            theta_r = self.drawTheta_pair( rnd[4], r_r, r0, dt )
            phi_r = rnd[5] * 2 * Pi
            newInterParticleS = numpy.array( [ r_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = self.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )

            dist = self.distance( newpos1, newpos2 )

            self.checkNewpos( newpos1, newpos2 )

            particle1.pos = newpos1
            particle2.pos = newpos2

        return ( self.single1, self.single2 )


    def drawR_single( self, rnd, t, shellSize ):

        gf = self.chooseSingleGreensFunction( t )
        gf.seta( shellSize )
        r = gf.drawR( rnd, t )
        while r > self.a_R: # redraw; shouldn't happen often
            print 'drawR_single: redraw'
            #self.sim.rejectedMoves += 1
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
            #self.sim.rejectedMoves += 1
            r = gf.drawR( rnd, r0, t )


        return r


    '''
    Draw theta for the pair inter-particle vector.
    '''
    def drawTheta_pair( self, rnd, r, r0, t ):

        gf = self.choosePairGreensFunction( r0, t )
        theta = gf.drawTheta( rnd, r, r0, t )

        return theta


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
            print 'DEBUG: dt, pos1, pos2, pos1, pos2',\
                  self.dt, pos1, pos2, pos1, pos2
            raise RuntimeError, 'New particles overlap'

        # check 2: particles within mobility radius.
        if self.distance( oldCoM, pos1 ) + species1.radius \
               > self.getShellSize() or \
               self.distance( oldCoM, pos2 ) + species2.radius \
               > self.getShellSize():
            raise RuntimeError, 'New particle(s) out of protective sphere.'


    def __str__( self ):
        buf = 'Pair( ' + str(self.single1.particle) +\
              ', ' + str(self.single2.particle) + ' )'
        if self.squeezed:
            buf += '; squeezed.'

        return buf



class DummySingle( object ):
    def __init__( self ):
        pass

    def getShellSize( self ):
        return 0.0

    def getRadius( self ):
        return 0.0

    def getD( self ):
        return 0.0

    def getPos( self ):
        return NOWHERE
    


class NoSpace( object ):
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

        self.stepCounter = 0

        self.smallT = 1e-8  # FIXME: is this ok?

        self.maxShellSize = INF

        self.lastEvent = None

        self.clearPopulationChanged()

        self.squeezed = 0


    def setMaxShellSize( self, maxShellSize ):

        self.maxShellSize = maxShellSize

    def getMaxShellSize( self ):

        return self.maxShellSize

    def initialize( self ):

        self.setAllRepulsive()

        self.scheduler.clear()

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
                self.addSingleEvent( single )


        #debug
        self.checkShellForAll()

        self.isDirty = False


    def stop( self, t ):

        if self.t == t:
            return

        if t >= self.scheduler.getTopEvent().getTime():
            raise RuntimeError, 'Stop time <= next event time.'

        self.t = t
        
        scheduler = self.scheduler
        
        pairList = []

        # first burst all Singles.
        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getArg()
            if obj.isPair():
                pairList.append( obj )
            else:
                self.burstSingle( obj )

        # then burst all Pairs.
        for obj in pairList:
            self.burstPair( obj )

        self.dt = 0.0


    def step( self ):

        self.clearPopulationChanged()

        if self.isDirty:
            self.initialize()

        if self.stepCounter % 100 == 0:
            self.checkInvariants()

        self.stepCounter += 1

        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getArg()

        print 't = ', self.t, ': event = ', self.lastEvent
        
        self.scheduler.step()

        nextEvent = self.scheduler.getTopEvent()
        nextTime, nextEventObject = nextEvent.getTime(), nextEvent.getArg()
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

        rt = self.getReactionType1( particle.species )
        single = Single( particle, rt )
        single.initialize( self.t )
        return single


    def createPair( self, single1, single2 ):

        assert single1.dt == 0.0
        assert single2.dt == 0.0
        assert single1.getMobilityRadius() == 0.0
        assert single2.getMobilityRadius() == 0.0

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.reactionTypeMap2.get( ( species1, species2 ) )

        pair = Pair( single1, single2, rt, self.distance, self.getCellSize() )
        pair.initialize( self.t )
        return pair


    def addEvent( self, t, func, arg ):

        return self.scheduler.addEvent( t, func, arg )


    def addSingleEvent( self, single ):

        eventID = self.addEvent( self.t + single.dt, 
                                 Delegate( self, EGFRDSimulator.fireSingle ), 
                                 single )
        single.eventID = eventID


    def addPairEvent( self, pair ):

        eventID = self.addEvent( self.t + pair.dt, 
                                 Delegate( self, EGFRDSimulator.firePair ), 
                                 pair )
        pair.eventID = eventID


    def removeEvent( self, event ):

        self.scheduler.removeEvent( event.eventID )


    def updateEvent( self, t, event ):
        self.scheduler.updateEventTime( event.eventID, t )


    def excludeVolume( self, pos, radius ):

        neighbors, distances = self.getNeighborShells( pos )
        n = numpy.searchsorted( distances, radius )
        neighbors = neighbors[:n]
        for neighbor in neighbors:
            print 'bursting', neighbor
            if isinstance( neighbor, Single ):
                self.burstSingle( neighbor )
            else:
                self.burstPair( neighbor )


    def fireSingleReaction( self, single ):

        reactantSpecies = single.particle.species
        rt = self.getReactionType1( reactantSpecies )
        oldpos = single.particle.pos.copy()
        
        if len( rt.products ) == 0:
            
            self.removeParticle( single.particle )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            single.particle.pos = NOWHERE

            if not self.checkOverlap( oldpos, productSpecies.radius ):
                print 'no space for product particle.'
                single.particle.pos = oldpos
                raise ''
                #raise NoSpace
                
            if reactantSpecies.radius < productSpecies.radius:
                self.excludeVolume( oldpos, productSpecies.radius )

            self.removeParticle( single.particle )
            newparticle = self.placeParticle( productSpecies, oldpos )
            newsingle = self.createSingle( newparticle )
            self.addSingleEvent( newsingle )
            print 'product;', newsingle

            
        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]
            
            D1 = productSpecies1.D
            D2 = productSpecies2.D
            D12 = D1 + D2
            
            single.particle.pos = NOWHERE

            #single.particle.setPos( NOWHERE )

            radius1 = productSpecies1.radius
            radius2 = productSpecies2.radius
            radius12 = radius1 + radius2

            for i in range( 100 ):
                unitVector = randomUnitVector()
                vector = unitVector * radius12 * (1.0 + 1e-10) # safety
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                newpos1 = oldpos + vector * ( D1 / D12 )
                newpos2 = oldpos - vector * ( D2 / D12 )

                #FIXME: check surfaces here
            
                newpos1 = self.applyBoundary( newpos1 )
                newpos2 = self.applyBoundary( newpos2 )

                # accept the new positions if there is enough space.
                if self.checkOverlap( newpos1, radius1 ) and \
                       self.checkOverlap( newpos2, radius2 ):
                    break
            else:
                print 'no space for product particles.'
                single.particle.pos = oldpos
                raise NoSpace

            self.excludeVolume( newpos1, productSpecies1.radius )
            self.excludeVolume( newpos2, productSpecies2.radius )

            self.removeParticle( single.particle )

            particle1 = self.placeParticle( productSpecies1, newpos1 )
            particle2 = self.placeParticle( productSpecies2, newpos2 )
            newsingle1 = self.createSingle( particle1 )
            newsingle2 = self.createSingle( particle2 )
            
            self.addSingleEvent( newsingle1 )
            self.addSingleEvent( newsingle2 )

            print 'products;', newsingle1, newsingle2

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1
        self.setPopulationChanged()



    def fireSingle( self, single ):

        #debug
        #self.checkShellForAll()

        # If this is immobile, don't move.
        D0 = single.getD()
        if D0 == 0.0:
            single.dt = numpy.inf
            single.eventType = EventType.ESCAPE
            return single.dt



        # Reaction.
        if single.eventType == EventType.REACTION:

            print 'single reaction', single
            single.gf.seta( single.getMobilityRadius() )
            r = single.gf.drawR( numpy.random.uniform(), single.dt )
            single.propagate( r, self.t )
            single.particle.pos = self.applyBoundary( single.particle.pos )


            try:
                self.fireSingleReaction( single )
            except NoSpace:
                self.rejectedMoves += 1
                self.updateEvent( self.t, single )
                single.dt = 0
                return single.dt

            single.dt = -1  # remove this Single from the Scheduler
            return single.dt

        # If not reaction, propagate.



        # (1) propagate
        #
        # Propagate this particle to the exit point on the shell.
        
        single.propagate( single.getMobilityRadius(), self.t )
        single.particle.pos = self.applyBoundary( single.particle.pos )


        self.updateEvent( self.t, single )

        # (2) Check shell size disparity.   Check if this Single needs
        #     to burst the closest Single or Pair.

        closest, distanceToClosestShell =\
                 self.getClosestShell( single.particle.pos, 
                                       ignore = [ single, ] )

        distanceToClosest = self.distance( single.particle.pos, 
                                           closest.getPos() )
        sqrtD0 = math.sqrt( D0 ) 
        radius0 = single.getRadius()

        realDistance = distanceToClosest - radius0 - closest.getRadius()
            
        SHELLSIZE_DISPARITY_FACTOR = 0.8

        criticalPoint = SHELLSIZE_DISPARITY_FACTOR * math.sqrt( D0 ) / \
                   ( math.sqrt( D0 ) + math.sqrt( closest.getD() ) ) *\
                   realDistance + radius0

        #print 'criticalPoint %g, closest shell %s, distance to shell %g' %\
        #(criticalPoint,closest,distanceToClosestShell)

        if distanceToClosestShell < criticalPoint or \
               distanceToClosestShell < radius0 * 10:

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

            #print 'partnerCandidates=',str(partnerCandidates)

            # try forming a Pair
            if len( partnerCandidates ) >= 1:

                pair = self.formPair( single, partnerCandidates[0] )

                if pair:
                    pair.determineNextEvent()

                    print pair, 'dt=', pair.dt, 'type=', pair.eventType
                    
                    self.addPairEvent( pair )
                    self.removeEvent( partnerCandidates[0] )
                    
                    for remainingCandidate in partnerCandidates[1:]:
                        self.updateSingle( remainingCandidate )
                        self.updateEvent( self.t + remainingCandidate.dt,
                                          remainingCandidate )
                        
                    single.dt = -1 # remove by rescheduling to past.
                    return single.dt

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

        return single.dt


    def updateSingle( self, single ):
        closest, distanceToClosestShell =\
                 self.getClosestShell( single.getPos(), ignore = [ single, ] )

        distanceToClosest = self.distance( single.getPos(), closest.getPos() )

        shellSize = single.calculateShellSize( closest, distanceToClosest,
                                               distanceToClosestShell )

        shellSize = min( shellSize, self.getCellSize(), self.maxShellSize )

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
        
        oldInterParticle = particle2.pos - particle1.pos

        oldCoM = self.applyBoundary( pair.getCoM() )

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
                newCoM = oldCoM + displacement_R
                
                #FIXME: SURFACE
                newPos = self.applyBoundary( newCoM )

                self.removeParticle( particle1 )
                self.removeParticle( particle2 )

                particle = self.createParticle( species3, newPos )
                newsingle = self.createSingle( particle )
                self.addSingleEvent( newsingle )

                self.reactionEvents += 1
                self.setPopulationChanged()
                
        
            else:
                raise NotImplementedError,\
                      'num products >= 2 not supported yet.'

            pair.dt = -1
            return pair.dt


        #
        # 2 Escape
        #

        r0 = self.distance( particle1.pos, particle2.pos )

        # temporary displace particles to do overlap check correctly.
        particle1.pos = NOWHERE
        particle2.pos = NOWHERE

        # 2.1 Escaping through a_r.
        if pair.eventType == EventType.ESCAPE:

            print 'escape r'

            print 'r0 = ', r0, 'dt = ', pair.dt, pair.pgf.dump()
            
            for i in range(100):

                rnd = numpy.random.uniform( size=5 )

                # calculate new R
            
                r_R = pair.drawR_single( rnd[0], pair.dt, pair.a_R )
                
                displacement_R_S = [ r_R, rnd[1] * Pi, rnd[2] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R

                # calculate new r
                theta_r = pair.drawTheta_pair( rnd[3], pair.a_r, r0,
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

                else:   # overlap check failed
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
                print 'r0 = ', r0, 'dt = ', pair.dt, pair.pgf.dump()
                r = pair.drawR_pair( rnd[0], r0, pair.dt, pair.a_r )
                print 'new r = ', r
                #assert r >= pair.sigma
            
                theta_r = pair.drawTheta_pair( rnd[1], r, r0, pair.dt )
                phi_r = rnd[2] * 2*Pi
                newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
                newInterParticle = sphericalToCartesian( newInterParticleS )
                
                # calculate new R
                displacement_R_S = [ pair.a_R, rnd[3] * Pi, rnd[4] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
            
                newCoM = oldCoM + displacement_R
                
                #print 'COM', oldCoM, newCoM

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


        if pair.squeezed:
            # make sure displaced particles don't intrude squeezer shells.
            self.excludeVolume( newpos1, radius1 )
            self.excludeVolume( newpos2, radius2 )


        pair.squeezed = False

        assert self.distance( newpos1, newpos2 ) >= pair.sigma

        particle1.pos = newpos1
        particle2.pos = newpos2

        single1, single2 = pair.single1, pair.single2

        single1.initialize( self.t )
        single2.initialize( self.t )
            
        self.addSingleEvent( single1 )
        self.addSingleEvent( single2 )

        pair.dt = -1
        return pair.dt


    def burstSingle( self, single ):
        single.burstShell( self.t )
        single.particle.pos = self.applyBoundary( single.particle.pos )
        self.updateEvent( self.t, single )


    def burstPair( self, pair ):
        single1, single2 = pair.breakUp( self.t )
        single1.particle.pos = self.applyBoundary( single1.particle.pos )
        single2.particle.pos = self.applyBoundary( single2.particle.pos )
        single1.initialize( self.t )
        single2.initialize( self.t )
        
        self.removeEvent( pair )
        self.addSingleEvent( single1 )
        self.addSingleEvent( single2 )


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

        shellSize = self.checkPairFormationCriteria( single1, single2,
                                                     pairClosest,
                                                     pairClosestShellDistance )
        if shellSize <= 0.0:  # Pair not formed
            return None

        pair = self.createPair( single1, single2 )

        # Squeezed; Pair must be formed but shell size bigger than given space.
        if shellSize > pairClosestShellDistance:
            print 'squeezed', shellSize, pairClosestShellDistance
            pair.squeezed = True
            self.squeezed += 1

        pair.setShellSize( shellSize )

        pairDistance = self.distance( single1.getPos(), single2.getPos() )
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

        # pairGap = real distance excluding radii
        pairGap = pairDistance - radius12


        #PairMakingFactor = 10
        #if pairDistance > radius12 * PairMakingFactor:
        #    return -0.0
            
        minShellSize = max( pairDistance * D1 / D12 + radius1,
                            pairDistance * D2 / D12 + radius2 )

        margin = math.sqrt( 6 * D12 * self.smallT ) # dummy
        shellSizeMargin = margin

        minShellSizeWithMargin = minShellSize + shellSizeMargin

        maxShellSize = min( self.getCellSize(), self.maxShellSize )


        # 0. Shell cannot be larger than max shell size or sim cell size.
        if minShellSizeWithMargin > maxShellSize:
            return -0.0

        # 1. Squeezed?
        if closestShellDistance <= minShellSizeWithMargin:
            print 'closest shell < minShellSize w/ margin; %g, %g' % \
                ( closestShellDistance, minShellSizeWithMargin )

            # there is enough pairGap; singles can do it better.
            if pairGap >= shellSizeMargin:
                print 'squeezed, but there is enough pairGap'
                return -0.0
            else:
                # real squeezing
                print 'squeezed; pairGap < shellSizeMargin'
                assert minShellSizeWithMargin < \
                    min( self.maxShellSize, self.getCellSize )
                return minShellSizeWithMargin

        # 2. Check if a Pair is better than two Singles.
        closestShell = closest.getShellSize()
        closestPos = closest.getPos()
        singleMobility = min( pairDistance - radius12,
                              self.distance( pos1, closestPos )
                              - closestShell - radius1,
                              self.distance( pos2, closestPos )
                              - closestShell - radius2 )
        
        pairMobility = min( closestShellDistance, maxShellSize ) - minShellSize
        if singleMobility >= pairMobility:
            print 'singleMobility %g >= pairMobility %g' %\
                  (singleMobility, pairMobility)
            return -0.0
        

        # 3. If not squeezed, and Singles are not better, then this couple
        #    must be a Pair.

        #FIXME: dummy?
        shellSize = minShellSize + ( closestShellDistance - minShellSize ) * .5
        shellSize = max( shellSize, minShellSizeWithMargin )

        shellSize = min( shellSize, self.getCellSize(), self.maxShellSize )

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



    '''
    Get neighbors simply by distance.

    This method does not take into account of particle radius or shell size.
    This method pick top n neighbors simply by the distance between given
    position pos to positions of objects (either Singles or Pairs) around.

    This method returns a tuple ( neighbors, distances ).
    '''

    def getNeighbors( self, pos, n=None ):

        scheduler = self.scheduler

        size = scheduler.getSize()
        if not n:
            n = size 

        neighbors = [DummySingle(),] * size
        positions = numpy.zeros( ( size, 3 ) )
        distances = numpy.zeros( size )

        for i in range( size ):
            obj = scheduler.getEventByIndex(i).getArg()
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

    def getNeighborShells( self, pos, n=None ):

        scheduler = self.scheduler

        size = scheduler.getSize()

        if not n:
            n = size

        neighbors = [DummySingle(),] * size
        distances = numpy.zeros( size )
        positions = numpy.zeros( ( size, 3 ) )
        shellSizes = numpy.zeros( size )

        for i in range( size ):
            obj = scheduler.getEventByIndex(i).getArg()
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
        return DummySingle(), numpy.inf


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
        return DummySingle(), numpy.inf


    #
    # consistency checkers
    #
    
    def checkShell( self, obj ):
        closest, distance = self.getClosestShell( obj.getPos(), [obj,] )
        shellSize = obj.getShellSize()

        if shellSize > self.getCellSize():
            raise RuntimeError, '%s shell size larger than simulator cell size'

        if shellSize > self.maxShellSize:
            raise RuntimeError, '%s shell size larger than maxShellSize'

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
            obj = scheduler.getEventByIndex(i).getArg()
            self.checkShell( obj )


    def checkEventStoichiometry( self ):

        population = 0
        for species in self.speciesList.values():
            population += species.pool.size
        
        eventPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
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
            print i, event.getTime(), event.getArg()

    def dump( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex(i)
            print i, event.getTime(), event.getArg(), event.getArg().getPos()



