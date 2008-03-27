#!/usr/env python

import weakref

import math

import numpy
#import scipy
#import scipy.optimize


from utils import *
from surface import *

from ObjectMatrix import *

from gfrdbase import *




class Delegate( object ):

    def __init__( self, obj, method ):
        self.obj = weakref.proxy( obj )
        self.method = method


    def __call__( self, arg ):
        return self.method( self.obj, arg )



class Single( object ):

    def __init__( self, particle, reactiontypes ):

        self.multiplicity = 1

        self.particle = particle
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.setRadius( self.getMinRadius() )
        self.eventID = None

        self.gf = FirstPassageGreensFunction( particle.species.D )

        self.updatek_tot()


    def getD( self ):

        return self.particle.species.D

    def getPos( self ):

        return self.particle.pos

    pos = property( getPos )


    def setRadius( self, radius ):

        if radius < self.getMinRadius():
            raise RuntimeError, 'Single radius < Particle radius; %g %g %g' % \
                  ( size, self.getMinRadius(), radius - self.getMinRadius() )

        self.radius = radius


    '''
    A radius of a Single is the distance from the current position
    of the particle to the shell of the Single.   The shell defines the
    farthest point in space that it can occupy when it makes the maximum
    displacement.
    '''

    def getRadius( self ):

        return self.radius


    def getMinRadius( self ):

        return self.particle.species.radius


    '''
    Initialize this Single.

    The radius (shell size) is shrunken to the radius of the particle
    it represents.   
    self.lastTime is reset to the current time, and self.dt
    is set to zero.
    '''

    def initialize( self, t ):

        self.reset()
        self.lastTime = t



    def calculateShellSize( self, closest, distance, shellDistance ):

        minRadius1 = self.getMinRadius()
        D1 = self.getD()

        if D1 == 0:
            return minRadius1

        D2 = closest.getD()
        minRadius2 = closest.getMinRadius()
        minRadius12 = minRadius1 + minRadius2
        sqrtD1 = math.sqrt( D1 )
            
        shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                         * ( distance - minRadius12 ) + minRadius1,
                         shellDistance )

        shellSize *= ( 1.0 - 1e-10 ) # safety
        shellSize = max( shellSize, minRadius1 ) # not smaller than the radius

        return shellSize

        


    '''
    A mobility radius indicates the maximum displacement this single
    particle can make.

    Mobility radius of a particle is calculated as follows;

    mobility radius = Single radius - Particle radius.

    '''
    
    def getMobilityRadius( self ):

        return self.radius - self.getMinRadius()


    def drawDisplacement( self, r ):

        rnd = numpy.random.uniform( size=2 )

        displacementS = [ r, rnd[0] * Pi, rnd[1] * 2 * Pi ]
        displacement = sphericalToCartesian( displacementS )

        return displacement


    '''
    Reset the Single.

    Radius (shell size) is shrunken to the actual radius of the particle.
    self.dt is reset to 0.0.  Do not forget to reschedule this Single
    after calling this method.
    '''

    def reset( self ):

        self.setRadius( self.getMinRadius() )
        self.dt = 0.0
        self.eventType = EventType.ESCAPE


    def isReset( self ):

        return self.radius == self.getMinRadius() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE
        

    def drawR( self, dt ):
        assert dt >= 0

        self.gf.seta( self.getMobilityRadius() )
        rnd = numpy.random.uniform()
        r = self.gf.drawR( rnd , dt )
        return r



    def determineNextEvent( self ):
        if self.getD() == 0:
            firstPassageTime = numpy.inf
        else:
            firstPassageTime = self.drawEscapeTime()
            
        reactionTime = self.drawReactionTime()

        if firstPassageTime <= reactionTime:
            self.dt = firstPassageTime
            self.eventType = EventType.ESCAPE
        else:
            self.dt = reactionTime
            self.eventType = EventType.REACTION


    def drawReactionTime( self ):
        
        if self.k_tot == 0:
            return numpy.inf

        rnd = numpy.random.uniform()
        dt = ( 1.0 / self.k_tot ) * math.log( 1.0 / rnd )

        return dt


    def drawEscapeTime( self ):
        
        rnd = numpy.random.uniform()
        self.gf.seta( self.getMobilityRadius() )
        dt = self.gf.drawTime( rnd )
        return dt


    def updatek_tot( self ):

        self.k_tot = 0

        if not self.reactiontypes:
            return

        for rt in self.reactiontypes:
            self.k_tot += rt.k


    def drawReactionType( self ):

        k_array = [ rt.k for rt in self.reactiontypes ]
        k_array = numpy.add.accumulate( k_array )
        k_max = k_array[-1]

        rnd = numpy.random.uniform()
        i = numpy.searchsorted( k_array, rnd * k_max )

        return self.reactiontypes[i]


    def __str__( self ):
        return 'Single' + str( self.particle )



'''
Just a free func ver of Pair.getCoM().
'''

def calculatePairCoM( pos1, pos2, D1, D2, worldSize ):

    pos2t = cyclicTranspose( pos2, pos1, worldSize )

    return ( ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 ) ) % worldSize


class Pair( object ):
    
    # CUTOFF_FACTOR is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__( self, single1, single2, rt, distFunc, worldSize ):

        self.multiplicity = 2

        # Order single1 and single2 so that D1 < D2.
        if single1.particle.species.D <= single2.particle.species.D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.rt = rt

        self.distance = distFunc
        self.worldSize = worldSize
        
        particle1 = self.single1.particle
        particle2 = self.single2.particle

        self.D1, self.D2 = particle1.species.D, particle2.species.D

        self.D_tot = self.D1 + self.D2
        self.D_geom = math.sqrt( self.D1 * self.D2 )  # geometric mean

        #self.minRadius = max( particle1.species.radius,
        #                      particle2.species.radius )
        self.sigma = particle1.species.radius + particle2.species.radius

        self.sgf = FirstPassageGreensFunction( self.D_geom )
        self.pgf = FirstPassagePairGreensFunction( self.D_tot, 
                                                   rt.k, self.sigma )
        self.pgf_free = FreePairGreensFunction( self.D_tot )
        self.pgf_basic = BasicPairGreensFunction( self.D_tot, rt.k, self.sigma )
        self.pgf_nocol = FirstPassageNoCollisionPairGreensFunction( self.D_tot )

        self.eventID = None

        self.radius = self.getMinRadius()
        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None


        self.squeezed = False


    #def __del__( self ):
    #pass
    #        print 'del', str( self )

    def initialize( self, t ):

        self.lastTime = t
        self.radius = self.getMinRadius()
        self.dt = 0
        self.eventType = None

    def getPos( self ):

        return self.getCoM()

    pos = property( getPos )

    def getD( self ):

        return self.D_tot #FIXME: is this correct?

    def setRadius( self, radius ):

        assert radius >= self.getMinRadius()
        self.radius = radius


    def getRadius( self ):

        return self.radius


    '''
    This method returns the radius from its CoM that this Pair must reserve
    to remain mobile.
    '''

    def getMinRadius( self ):

        pairDistance = self.distance( self.single1.pos,
                                      self.single2.pos )
        minRadius = max( pairDistance * self.D1 /
                         self.D_tot + self.single1.getMinRadius(),
                         pairDistance * self.D2 /
                         self.D_tot + self.single2.getMinRadius() )
        return minRadius


    '''
    Calculate and return the "Center of Mass" (== CoM) of this pair.
    '''

    def getCoM( self ):

        particle1 = self.single1.particle
        particle2 = self.single2.particle
        
        pos1 = particle1.pos
        pos2 = particle2.pos

        pos2t = cyclicTranspose( pos2, pos1, self.worldSize ) #FIXME:
        
        com = ( pos1 * self.D2 + pos2t * self.D1 ) / self.D_tot

        return com % self.worldSize


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
                # near sigma; use BasicPairGreensFunction
                print 'near only sigma'
                return self.pgf_basic
                #return self.pgf
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
            #rotated = newInterParticle * numpy.array( [ 1, 1, -1 ] )

        newpos1 = CoM - rotated * ( self.D1 / self.D_tot )
        newpos2 = CoM + rotated * ( self.D2 / self.D_tot )

        return newpos1, newpos2
        

    def determineNextEvent( self ):

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        species1 = particle1.species
        species2 = particle2.species
        particleRadius1 = species1.radius
        particleRadius2 = species2.radius

        pos1 = particle1.pos
        pos2 = particle2.pos

        D1 = self.D1
        D2 = self.D2

        D1_factor = D1 / self.D_tot
        D2_factor = D2 / self.D_tot

        shellSize = self.radius

        sqrtD_tot = math.sqrt( self.D_tot )
        sqrtD_geom = math.sqrt( self.D_geom )

        r0 = self.distance( pos1, pos2 )

        assert r0 >= self.sigma, \
            '%s;  r0 %g < sigma %g' % ( self, r0, self.sigma )

        # equalize expected mean t_r and t_R.

        r0_1 = r0 * D1_factor
        r0_2 = r0 * D2_factor

        D_factor = sqrtD_tot + sqrtD_geom

        qrrtD1D25 = ( D1    * D2**5 ) ** 0.25
        qrrtD15D2 = ( D1**5 * D2 ) ** 0.25

        if qrrtD15D2 * r0 + ( qrrtD15D2 + qrrtD1D25 ) * particleRadius1 \
                + D1 * ( sqrtD_tot * ( shellSize - particleRadius2 ) 
                         - sqrtD_geom * particleRadius2 )\
                - D2 * ( sqrtD_geom * r0 + sqrtD_tot * 
                         ( shellSize - particleRadius1 ) )\
                         - qrrtD1D25 * particleRadius2 >= 0:

            den1 = qrrtD1D25 + D1 * ( sqrtD_geom + sqrtD_tot )

            a_R_1 = sqrtD_geom * ( D2 * ( shellSize - particleRadius1) + 
                                   D1 * ( shellSize - r0 - particleRadius1 ) ) / den1

            a_r_1 = self.D_tot * ( sqrtD_geom * r0 + sqrtD_tot * 
                                   ( shellSize - particleRadius1 ) ) / den1

            assert a_R_1 + a_r_1 * D1_factor + particleRadius1 >= \
                a_R_1 + a_r_1 * D2_factor + particleRadius2

            assert abs( a_R_1 + a_r_1 * D1_factor + particleRadius1 - shellSize ) \
                < 1e-12 * shellSize

            self.a_r = a_r_1
            self.a_R = a_R_1
        else:
            den2 = qrrtD15D2 + D2 * ( sqrtD_geom + sqrtD_tot )

            a_R_2 = sqrtD_geom * ( D1 * ( shellSize - particleRadius2 ) + 
                                   D2 * ( shellSize - r0 - particleRadius2 ) ) / den2

            a_r_2 = self.D_tot * ( sqrtD_geom * r0 + sqrtD_tot * 
                                   ( shellSize - particleRadius2 ) ) / den2

            assert a_R_2 + a_r_2 * D2_factor + particleRadius2 >= \
                a_R_2 + a_r_2 * D1_factor + particleRadius1

            assert abs( a_R_2 + a_r_2 * D2_factor + particleRadius2 - shellSize ) \
                < 1e-12 * shellSize

            self.a_r = a_r_2
            self.a_R = a_R_2

        #print 'r R', self.a_r, self.a_R
        #print 'tr, tR', (( self.a_r - r0 ) / math.sqrt(6 * self.D_tot))**2,\
        #      (self.a_R / math.sqrt( 6*self.D_geom ))**2

        #print 'a a_r a_R', shellSize, self.a_r, self.a_R
        assert self.a_r > 0
        assert self.a_R > 0 or ( self.a_R == 0 and ( D1 == 0 or D2 == 0 ) )
        assert self.a_r > r0, '%g %g' % ( self.a_r, r0 )

        rnd = numpy.random.uniform( size=3 )

        # draw t_R
        self.sgf.seta( self.a_R )
        self.t_R = self.sgf.drawTime( rnd[0] )

        # draw t_r
        self.pgf.seta( self.a_r )
        self.t_r = self.pgf.drawTime( rnd[1], r0 )

        # draw t_reaction
        t_reaction1 = self.single1.drawReactionTime()
        t_reaction2 = self.single2.drawReactionTime()

        if t_reaction1 < t_reaction2:
            self.t_single_reaction = t_reaction1
            self.reactingsingle = self.single1
        else:
            self.t_single_reaction = t_reaction2
            self.reactingsingle = self.single2

        #print ( self.t_R, self.t_r, self.t_single_reaction )
        self.dt = min( self.t_R, self.t_r, self.t_single_reaction )

        assert self.dt >= 0
        print 'dt ', self.dt, 't_R', self.t_R, 't_r', self.t_r
        print self.pgf.dump()
        if self.dt == self.t_r:  # type = 0 (REACTION) or 1 (ESCAPE_r)
            self.eventType = self.pgf.drawEventType( rnd[2],
                                                     r0, self.t_r )
        elif self.dt == self.t_R: # type = ESCAPE_R (2)
            self.eventType = 2
        elif self.dt == self.t_single_reaction:  # type = single reaction (3)
            self.eventType = 3 
        else:
            raise 'never get here'

        #assert False




    def drawR_single( self, t, a ):

        self.sgf.seta( a )
        r = self.sgf.drawR( numpy.random.uniform(), t )
        while r > self.a_R: # redraw; shouldn't happen often
            print 'drawR_single: redraw'
            r = self.sgf.drawR( numpy.random.uniform(), t )

        return r


    '''
    Draw r for the pair inter-particle vector.
    '''
    def drawR_pair( self, r0, t, a ):

        gf = self.choosePairGreensFunction( r0, t )

        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        r = gf.drawR( numpy.random.uniform(), r0, t )

        while r >= self.a_r or r <= self.sigma: # redraw; shouldn't happen often
            print 'drawR_pair: redraw'
            self.sim.rejectedMoves += 1
            r = gf.drawR( numpy.random.uniform(), r0, t )

        return r


    '''
    Draw theta for the pair inter-particle vector.
    '''
    def drawTheta_pair( self, rnd, r, r0, t ):

        gf = self.choosePairGreensFunction( r0, t )
        theta = gf.drawTheta( rnd, r, r0, t )

        return theta


    def checkNewpos( self, pos1, pos2, com ):

        species1 = self.single1.particle.species
        species2 = self.single2.particle.species

        oldCoM = com
        
        # debug: check if the new positions are valid:
        newDistance = distance( pos1, pos2 )
        particleRadius12 = species1.radius + species2.radius

        # check 1: particles don't overlap.
        if newDistance <= particleRadius12:
            print 'rejected move: ', 'radii, interp',\
                  species1.radius + species2.radius, newDistance
            print 'DEBUG: dt, pos1, pos2, pos1, pos2',\
                  self.dt, pos1, pos2, pos1, pos2
            raise RuntimeError, 'New particles overlap'

        # check 2: particles within mobility radius.
        if self.distance( oldCoM, pos1 ) + species1.radius \
               > self.radius or \
               self.distance( oldCoM, pos2 ) + species2.radius \
               > self.radius:
            raise RuntimeError, 'New particle(s) out of protective sphere.'


    def __str__( self ):
        buf = 'Pair( ' + str(self.single1.particle) +\
              ', ' + str(self.single2.particle) + ' )'
        if self.squeezed:
            buf += '; squeezed.'

        return buf


class Multi( object ):
    def __init__( self ):
        self.multiplicity = 0
        self.singleList = []

    def getMinRadius( self ):
        assert False  #FIXME:
        #return 0.0

    def getD( self ):
        return 0.0

    def getPos( self ):
        return NOWHERE

    pos = property( getPos )

    def add( self, single ):
        assert isinstance( single, Single )
        self.singleList.append( single )

    def remove( self, single ):
        self.singleList.remove( single )

    def run( self, T ):
        pass

    def check( self ):

        for s in self.singleList:
            assert not s.isReset()
            #s.radius
        

    def __str__( self ):
        buf = 'Multi( '
        for s in self.singleList:
            buf += str( s.particle )
        return buf


class DummySingle( object ):
    def __init__( self ):
        self.multiplicity = 1

        self.radius = 0.0

    def getMinRadius( self ):
        return 0.0

    def getD( self ):
        return 0.0

    def getPos( self ):
        return NOWHERE

    pos = property( getPos )

    


class EGFRDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        #self.shellMatrix = ObjectMatrix()
        self.shellMatrix = SimpleObjectMatrix()

        GFRDSimulatorBase.__init__( self )

        self.isDirty = True
        self.scheduler = EventScheduler()

        self.smallT = 1e-8  # FIXME: is this ok?

        self.maxShellSize = INF

        self.reset()

    def setWorldSize( self, size ):
        GFRDSimulatorBase.setWorldSize( self, size )
        self.shellMatrix.setWorldSize( size )

    def setMatrixSize( self, size ):
        GFRDSimulatorBase.setMatrixSize( self, size )
        self.shellMatrix.setMatrixSize( size )

    def getMatrixCellSize( self ):

        return self.shellMatrix.cellSize

    def getNextTime( self ):
        return self.scheduler.getNextTime()

    def setMaxShellSize( self, maxShellSize ):

        self.maxShellSize = maxShellSize

    def getMaxShellSize( self ):

        return self.maxShellSize

    def reset( self ):

        self.t = 0.0
        self.dt = INF
        self.stepCounter = 0
        self.zeroSteps = 0
        self.rejectedMoves = 0
        self.reactionEvents = 0
        self.lastEvent = None
        self.clearPopulationChanged()
        self.squeezed = 0

        self.isDirty = True
        #self.initialize()
        

    def initialize( self ):

        GFRDSimulatorBase.initialize( self )

        self.setAllRepulsive()

        self.scheduler.clear()
        self.shellMatrix.clear()

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
                self.addSingleEvent( single )

        self.isDirty = False


    def stop( self, t ):

        print 'stop at', t

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
            if isinstance( obj, Pair ):
                pairList.append( obj )
            elif isinstance( obj, Single ):
                try:
                    self.burstSingle( obj )
                except NoSpace:
                    self.rejectedMoves += 1
            else:
                raise NotImplementedError


        # then burst all Pairs.
        for obj in pairList:
            single1, single2 = self.burstPair( obj )
            self.removeEvent( obj )
            self.addSingleEvent( single1 )
            self.addSingleEvent( single2 )

        self.dt = 0.0


    def step( self ):

        self.clearPopulationChanged()

        if self.isDirty:
            self.initialize()

        if self.stepCounter % 100 == 0:
            self.check()

        self.stepCounter += 1

        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getArg()

        print self.stepCounter, ': t = ', self.t, ': event = ', self.lastEvent
        
        self.scheduler.step()

        nextEvent = self.scheduler.getTopEvent()
        nextTime, nextEventObject = nextEvent.getTime(), nextEvent.getArg()
        self.dt = nextTime - self.t


        # assert if not too many successive dt=0 steps occur.
        if self.dt == 0:
            self.zeroSteps += 1
        else:
            self.zeroSteps = 0
        assert self.zeroSteps < self.scheduler.getSize() * 2,\
            'too many dt=zero steps.  simulator halted?'


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

        self.shellMatrix.add( single, single.pos, single.radius )

        return single


    def createPair( self, single1, single2 ):

        assert single1.dt == 0.0
        assert single2.dt == 0.0
        assert single1.getMobilityRadius() == 0.0
        assert single2.getMobilityRadius() == 0.0

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.reactionTypeMap2.get( ( species1, species2 ) )

        pair = Pair( single1, single2, rt, self.distance, self.getWorldSize() )
        pair.initialize( self.t )

        self.shellMatrix.remove( single1 )
        self.shellMatrix.remove( single2 )

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
                try:
                    self.burstSingle( neighbor )
                except NoSpace:
                    self.rejectedMoves += 1
            else:  # Pair
                single1, single2 = self.burstPair( neighbor )
                self.removeEvent( neighbor )
                self.addSingleEvent( single1 )
                self.addSingleEvent( single2 )

    def excludeSingleVolume( self, pos, radius ): # to be removed

        neighbors, distances = self.getNeighborShells( pos )
        n = numpy.searchsorted( distances, radius )
        neighbors = neighbors[:n]
        for neighbor in neighbors:
            print 'bursting', neighbor
            if isinstance( neighbor, Single ):
                try:
                    self.burstSingle( neighbor )
                except NoSpace:
                    self.rejectedMoves += 1
            else:  # Pair
                neighbor.squeezed = True


    def fireSingleReaction( self, single ):

        reactantSpecies = single.particle.species
        oldpos = single.particle.pos.copy()
        
        rt = single.drawReactionType()

        if len( rt.products ) == 0:
            
            self.removeParticle( single.particle )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            if not self.checkOverlap( oldpos, productSpecies.radius,
                                      ignore = [ single.particle, ] ):
                print 'no space for product particle.'
                raise NoSpace()
                
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
            
            particleRadius1 = productSpecies1.radius
            particleRadius2 = productSpecies2.radius
            particleRadius12 = particleRadius1 + particleRadius2

            # clean up space.
            rad = max( particleRadius12 * ( D1 / D12 ) + particleRadius1,
                       particleRadius12 * ( D2 / D12 ) + particleRadius2 )
            self.excludeVolume( oldpos, rad )


            for i in range( 100 ):
                unitVector = randomUnitVector()
                vector = unitVector * particleRadius12 * (1.0 + 1e-10) # safety
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                newpos1 = oldpos + vector * ( D1 / D12 )
                newpos2 = oldpos - vector * ( D2 / D12 )

                #FIXME: check surfaces here
            
                newpos1 = self.applyBoundary( newpos1 )
                newpos2 = self.applyBoundary( newpos2 )

                # accept the new positions if there is enough space.
                if self.checkOverlap( newpos1, particleRadius1,
                                      ignore = [ single.particle, ] ) \
                and \
                    self.checkOverlap( newpos2, particleRadius2,
                                       ignore = [ single.particle, ] ):
                    break
            else:
                print 'no space for product particles.'
                raise NoSpace()

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


    def propagateSingle( self, single, r ):
        
        closest, distanceToClosestShell =\
                 self.getClosestShell( single.particle.pos, 
                                       ignore = [ single, ] )
        squeezed = False;


        # check if this Single is squeezed.
        if closest.multiplicity != 1 and distanceToClosestShell < single.radius:
            print 'single ', single, ' squeezed by ', closest, 'distance ', \
                distanceToClosestShell
            #assert closest.squeezed,\
            #'When Single is squeezed, the closest must be a squeezed Pair'
            squeezed = True


        particleRadius = single.getMinRadius()
        oldpos = single.particle.pos.copy()
        
        for i in range( 100 ):
            
            displacement = single.drawDisplacement( r )
            
            newpos = oldpos + displacement
            
            if self.checkOverlap( newpos, particleRadius,
                                  ignore = [ single.particle, ] ):
                break
            
        else:
            single.initialize( self.t )
            self.shellMatrix.update( single, single.pos, single.radius )
            raise NoSpace()
        
        newpos = self.applyBoundary( newpos )
        self.moveParticle( single.particle, newpos )

        single.initialize( self.t )

        self.shellMatrix.update( single, single.pos, single.radius )

        return squeezed



    def fireSingle( self, single ):

        print single, single.dt

        # Reaction.
        if single.eventType == EventType.REACTION:

            print 'single reaction', single
            r = single.drawR( single.dt )

            try:
                self.propagateSingle( single, r )
            except NoSpace:
                # just reject the move, still try reaction.
                print 'single reaction; pre-reaction propagation failed'
                self.rejectedMoves += 1

            try:
                self.shellMatrix.remove( single )
                self.fireSingleReaction( single )
            except NoSpace:
                self.shellMatrix.add( single, single.pos, single.radius )
                self.rejectedMoves += 1
                #self.updateEvent( self.t, single )
                single.dt = self.smallT
                return single.dt

            single.dt = -numpy.inf  # remove this Single from the Scheduler
            return single.dt

        # If not reaction, propagate.

        D0 = single.getD()

        if D0 == 0:
            # no propagation, just calculate reaction times.
            single.determineNextEvent() 
            return single.dt
        

        # (1) propagate
        #
        # Propagate this particle to the exit point on the shell.
        
        squeezed = False
        try:
            squeezed = self.propagateSingle( single, 
                                             single.getMobilityRadius() )
        except NoSpace:
            print 'single propagation move failed'
            self.rejectedMoves += 1
            squeezed = True

        closestShell, distanceToClosestShell =\
                 self.getClosestShell( single.particle.pos, 
                                       ignore = [ single, ] )

        if squeezed:  
            # If this single was squeezed, update the shell, and just return.
            # Because the squeezer Pair was broken up in the propagateSingle()
            # above, one of the Sinlge in the Pair will step in the next step.
            # It is important for this Single to take non-zero step size,
            # otherwise this Single and the squeezer Pair will burst each
            # other indefinitely.

            self.updateSingle( single, closestShell, distanceToClosestShell,
                               math.sqrt( 6 * single.getD() *
                                          self.smallT ) + single.getMinRadius() )
            print 'squeezed single; shell', single.radius, \
                'dt', single.dt

            return single.dt

            
        # (2) Check shell size disparity.   Check if this Single needs
        #     to burst the closest Single or Pair.

        distanceToClosest = self.distance( single.particle.pos, 
                                           closestShell.pos )
        particleRadius0 = single.getMinRadius()

        #print 'criticalPoint %g, closest shell %s, distance to shell %g' %\
        #(criticalPoint,closestShell,distanceToClosestShell)

        somethingBursted = False
        partnerCandidates = []

        # (2-1) Burst the closest and do pair check with that.

        if isinstance( closestShell, Pair ): # pair

            if distanceToClosestShell < particleRadius0 * 2:

                # If the partner was a Pair, and it is bursted into two
                # singles.  Find the closer one.
                single1, single2 = self.burstPair( closestShell )
                somethingBursted = True
                self.removeEvent( closestShell )
                self.addSingleEvent( single1 )
                self.addSingleEvent( single2 )
                
                candidate1 = closestShell.single1
                candidate2 = closestShell.single2
                D1 = candidate1.getD()
                D2 = candidate2.getD()
                pos0 = single.pos
                pos1 = candidate1.pos
                pos2 = candidate2.pos

                dist01 = self.distance( pos0, pos1 ) # radius?
                dist02 = self.distance( pos0, pos2 )
                dist12 = self.distance( pos1, pos2 )

                particleRadius1 = candidate1.getMinRadius()
                particleRadius2 = candidate2.getMinRadius()
                
                # MTTC = Mean Time to Correlate
                MTTC01 = meanArrivalTime( dist01 - particleRadius0 
                                          - particleRadius1, D0 + D1 )
                MTTC02 = meanArrivalTime( dist02 - particleRadius0 
                                          - particleRadius2, 
                                          D0 + D2 )
                MTTC12 = meanArrivalTime( dist12 - particleRadius1 
                                          - particleRadius2, D1 + D2 )

                if dist01 < dist02 and MTTC01 < MTTC12:
                    partnerCandidates = [ candidate1, candidate2 ]
                elif dist02 < dist01 and MTTC02 < MTTC12:
                    partnerCandidates = [ candidate2, candidate1 ]
                # else, partnerCandidates = empty

        else:  
            # If the closest was a Single, that is a potential partnerCandidate.
            
            realDistance = distanceToClosest - particleRadius0 \
                - closestShell.getMinRadius()
            SHELLSIZE_DISPARITY_FACTOR = 0.5
            sqrtD0 = math.sqrt( D0 ) 
            criticalPoint = SHELLSIZE_DISPARITY_FACTOR * sqrtD0 / \
                ( sqrtD0 + math.sqrt( closestShell.getD() ) ) *\
                realDistance + particleRadius0

            #print closestShell, distanceToClosestShell, criticalPoint

            if distanceToClosestShell < criticalPoint or \
                    distanceToClosestShell < particleRadius0 * 2:
                somethingBursted = True
                try:
                    self.burstSingle( closestShell )
                except NoSpace:
                    self.rejectedMoves += 1
                partnerCandidates = [closestShell,]

        #print partnerCandidates
        # try forming a Pair
        if len( partnerCandidates ) >= 1:

            pair = self.formPair( single, partnerCandidates[0] )

            if pair:
                pair.determineNextEvent()

                print pair, 'dt=', pair.dt, 'type=', pair.eventType
                    
                self.addPairEvent( pair )

                self.removeEvent( partnerCandidates[0] )
                    
                for remainingCandidate in partnerCandidates[1:]:
                    c, d =\
                        self.getClosestShell( remainingCandidate.pos,
                                              ignore = 
                                              [ remainingCandidate, ] )
                    self.updateSingle( remainingCandidate, c, d )
                    self.updateEvent( self.t + remainingCandidate.dt,
                                      remainingCandidate )

                single.dt = -numpy.inf # remove by rescheduling to past.
                return single.dt

            else:
                for remainingCandidate in partnerCandidates:
                    c, d =\
                        self.getClosestShell( remainingCandidate.pos,
                                              ignore = 
                                              [ remainingCandidate, ] )
                    self.updateSingle( remainingCandidate, c, d )
                    self.updateEvent( self.t + remainingCandidate.dt,
                                      remainingCandidate )


        # If this Single bursted something,
        # Recheck closest and closest shell distance.
        if somethingBursted:
            closestShell, distanceToClosestShell =\
                self.getClosestShell( single.pos, ignore = [ single, ] )

        # (3) If a new Pair was not formed, this Single continues.
        #     Determine a new shell size and dt.

        # recheck the closest and distance to it.
        self.updateSingle( single, closestShell, distanceToClosestShell )

        #print 'single shell', single.radius, 'dt', single.dt

        return single.dt


    def updateSingle( self, single, closest, distanceToShell, 
                      minShellSize=0.0 ):

        distanceToClosest = self.distance( single.pos, closest.pos )

        shellSize = single.calculateShellSize( closest, distanceToClosest,
                                               distanceToShell )

        shellSize = min( shellSize, self.getMatrixCellSize(),
                         self.maxShellSize )
        shellSize = max( shellSize, minShellSize )

        single.setRadius( shellSize )
        single.determineNextEvent()

        self.shellMatrix.update( single, single.pos, single.radius )




    def firePair( self, pair ):

        print 'fire:', pair, pair.eventType

        particle1 = pair.single1.particle
        particle2 = pair.single2.particle
        species1 = particle1.species
        species2 = particle2.species
        particleRadius1 = species1.radius
        particleRadius2 = species2.radius
        
        oldInterParticle = particle2.pos - particle1.pos

        oldCoM = self.applyBoundary( pair.getCoM() )

        # Three cases:
        #  0. Reaction
        #  1. Escaping through a_r.
        #  2. Escaping through a_R.
        #  3. Single reaction 

        # First handle single reaction case.
        if pair.eventType == 3:

            reactingsingle = pair.reactingsingle

            print 'pair: single reaction', reactingsingle

            if reactingsingle == pair.single1:
                theothersingle = pair.single2
            else:
                theothersingle = pair.single1

            self.burstPair( pair )

            self.addSingleEvent( theothersingle )

            try:
                self.shellMatrix.remove( reactingsingle )
                self.fireSingleReaction( reactingsingle )
            except NoSpace:
                self.shellMatrix.add( reactingsingle, 
                                      reactingsingle.pos, 
                                      reactingsingle.radius )
                self.rejectedMoves += 1
                reactingsingle.dt = 0
                self.addSingleEvent( reactingsingle )

            pair.dt = -numpy.inf
            return pair.dt
        


        #
        # 0. Reaction
        #
        if pair.eventType == EventType.REACTION:

            print 'reaction'

            if len( pair.rt.products ) == 1:
                
                species3 = pair.rt.products[0]

                rnd = numpy.random.uniform( size=2 )

                # calculate new R
            
                r_R = pair.drawR_single( pair.dt, pair.a_R )
            
                displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R
                
                assert self.distance( oldCoM, newCoM ) + species3.radius <\
                    pair.radius

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
                      'num products >= 2 not supported.'

            self.shellMatrix.remove( pair )

            pair.dt = -numpy.inf
            return pair.dt


        #
        # Escape 
        #

        r0 = self.distance( particle1.pos, particle2.pos )

        # 1 Escaping through a_r.
        if pair.eventType == EventType.ESCAPE:

            print 'escape r'

            print 'r0 = ', r0, 'dt = ', pair.dt, pair.pgf.dump()
            
            for i in range(100):

                rnd = numpy.random.uniform( size=4 )

                # calculate new R
            
                r_R = pair.drawR_single( pair.dt, pair.a_R )
                
                displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R

                # calculate new r
                theta_r = pair.drawTheta_pair( rnd[2], pair.a_r, r0, pair.dt )
                phi_r = rnd[3] * 2 * Pi
                newInterParticleS = numpy.array( [ pair.a_r, theta_r, phi_r ] )
                newInterParticle = sphericalToCartesian( newInterParticleS )
                
                newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                      oldInterParticle )
                newpos1 = self.applyBoundary( newpos1 )
                newpos2 = self.applyBoundary( newpos2 )

                if not pair.squeezed or \
                        ( self.checkOverlap( newpos1, particleRadius1,
                                             ignore = [ particle1, 
                                                        particle2 ] ) and \
                              self.checkOverlap( newpos2, particleRadius2,
                                                 ignore = [ particle1, 
                                                            particle2 ] ) ):
                    break

                else:   # overlap check failed
                    self.rejectedMoves += 1
                    print '%s:ESCAPE_r: rejected move. redrawing..' % pair
            else:
                print 'redrawing limit reached.  hanging up..'
                raise RuntimeError,\
                      'redrawing limit reached under squeezing in Pair.ESCAPE_r'


        # 2 escaping through a_R.
        elif pair.eventType == 2:

            print 'escape R'

            for i in range(1000):
                
                rnd = numpy.random.uniform( size = 4 )

                # calculate new r
                print 'r0 = ', r0, 'dt = ', pair.dt, pair.pgf.dump()
                r = pair.drawR_pair( r0, pair.dt, pair.a_r )
                print 'new r = ', r
                #assert r >= pair.sigma
            
                theta_r = pair.drawTheta_pair( rnd[0], r, r0, pair.dt )
                phi_r = rnd[1] * 2*Pi
                newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
                newInterParticle = sphericalToCartesian( newInterParticleS )
                
                # calculate new R
                displacement_R_S = [ pair.a_R, rnd[2] * Pi, rnd[3] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
            
                newCoM = oldCoM + displacement_R
                
                #print 'COM', oldCoM, newCoM

                newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                      oldInterParticle )
                newpos1 = self.applyBoundary( newpos1 )
                newpos2 = self.applyBoundary( newpos2 )

                if not pair.squeezed or \
                        ( self.checkOverlap( newpos1, particleRadius1,
                                             ignore = [ particle1, 
                                                        particle2 ] ) and \
                              self.checkOverlap( newpos2, particleRadius2,
                                                 ignore = [ particle1, 
                                                            particle2 ]) ):
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

        # this has to be done before the following excludeVolume()
        self.shellMatrix.remove( pair )

        assert self.distance( newpos1, newpos2 ) >= pair.sigma

        self.moveParticle( particle1, newpos1 )
        self.moveParticle( particle2, newpos2 )

        if pair.squeezed:
            # make sure displaced particles don't intrude squeezer shells.
            self.excludeSingleVolume( newpos1, particleRadius1 )
            self.excludeSingleVolume( newpos2, particleRadius2 )


        assert self.checkOverlap( newpos1, particleRadius1,
                                  ignore = [ particle1, particle2 ] )
        assert self.checkOverlap( newpos2, particleRadius2,
                                  ignore = [ particle1, particle2 ] )

        single1, single2 = pair.single1, pair.single2

        single1.initialize( self.t )
        single2.initialize( self.t )
            
        self.addSingleEvent( single1 )
        self.addSingleEvent( single2 )

        self.shellMatrix.add( single1, single1.pos, single1.radius )
        self.shellMatrix.add( single2, single2.pos, single2.radius )

        pair.dt = -numpy.inf
        return pair.dt

    def fireMulti( self, multi ):
        print 'fire', multi

        bdsim = BDSimulator()
        #for single in multi.singleList:
        #    bdsim.



    def burstSingle( self, single ):

        #print 'b', self, 't ', t, 'last ', self.lastTime, 'dt ', self.dt
        assert self.t >= single.lastTime
        assert self.t <= single.lastTime + single.dt
        assert single.radius >= single.getMinRadius()

        dt = self.t - single.lastTime

        particleRadius = single.getMinRadius()
        oldpos = single.particle.pos.copy()

        for i in range( 100 ):
            
            r = single.drawR( dt )

            displacement = single.drawDisplacement( r )
            
            newpos = oldpos + displacement
            
            if self.checkOverlap( newpos, particleRadius,
                                  ignore = [ single.particle, ] ):
                break
            
        else:
            #single.particle.pos = oldpos
            #self.moveParticle( single.particle, oldpos )
            single.initialize( self.t )
            self.shellMatrix.update( single, single.pos, single.radius )
            raise NoSpace()

        single.initialize( self.t )
        newpos = self.applyBoundary( newpos )
        self.moveParticle( single.particle, newpos )

        self.shellMatrix.update( single, single.pos, single.radius )
        self.updateEvent( self.t, single )


    def breakUpPair( self, pair ):

        assert self.t >= pair.lastTime

        dt = self.t - pair.lastTime 

        if dt != 0.0:

            particle1 = pair.single1.particle
            particle2 = pair.single2.particle

            particleRadius1 = particle1.species.radius
            particleRadius2 = particle2.species.radius
            
            oldpos1 = particle1.pos.copy()
            oldpos2 = particle2.pos.copy()

            oldInterParticle = particle2.pos - particle1.pos
            oldCoM = pair.getCoM()
            r0 = pair.distance( particle1.pos, particle2.pos )
            
            for i in range(100):
                rnd = numpy.random.uniform( size = 4 )

                # calculate new CoM
                r_R = pair.drawR_single( dt, pair.a_R )
            
                displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R
            
                # calculate new interparticle
                r_r = pair.drawR_pair( r0, dt, pair.a_r )
                theta_r = pair.drawTheta_pair( rnd[2], r_r, r0, dt )
                phi_r = rnd[3] * 2 * Pi
                newInterParticleS = numpy.array( [ r_r, theta_r, phi_r ] )
                newInterParticle = sphericalToCartesian( newInterParticleS )

                newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                      oldInterParticle )
                if not pair.squeezed or \
                        ( self.checkOverlap( newpos1, particleRadius1,
                                             ignore = [ particle1, 
                                                        particle2 ] ) and \
                              self.checkOverlap( newpos2, particleRadius2,
                                                 ignore = [ particle1, 
                                                            particle2 ] ) ):
                    pair.checkNewpos( newpos1, newpos2, oldCoM )

                    newpos1 = self.applyBoundary( newpos1 )
                    newpos2 = self.applyBoundary( newpos2 )
                    self.moveParticle( particle1, newpos1 )
                    self.moveParticle( particle2, newpos2 )
                    break

            else:
                print 'redrawing limit reached.  giving up displacement.'
                self.rejectedMoves += 1

        return ( pair.single1, pair.single2 )


    def burstPair( self, pair ):
        print 'burst', pair

        single1, single2 = self.breakUpPair( pair )
        single1.initialize( self.t )
        single2.initialize( self.t )
        
        self.shellMatrix.remove( pair )
        self.shellMatrix.add( single1, single1.pos, single1.radius )
        self.shellMatrix.add( single2, single2.pos, single2.radius )


        #self.removeEvent( pair )
        #self.addSingleEvent( single1 )
        #self.addSingleEvent( single2 )
        return single1, single2


    def formPair( self, single1, single2 ):

        assert single1.isReset()
        assert single2.isReset()

        # Then, check if this pair of singles meets the pair formation
        # criteria defined in self.checkPairFormationCriteria().
        
        com = calculatePairCoM( single1.pos, single2.pos,\
                                single1.getD(), single2.getD(),\
                                self.getWorldSize() )
        com = self.applyBoundary( com )
        pairClosest, pairClosestShellDistance =\
                     self.getClosestShell( com, ignore = ( single1, single2 ) )
        
        particleRadius1 = single1.getMinRadius()
        particleRadius2 = single2.getMinRadius()
        particleRadius12 = particleRadius1 + particleRadius2

        shellSize = self.checkPairFormationCriteria( single1, single2,
                                                     pairClosest,
                                                     pairClosestShellDistance )
        if shellSize <= 0.0:  # Pair not formed
            return None

        pair = self.createPair( single1, single2 )

        # Squeezed; Pair must be formed but shell size bigger than given space.
        if shellSize > pairClosestShellDistance:
            print 'squeezed shell=', shellSize, 'closest shell distance =',\
                pairClosestShellDistance
            pair.squeezed = True
            self.squeezed += 1

        pair.setRadius( shellSize )
        self.shellMatrix.add( pair, pair.pos, pair.radius )

        pairDistance = self.distance( single1.pos, single2.pos )
        print 'Pair formed: ', pair, 'pair distance', pairDistance,\
              'shell size=', pair.radius,\
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

        particleRadius1 = species1.radius
        particleRadius2 = species2.radius
        particleRadius12 = particleRadius1 + particleRadius2

        pos1, pos2 = single1.pos, single2.pos
        pairDistance = self.distance( pos1, pos2 )

        # pairGap = real distance including radii
        pairGap = pairDistance - particleRadius12
        assert pairGap >= 0, 'pairGap between %s and %s = %g < 0' \
            % ( single1, single2, pairGap )


        minShellSize = max( pairDistance * D1 / D12 + particleRadius1,
                            pairDistance * D2 / D12 + particleRadius2 )

        # consider both D_IV and D_CoM?
        #tau = particleRadius12 * particleRadius12 / D12
        #shellSizeMargin = math.sqrt( 6 * D12 * self.smallT )
        shellSizeMargin = particleRadius12 / 2
        #shellSizeMargin = 5e-1 * ( D1 * particleRadius1 + D2 * particleRadius2 ) / D12
        #print 'margin', shellSizeMargin

        minShellSizeWithMargin = minShellSize + shellSizeMargin

        maxShellSize = min( self.getMatrixCellSize(), self.maxShellSize )


        # 0. Shell cannot be larger than max shell size or sim cell size.
        if minShellSizeWithMargin > maxShellSize:
            return -0.0

        # 1. Squeezed?
        if closestShellDistance <= minShellSizeWithMargin:
            print 'closest shell < minShellSize w/ margin; %g, %g' % \
                ( closestShellDistance, minShellSizeWithMargin )

            # there is enough pairGap; singles can do it better.
            if pairGap >= shellSizeMargin:
                print 'squeezed, but there is enough pairGap.  Pair not formed.'
                return -0.0
            else:
                # real squeezing
                print 'squeezed; pairGap < shellSizeMargin'
                assert minShellSizeWithMargin < \
                    min( self.maxShellSize, self.getMatrixCellSize() )
                return minShellSizeWithMargin

        # 2. Check if a Pair is better than two Singles.
        closestShell = closest.radius
        closestPos = closest.pos
        singleMobility = min( pairDistance - particleRadius12,
                              self.distance( pos1, closestPos )
                              - closestShell - particleRadius1,
                              self.distance( pos2, closestPos )
                              - closestShell - particleRadius2 )
        
        pairMobility = min( closestShellDistance, maxShellSize ) - minShellSize
        if singleMobility >= pairMobility:
            #print 'singleMobility %g >= pairMobility %g' %\
            #      (singleMobility, pairMobility)
            return -0.0
        

        # 3. If not squeezed, and Singles are not better, then this couple
        #    must be a Pair.

        #FIXME: dummy?
        shellSize = minShellSize + ( closestShellDistance - minShellSize ) * .5
        shellSize = max( shellSize, minShellSizeWithMargin )

        shellSize = min( shellSize, self.getMatrixCellSize(), 
                         self.maxShellSize )

        return shellSize
    


    '''
    Find closest n shells.

    This method returns a tuple ( neighbors, distances ).
    '''

    def getNeighborShells( self, pos, n=None ):

        return self.shellMatrix.getNeighbors( pos, n, dummy=DummySingle() )


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

                #assert not closest in ignore
                return closest, distance

        # default case: none left.
        return DummySingle(), numpy.inf


    '''
    Get neighbors simply by distance.

    This method does not take into account of particle radius or shell size.
    This method picks top n neighbors simply by the distance between given
    position pos to positions of objects (either Singles or Pairs) around.

    This method returns a tuple ( neighbors, distances ).

    def getNeighbors( self, pos, n=None ):
        return self.shellMatrix.getNeighbors( pos, n )

    def getNeighbors( self, pos, n=None ):

        objMatrix = self.matrix

        size = objMatrix.size
        if not n:
            n = size 

        distances = self.distanceSqArray( objMatrix.positions, pos )

        topargs = distances.argsort()[:n]
        distances = distances.take( topargs )
        distances = numpy.sqrt( distances )
        neighbors = [ objMatrix.objList[arg] for arg in topargs ]

        return neighbors, distances

    def getClosestNeighbor( self, pos, ignore=[] ):

        neighbors, distances = self.getNeighbors( pos, len( ignore ) + 1 )

        for i in range( len( neighbors ) ): 
            if neighbors[i] not in ignore:
                closest, distance = neighbors[i], distances[i]

                #assert not closest in ignore
                return closest, distance

        # default case: none left.
        return DummySingle(), numpy.inf
    '''


    #
    # consistency checkers
    #
    
    def checkShell( self, obj ):
        closest, distance = self.getClosestShell( obj.pos, [obj,] )
        radius = obj.radius

        if radius > self.getMatrixCellSize():
            raise RuntimeError, '%s shell size larger than simulator cell size'

        if radius > self.maxShellSize:
            raise RuntimeError, '%s shell size larger than maxShellSize'

        if distance - radius < 0.0:
            if ( isinstance( obj, Pair ) and obj.squeezed ) or \
                   ( isinstance( closest, Pair ) and closest.squeezed ):
                print '%s overlaps with %s.  ignoring because squeezed.' \
                          % ( str( obj ), str( closest ) )
            else:
                raise RuntimeError,\
                      '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                      % ( str( obj ), str( closest ), radius, distance,\
                          distance - radius )


    def checkShellForAll( self ):

        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            self.checkShell( obj )


    def checkEventStoichiometry( self ):

        population = 0
        for species in self.speciesList.values():
            population += species.pool.size
        
        eventPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            eventPopulation += obj.multiplicity

        if population != eventPopulation:
            raise RuntimeError, 'population %d != eventPopulation %d' %\
                  ( population, eventPopulation )

    def checkShellMatrix( self ):

        if self.worldSize != self.shellMatrix.worldSize:
            raise RuntimeError,\
                'self.worldSize != self.shellMatrix.worldSize'

        if self.scheduler.getSize() != self.shellMatrix.size:
            raise RuntimeError,\
                'self.scheduler.getSize() != self.shellMatrix.size'
        
        self.shellMatrix.check()

        for key in self.shellMatrix.keyList:
            pos, radius = self.shellMatrix.get( key )
            if ( key.pos - pos ).sum() != 0:
                raise RuntimeError, 'shellMatrix positions consistency broken'
            if key.radius != radius:
                raise RuntimeError, 'shellMatrix radii consistency broken'



    def check( self ):

        GFRDSimulatorBase.check( self )

        assert self.t >= 0.0
        assert self.dt >= 0.0

        self.checkShellMatrix()

        self.checkEventStoichiometry()
        
        self.checkShellForAll()




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
            print i, event.getTime(), event.getArg(), event.getArg().pos



