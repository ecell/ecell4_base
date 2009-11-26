#!/usr/env python


from weakref import ref
import math

import numpy

from _gfrd import EventScheduler, FirstPassageGreensFunction, FirstPassagePairGreensFunction, FirstPassageNoCollisionPairGreensFunction, BasicPairGreensFunction, FreePairGreensFunction, EventType

from utils import *
from surface import CuboidalSurface

from gfrdbase import *
from bd import BDSimulatorCoreBase

from numpy.random import uniform

import logging

log = logging.getLogger( 'epdp' )

SAFETY = 1.0 + 1e-5

MULTI_SHELL_FACTOR = 0.05
SINGLE_SHELL_FACTOR = 0.1

class Delegate( object ):
    def __init__( self, obj, method ):
        self.ref = ref( obj )
        self.method = method

    def __call__( self, *arg ):
        return self.method( self.ref(), *arg )


class MultiBDCore( BDSimulatorCoreBase ):
    '''
    Used internally by Multi.
    '''
    def __init__( self, main, multi ):

        BDSimulatorCoreBase.__init__( self, main )

        # this has to be ref, not proxy, since it is used for comparison.
        self.multiref = ref( multi )

        self.particleMatrix = ObjectMatrix()
        self.particleMatrix.setWorldSize( self.main.worldSize )

        self.shellMatrix = ObjectMatrix()
        self.shellMatrix.setWorldSize( self.main.worldSize )

        self.escaped = False

    def updateParticle( self, particle, pos ):

        self.particleMatrix.update( particle, pos, particle.radius )
        self.main.updateOnParticleMatrix( particle, pos )

    def initialize( self ):

        BDSimulatorCoreBase.initialize( self )

    def step( self ):
        self.escaped = False
        BDSimulatorCoreBase.step( self )

    def addParticle( self, particle ):
        self.addToParticleList( particle )
        self.particleMatrix.update( particle,
                                    particle.pos, particle.radius )

    def removeParticle( self, particle ):
        self.main.removeParticle( particle )
        self.removeFromParticleList( particle )
        self.particleMatrix.remove( particle )

    def createParticle( self, species, pos ):
        particle = self.main.createParticle( species, pos )
        self.addParticle( particle )
        return particle

    def moveParticle( self, particle, pos ):
        particle.pos = pos
        self.updateParticle( particle, pos )

    def clearVolume( self, pos, radius, ignore=[] ):
        if not self.withinShell( pos, radius ):
            self.escaped = True
            self.clearOuterVolume( pos, radius, ignore )

    def clearOuterVolume( self, pos, radius, ignore=[] ):
        self.main.clearVolume( pos, radius, ignore=[self.multiref(),] )
        if not self.main.checkOverlap( pos, radius, ignore ):
            raise NoSpace()

    def withinShell( self, pos, radius ):
        n, _ = self.shellMatrix.getNeighborsWithinRadiusNoSort( pos, - radius )
        return n
        
    def checkOverlap( self, pos, radius, ignore=[] ):
        n, _ = self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        n = list( n )
        for particle in ignore:
            if particle in n:
                n.remove( particle )

        return not n

    def getParticlesWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        neighbors, _ = \
            self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
        return [ n for n in neighbors if n not in ignore ]

    def check( self ):
        BDSimulatorCoreBase.check( self )

        # shells are contiguous
        for shellID in self.multiref().shellIDList:
            pos, radius = self.shellMatrix.get(shellID)
            _, d = self.shellMatrix.getNeighborsCyclicNoSort(pos)
            d.sort()
            assert d[1] - radius < 0.0,\
                'shells of %s are not contiguous.' % str(self.multi())

        # all particles within the shell.
        for p in self.particleList:
            assert self.withinShell( p.pos, p.species.radius ),\
                'not all particles within the shell.'


class Single( object ):
    def __init__( self, particle, reactiontypes ):
        self.multiplicity = 1

        self.particle = particle
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.shellID = (self, 0)
        self.shellIDList = [self.shellID]

        self.eventID = None

        self.updatek_tot()

    def getD( self ):
        return self.particle.species.D

    D = property( getD )

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


    '''
    Reset the Single.

    Radius (shell size) is shrunken to the actual radius of the particle.
    self.dt is reset to 0.0.  Do not forget to reschedule this Single
    after calling this method.
    '''

    def reset( self ):

        self.dt = 0.0
        self.eventType = EventType.ESCAPE


    def isReset( self ):

        return self.dt == 0.0 and self.eventType == EventType.ESCAPE
        

    def drawR( self, dt, shellSize):
        assert dt >= 0

        a = shellSize - self.particle.species.radius

        rnd = uniform()
        gf = FirstPassageGreensFunction( self.particle.species.D )
        gf.seta(a)

        try:
            r = gf.drawR(rnd , dt)
        except Exception, e:
            raise Exception, 'gf.drawR failed; %s; rnd=%g, t=%g, %s' %\
                (str(e), rnd, dt, gf.dump())

        return r




    def drawReactionTime( self ):
        
        if self.k_tot == 0:
            return numpy.inf

        if self.k_tot == numpy.inf:
            return 0.0

        rnd = uniform()
        dt = ( 1.0 / self.k_tot ) * math.log( 1.0 / rnd )

        return dt


    def drawEscapeTime(self, a):
        
        gf = FirstPassageGreensFunction( self.particle.species.D )
        gf.seta(a)

        try:
            rnd = uniform()
            return gf.drawTime( rnd )
        except Exception, e:
            raise Exception, 'gf.drawTime() failed; %s; rnd=%g, %s' %\
                ( str( e ), rnd, gf.dump() )


    def updatek_tot( self ):

        self.k_tot = 0

        if not self.reactiontypes:
            return

        for rt in self.reactiontypes:
            self.k_tot += rt.k


    def drawReactionRule( self ):

        k_array = [ rt.k for rt in self.reactiontypes ]
        k_array = numpy.add.accumulate( k_array )
        k_max = k_array[-1]

        rnd = uniform()
        i = numpy.searchsorted( k_array, rnd * k_max )

        return self.reactiontypes[i]


    def check( self ):
        pass

    def __repr__( self ):
        return 'Single[%s: eventID=%s]' % ( self.particle, self.eventID )


# def calculatePairCoM( pos1, pos2, D1, D2, worldSize ):
#     '''
#     Calculate and return the center-of-mass of a Pair.
#     '''
#     #pos2t = cyclicTranspose( pos2, pos1, worldSize )
#     pos2t = cyclic_transpose( pos2, pos1, worldSize )
#     return ( ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 ) ) % worldSize


class Pair( object ):
    
    # CUTOFF_FACTOR is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__(self, single1, single2, rt):
        self.multiplicity = 2

        # Order single1 and single2 so that D1 < D2.
        if single1.particle.species.D <= single2.particle.species.D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.rt = rt

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        self.D1, self.D2 = particle1.species.D, particle2.species.D

        self.D_tot = self.D1 + self.D2
        self.D_R = ( self.D1 * self.D2 ) / self.D_tot

        #self.minRadius = max( particle1.species.radius,
        #                      particle2.species.radius )
        self.sigma = particle1.species.radius + particle2.species.radius

        self.eventID = None

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.shellID = (self, 0)
        self.shellIDList = [self.shellID]

    def __del__( self ):
        if __debug__:
            log.debug( 'del %s' % str( self ) )

    def initialize( self, t ):

        self.lastTime = t
        self.dt = 0
        self.eventType = None

    def getD( self ):
        return self.D_tot #FIXME: is this correct?

    def choosePairGreensFunction( self, r0, t ):
        distanceFromSigma = r0 - self.sigma
        distanceFromShell = self.a_r - r0

        thresholdDistance = Pair.CUTOFF_FACTOR * \
            math.sqrt( 6.0 * self.D_tot * t )

        if distanceFromSigma < thresholdDistance:
        
            if distanceFromShell < thresholdDistance:
                # near both a and sigma;
                # use FirstPassagePairGreensFunction
                if __debug__:
                    log.debug( 'GF: normal' )
                pgf = FirstPassagePairGreensFunction( self.D_tot, 
                                                      self.rt.k, self.sigma )

                return pgf
            else:
                # near sigma; use BasicPairGreensFunction
                if __debug__:
                    log.debug( 'GF: only sigma' )
                pgf = BasicPairGreensFunction( self.D_tot, self.rt.k, 
                                               self.sigma )
                return pgf
        else:
            if distanceFromShell < thresholdDistance:
                # near a;
                if __debug__:
                    log.debug( 'GF: only a' )
                pgf = FirstPassageNoCollisionPairGreensFunction( self.D_tot )
                return pgf
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug( 'GF: free' )
                pgf = FreePairGreensFunction( self.D_tot )
                return pgf

    def newPositions( self, CoM, newInterParticle, oldInterParticle ):
        '''
        Calculate new positions of the pair particles using
        a new center-of-mass, a new inter-particle vector, and
        an old inter-particle vector.
        '''
        #FIXME: need better handling of angles near zero and pi.

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


    def drawTime_single( self, sgf ):
        rnd = uniform()
        return sgf.drawTime( rnd )

    def drawTime_pair( self, pgf, r0 ):
        rnd = uniform()
        #print 'r0 = ', r0, ', rnd = ', rnd[1],\
        #    pgf.dump()
        return pgf.drawTime( rnd, r0 )

    def drawEventType( self, pgf, r0, t ):
        rnd = uniform()
        return pgf.drawEventType( rnd, r0, t )

    def drawR_single( self, sgf, t ):
        rnd = uniform()
        try:
            r = sgf.drawR( rnd, t )
            while r > self.a_R: # redraw; shouldn't happen often
                if __debug__:
                    log.info( 'drawR_single: redraw' )
                rnd = uniform()
                r = sgf.drawR( rnd, t )
        except Exception, e:
            raise Exception,\
                'gf.drawR_single() failed; %s; rnd= %g, t= %g, %s' %\
                ( str( e ), rnd, t, sgf.dump() )

        return r

    def drawR_pair( self, r0, t, a ):
        '''
        Draw r for the pair inter-particle vector.
        '''
        gf = self.choosePairGreensFunction( r0, t )

        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        rnd = uniform()
        try:
            r = gf.drawR( rnd, r0, t )
            # redraw; shouldn't happen often
            while r >= self.a_r or r <= self.sigma: 
                if __debug__:
                    log.info( 'drawR_pair: redraw' )
                #self.sim.rejectedMoves += 1  #FIXME:
                rnd = uniform()
                r = gf.drawR( rnd, r0, t )
        except Exception, e:
            raise Exception,\
                'gf.drawR_pair() failed; %s; rnd= %g, r0= %g, t= %g, %s' %\
                ( str( e ), rnd, r0, t, gf.dump() )


        return r

    def drawTheta_pair( self, rnd, r, r0, t, a ):
        '''
        Draw theta for the pair inter-particle vector.
        '''
        gf = self.choosePairGreensFunction( r0, t )

        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        try:
            theta = gf.drawTheta( rnd, r, r0, t )
        except Exception, e:
            raise Exception,\
                'gf.drawTheta() failed; %s; rnd= %g, r= %g, r0= %g, t=%g, %s' %\
                ( str( e ), rnd, r, r0, t, gf.dump() )

        return theta

    def check( self ):
        pass

    def __repr__( self ):
        return 'Pair[%s, %s: eventID=%s]' % (
            self.single1.particle,
            self.single2.particle,
            self.eventID )


class Multi( object ):
    def __init__( self, main ):
        self.shellIDList = []
        self.eventID = None

        self.sim = MultiBDCore( main, self )

    def initialize( self, t ):
        self.lastTime = t
        self.startTime = t

        self.sim.initialize() # ??

    def getDt( self ):
        return self.sim.dt

    dt = property( getDt )

    def getMultiplicity( self ):
        return len( self.sim.particleList )

    multiplicity = property( getMultiplicity )

    def addParticle( self, particle ):
        self.sim.addParticle( particle )

    def addShell( self, pos, radius ):
        shellID = (self,len(self.shellIDList))
        self.shellIDList.append(shellID)
        self.sim.main.shellMatrix.update(shellID, pos, radius)
        self.sim.shellMatrix.update(shellID, pos, radius)

    def check( self ):
        self.sim.check()

    def __repr__( self ):
        return 'Multi[%s: eventID=%s]' % (
            ', '.join( repr( p ) for p in self.sim.particleList ),
            self.eventID )


class DummySingle( object ):
    def __init__( self ):
        self.radius = 0.0

    def getMinRadius( self ):
        return 0.0

    def getD( self ):
        return 0.0

    def __getPos( self ):
        return NOWHERE

    pos = property( __getPos )

    def __repr__( self ):
        return 'DummySingle()'


class EGFRDSimulator( ParticleSimulatorBase ):
    def __init__( self ):
        self.shellMatrix = ObjectMatrix()

        ParticleSimulatorBase.__init__( self )

        self.MULTI_SHELL_FACTOR = 0.05
        self.SINGLE_SHELL_FACTOR = 0.1

        self.isDirty = True
        self.scheduler = EventScheduler()

        self.smallT = 1e-8  # FIXME: is this ok?

        self.userMaxShellSize = numpy.inf

        self.reset()

    def setWorldSize( self, size ):
        if isinstance( size, list ) or isinstance( size, tuple ):
            size = numpy.array( size )

        ParticleSimulatorBase.setWorldSize( self, size )
        self.shellMatrix.setWorldSize( size )

    def setMatrixSize( self, size ):
        ParticleSimulatorBase.setMatrixSize( self, size )

        if self.maxMatrixSize == 0:
            self.shellMatrix.setMatrixSize( size )
        else:
            self.shellMatrix.setMatrixSize( max( size, self.maxMatrixSize ) )

        self.shellMatrix.setMatrixSize( size )

    def getMatrixCellSize( self ):
        return self.shellMatrix.cellSize

    def getNextTime( self ):
        if self.scheduler.getSize() == 0:
            return self.t

        return self.scheduler.getTopTime()

    def setUserMaxShellSize( self, size ):
        self.userMaxShellSize = size

    def getUserMaxShellSize( self ):
        return self.userMaxShellSize

    def getMaxShellSize( self ):
        return min( self.getMatrixCellSize() * .5 / SAFETY,
                    self.userMaxShellSize )

    def reset( self ):
        self.t = 0.0
        self.dt = 0.0
        self.stepCounter = 0
        self.zeroSteps = 0
        self.rejectedMoves = 0
        self.reactionEvents = 0
        self.lastEvent = None
        self.lastReaction = None

        self.isDirty = True
        #self.initialize()

    def initialize( self ):
        ParticleSimulatorBase.initialize( self )

        self.scheduler.clear()
        self.shellMatrix.clear()

        singles = []
        for p in self.particleMatrix:
            particle = p[0]
            singles.append(self.createSingle(particle))
        for single in singles:
            self.addSingleEvent( single )

        self.isDirty = False

    def stop( self, t ):
        if __debug__:
            log.info( 'stop at %g' % t )

        if self.t == t:
            return

        if t >= self.scheduler.getTopEvent().getTime():
            raise RuntimeError, 'Stop time >= next event time.'

        if t < self.t:
            raise RuntimeError, 'Stop time < current time.'

        self.t = t
        
        scheduler = self.scheduler
        
        nonSingleList = []

        # first burst all Singles.
        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getArg()
            if isinstance( obj, Pair ) or isinstance( obj, Multi ):
                nonSingleList.append( obj )
            elif isinstance( obj, Single ):
                if __debug__:
                    log.debug( 'burst %s, lastTime= %g' % 
                           ( str( obj ), obj.lastTime ) )
                self.burstSingle( obj )
            else:
                assert False, 'do not reach here'


        # then burst all Pairs and Multis.
        if __debug__:
            log.debug( 'burst %s' % nonSingleList )
        self.burstObjs( nonSingleList )

        self.dt = 0.0

    def step( self ):
        self.lastReaction = None

        if self.isDirty:
            self.initialize()
            
        #if self.stepCounter % 100 == 0:
        #self.check()
        
        self.stepCounter += 1

        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getArg()

        if __debug__:
            log.info( '\n%d: t=%g dt=%g\nevent=%s reactions=%d rejectedmoves=%d' 
                      % ( self.stepCounter, self.t, self.dt, self.lastEvent, 
                          self.reactionEvents, self.rejectedMoves ) )
        
        self.scheduler.step()

        nextTime = self.scheduler.getTopTime()
        self.dt = nextTime - self.t


        # assert if not too many successive dt=0 steps occur.
        if __debug__:
            if self.dt == 0:
                self.zeroSteps += 1
                if self.zeroSteps >= max( self.scheduler.getSize() * 3, 10 ):
                    raise RuntimeError, 'too many dt=zero steps.  simulator halted?'
            else:
                self.zeroSteps = 0


        assert self.scheduler.getSize() != 0

    def createSingle( self, particle ):
        rt = self.getReactionRule1( particle.species )
        single = Single( particle, rt )
        single.initialize( self.t )
        self.moveShell(single.shellID, particle.pos,\
                           single.particle.species.radius)
 
        return single


    def determineSingleEvent(self, single, t, shellSize):

        a = shellSize - single.particle.species.radius
        if single.getD() == 0:
            firstPassageTime = numpy.inf
        else:
            firstPassageTime = single.drawEscapeTime(a)
            
        reactionTime = single.drawReactionTime()

        if firstPassageTime <= reactionTime:
            single.dt = firstPassageTime
            single.eventType = EventType.ESCAPE
        else:
            single.dt = reactionTime
            single.eventType = EventType.REACTION

        single.lastTime = t

        if __debug__:
            log.info( "determineSingleEvent: type=%s, dt=%g" %\
                          ( single.eventType, single.dt ) )


    def createPair(self, single1, single2, shellSize):
        assert single1.dt == 0.0
        assert single2.dt == 0.0

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.getReactionRule2( species1, species2 )[ 0 ]

        pair = Pair( single1, single2, rt ) 
        pair.initialize( self.t )

        self.moveShell(pair.shellID, self.calculatePairCoM(pair),\
                           shellSize)

        return pair

    def calculatePairCoM(self, pair):
        return calculate_pair_CoM(pair.single1.particle.pos,
                                  pair.single2.particle.pos,
                                  pair.D1, pair.D2, self.worldSize)


    def determinePairEvent(self, pair, t, pos1, pos2, shellSize):
        pair.lastTime = t

        single1 = pair.single1
        single2 = pair.single2
        radius1 = single1.particle.species.radius
        radius2 = single2.particle.species.radius

        D1 = pair.D1
        D2 = pair.D2

        shellSize /= SAFETY  # FIXME:

        D_tot = D1 + D2
        D_geom = math.sqrt(D1 * D2)

        r0 = self.distance(pos1, pos2)

        assert r0 >= pair.sigma, \
            '%s;  r0 %g < sigma %g' % ( pair, r0, pair.sigma )

        # equalize expected mean t_r and t_R.
        if ((D_geom - D2) * r0) / D_tot + shellSize +\
                math.sqrt(D2 / D1) * (radius1 - shellSize) - radius2 >= 0:
            Da = D1
            Db = D2
            radiusa = radius1
            radiusb = radius2
        else:
            Da = D2
            Db = D1
            radiusa = radius2
            radiusb = radius1


        #aR
        pair.a_R = (D_geom * (Db * (shellSize - radiusa) + \
                               Da * (shellSize - r0 - radiusa))) /\
                               (Da * Da + Da * Db + D_geom * D_tot)

        #ar
        pair.a_r = (D_geom * r0 + D_tot * (shellSize - radiusa)) /\
            (Da + D_geom)

        assert pair.a_R + pair.a_r * Da / D_tot + radius1 >= \
            pair.a_R + pair.a_r * Db / D_tot + radius2

        assert abs( pair.a_R + pair.a_r * Da / D_tot + radiusa - shellSize ) \
            < 1e-12 * shellSize


        if __debug__:
          log.debug( 'a %g, r %g, R %g r0 %g' % 
                  ( shellSize, pair.a_r, pair.a_R, r0 ) )
        if __debug__:
          log.debug( 'tr %g, tR %g' % 
                     ( ( ( pair.a_r - r0 ) / math.sqrt(6 * pair.D_tot))**2,\
                           (pair.a_R / math.sqrt( 6*pair.D_R ))**2 ) )
        assert pair.a_r > 0
        assert pair.a_r > r0, '%g %g' % ( pair.a_r, r0 )
        assert pair.a_R > 0 or ( pair.a_R == 0 and ( D1 == 0 or D2 == 0 ) )

        sgf = FirstPassageGreensFunction(pair.D_R)
        sgf.seta(pair.a_R)


        # draw t_R
        try:
            pair.t_R = pair.drawTime_single( sgf )
        except Exception, e:
            raise Exception, 'sgf.drawTime() failed; %s; %s' %\
                ( str( e ), sgf.dump() )

        pgf = FirstPassagePairGreensFunction( pair.D_tot, 
                                              pair.rt.k, pair.sigma )
        pgf.seta( pair.a_r )

        # draw t_r
        try:
            pair.t_r = pair.drawTime_pair(pgf, r0)
        except Exception, e:
            raise Exception, \
                'pgf.drawTime() failed; %s; r0=%g, %s' % \
                ( str( e ), r0, pgf.dump() )


        # draw t_reaction
        t_reaction1 = pair.single1.drawReactionTime()
        t_reaction2 = pair.single2.drawReactionTime()

        if t_reaction1 < t_reaction2:
            pair.t_single_reaction = t_reaction1
            pair.reactingsingle = pair.single1
        else:
            pair.t_single_reaction = t_reaction2
            pair.reactingsingle = pair.single2

        pair.dt = min( pair.t_R, pair.t_r, pair.t_single_reaction )

        assert pair.dt >= 0
        if __debug__:
            log.debug( 'dt %g, t_R %g, t_r %g' % 
                     ( pair.dt, pair.t_R, pair.t_r ) )

        if pair.dt == pair.t_r:  # type = 0 (REACTION) or 1 (ESCAPE_r)
            try:
                pair.eventType = pair.drawEventType(pgf, r0, pair.t_r)
            except Exception, e:
                raise Exception,\
                    'pgf.drawEventType() failed; %s; r0=%g, %s' %\
                    ( str( e ), r0, pgf.dump() )

        elif pair.dt == pair.t_R: # type = ESCAPE_R (2)
            pair.eventType = 2
        elif pair.dt == pair.t_single_reaction:  # type = single reaction (3)
            pair.eventType = 3 
        else:
            raise AssertionError, "Never get here"




    def createMulti( self ):
        multi = Multi( self )
        return multi

    def moveSingle(self, single, pos, radius):
        self.moveShell(single.shellID, pos, radius)
        self.moveParticle(single.particle, pos)

    def removeFromShellMatrix( self, obj ):
        for shellID in obj.shellIDList:
            self.shellMatrix.remove(shellID)
    
    def getShell(self, shellID):
        return self.shellMatrix.get(shellID)

    def moveShell(self, shellID, pos, radius):
        self.shellMatrix.update(shellID, pos, radius)

    def calculateSingleMobilityRadius(self, single):
        radius = self.shellMatrix.get(single.shellID)[1]
        return radius - single.particle.species.radius

    def addEvent( self, t, func, arg ):
        return self.scheduler.addEvent( t, func, arg )

    def addSingleEvent( self, single ):
        eventID = self.addEvent( self.t + single.dt, 
                                 Delegate( self, EGFRDSimulator.fireSingle ), 
                                 single )
        if __debug__:
            log.info( 'addSingleEvent: #%d (t=%g)' % (
                eventID, self.t + single.dt ) )
        single.eventID = eventID

    def addPairEvent( self, pair ):
        eventID = self.addEvent( self.t + pair.dt, 
                                 Delegate( self, EGFRDSimulator.firePair ), 
                                 pair )
        if __debug__:
            log.info( 'addPairEvent: #%d (t=%g)' % (
                eventID, self.t + pair.dt ) )
        pair.eventID = eventID

    def addMultiEvent( self, multi ):
        eventID = self.addEvent( self.t + multi.dt, 
                                 Delegate( self, EGFRDSimulator.fireMulti ), 
                                 multi )
        if __debug__:
            log.info( 'addMultiEvent: #%d (t=%g)' % (
                eventID, self.t + multi.dt ) )
        multi.eventID = eventID

    def removeEvent( self, event ):
        if __debug__:
            log.info( 'removeEvent: #%d' % event.eventID )
        self.scheduler.removeEvent( event.eventID )

    def updateEvent( self, t, event ):
        if __debug__:
            log.info( 'updateEvent: #%d (t=%g)' % ( event.eventID, t ) )
        self.scheduler.updateEventTime( event.eventID, t )

    def burstObj( self, obj ):
        if __debug__:
            log.info( 'burstObj: bursting %s' % obj )

        if isinstance( obj, Single ):
            self.burstSingle( obj )
            bursted = [obj,]
        elif isinstance( obj, Pair ):  # Pair
            single1, single2 = self.burstPair( obj )
            self.removeEvent( obj )
            self.addSingleEvent( single1 )
            self.addSingleEvent( single2 )
            bursted = [ single1, single2 ]
        else:  # Multi
            bursted = self.burstMulti( obj )
            self.removeEvent( obj )

        if __debug__:
            log.info( 'burstObj: bursted=%s' % bursted )

        return bursted

    def burstObjs( self, objs ):
        bursted = []
        for obj in objs:
            b = self.burstObj( obj )
            bursted.extend( b )

        return bursted

    def clearVolume( self, pos, radius, ignore=[] ):
        neighbors = self.getNeighborsWithinRadiusNoSort( pos, radius, ignore )
        return self.burstObjs( neighbors )

    def burstNonMultis( self, neighbors ):
        bursted = []

        for obj in neighbors:
            if not isinstance( obj, Multi ):
                b = self.burstObj( obj )
                bursted.extend( b )
            else:
                bursted.append( obj )

        return bursted

    def fireSingleReaction( self, single ):
        reactantSpecies = single.particle.species
        oldpos = single.particle.pos.copy()
        
        rt = single.drawReactionRule()

        if len( rt.products ) == 0:
            
            self.removeParticle( single.particle )

            self.lastReaction = Reaction( rt, [single.particle], [] )

            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            if reactantSpecies.radius < productSpecies.radius:
                self.clearVolume( oldpos, productSpecies.radius )

            if not self.checkOverlap( oldpos, productSpecies.radius,
                                      ignore = [ single.particle, ] ):
                if __debug__:
                    log.info( 'no space for product particle.' )
                raise NoSpace()

            self.removeParticle( single.particle )
            newparticle = self.createParticle( productSpecies, oldpos )
            newsingle = self.createSingle( newparticle )
            self.addSingleEvent( newsingle )

            self.lastReaction = Reaction( rt, [single.particle], [newparticle] )

            if __debug__:
                log.info( 'product; %s' % str( newsingle ) )

            
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

            self.clearVolume( oldpos, rad )

            for _ in range( 100 ):
                unitVector = randomUnitVector()
                vector = unitVector * particleRadius12 * ( 1.0 + 1e-7 )
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?

                while 1:
                    newpos1 = oldpos + vector * ( D1 / D12 )
                    newpos2 = oldpos - vector * ( D2 / D12 )
                    self.applyBoundary( newpos1 )
                    self.applyBoundary( newpos2 )

                    if self.distance( newpos1, newpos2 ) >= particleRadius12:
                        break

                    vector *= 1.0 + 1e-7


                # accept the new positions if there is enough space.
                if ( self.checkOverlap( newpos1, particleRadius1,
                                        ignore = [ single.particle, ] ) and
                     self.checkOverlap( newpos2, particleRadius2,
                                        ignore = [ single.particle, ] ) ):
                    break
            else:
                if __debug__:
                    log.info( 'no space for product particles.' )
                raise NoSpace()

            self.removeParticle( single.particle )

            particle1 = self.createParticle( productSpecies1, newpos1 )
            particle2 = self.createParticle( productSpecies2, newpos2 )
            newsingle1 = self.createSingle( particle1 )
            newsingle2 = self.createSingle( particle2 )

            self.addSingleEvent( newsingle1 )
            self.addSingleEvent( newsingle2 )

            self.lastReaction = Reaction( rt, [single.particle], 
                                          [particle1, particle2] )

            if __debug__:
                log.info( 'products; %s %s' % 
                      ( str( newsingle1 ), str( newsingle2 ) ) )

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1

    def propagateSingle( self, single, r ):
        if __debug__:
            log.debug( "single.dt=%f, single.lastTime=%f, self.t=%f" % (
                single.dt, single.lastTime, self.t ) )
        assert abs( single.dt + single.lastTime - self.t ) <= 1e-18 * self.t
        
        displacement = randomVector( r )

        assert abs( length( displacement ) - r ) <= 1e-15 * r
            
        newpos = single.particle.pos + displacement
        self.applyBoundary( newpos )
            
        assert self.checkOverlap( newpos, single.particle.species.radius,\
                                  ignore = [ single.particle, ] )

        single.initialize( self.t )

        self.moveSingle(single, newpos, single.particle.species.radius)

        return newpos


    def fireSingle( self, single ):

        shellSize = self.getShell(single.shellID)[1]

        # Reaction.
        if single.eventType == EventType.REACTION:

            if __debug__:
                log.info( 'single reaction %s' % str( single ) )
            r = single.drawR( single.dt, shellSize )

            self.propagateSingle( single, r )

            try:
                self.removeFromShellMatrix( single )
                self.fireSingleReaction( single )
                return
            except NoSpace:
                if __debug__:
                    log.info( 'single reaction; placing product failed.' )
                self.addToShellMatrix( single )
                self.rejectedMoves += 1
                single.reset()
                self.addSingleEvent(single)
                return

        # Propagate, if not reaction.

        # Handle immobile case first.
        if single.getD() == 0:
            # no propagation, just calculate next reaction time.
            self.determineSingleEvent(single.self.t, shellSize) 
            self.addSingleEvent(single)
            return
        
        # Propagate this particle to the exit point on the shell.
        a = shellSize - single.particle.species.radius
        singlepos = self.propagateSingle(single, a)

        # (2) Clear volume.

        minShell = single.particle.species.radius * ( 1.0 + SINGLE_SHELL_FACTOR )

        closeNeighbors, distances = self.getNeighbors( singlepos, minShell,
                                                       ignore=[single,] )
        if __debug__:
            log.debug( "closeNeighbors: %s" % closeNeighbors )
        # This is a bit tricky, but the last one in closeNeighbors
        # is the closest object to this Single.
        # getNeighbors() returns closeNeighbors within minShell *plus* one.
        closest = closeNeighbors.pop()
        closestShellDistance = distances[-1]

        bursted = []
        
        if closeNeighbors:
            bursted = self.burstNonMultis( closeNeighbors )
            obj = self.formPairOrMulti(single, singlepos, bursted)

            if obj:
                return

            # if nothing was formed, recheck closest and restore shells.
            closest, closestShellDistance = \
                self.getClosestObj( singlepos, ignore = [ single, ] )

        self.updateSingle( single, singlepos, closest, closestShellDistance )

        bursted = uniq( bursted )
        burstedSingles = [ s for s in bursted if isinstance( s, Single ) ]
        self.restoreSingleShells( burstedSingles )
            
        if __debug__:
            log.info( 'single shell %s dt %g.' %\
                          ( self.getShell(single.shellID), single.dt ) )

        self.addSingleEvent(single)
        return

    def restoreSingleShells( self, singles ):
        for single in singles:
            assert single.isReset()
            singlepos, singleradius = self.getShell(single.shellID)
            c, d = self.getClosestObj( singlepos, ignore = [single,] )

            self.updateSingle( single, singlepos, c, d )
            self.updateEvent( self.t + single.dt, single )
            if __debug__:
                log.debug( 'restore shell %s %g dt %g closest %s %g' %
                       ( single, singleradius, single.dt, c, d ) )

    def calculateSingleShellSize( self, single, closest, 
                                  distance, shellDistance ):
        assert isinstance( closest, Single )

        minRadius1 = single.particle.species.radius
        D1 = single.getD()

        if D1 == 0:
            return minRadius1

        D2 = closest.getD()
        minRadius2 = closest.particle.species.radius
        minRadius12 = minRadius1 + minRadius2
        sqrtD1 = math.sqrt( D1 )
            
        shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                         * ( distance - minRadius12 ) + minRadius1,
                         shellDistance )
        shellSize /= SAFETY
        shellSize = max( shellSize, minRadius1 ) # not smaller than the radius

        return shellSize

    def updateSingle( self, single, singlepos, closest, distanceToShell ): 
        if isinstance( closest, Single ):
            closestpos = self.getShell(closest.shellID)[0]
            distanceToClosest = self.distance(singlepos, closestpos)
            shellSize = self.calculateSingleShellSize( single, closest, 
                                                       distanceToClosest,
                                                       distanceToShell )
        else:  # Pair or Multi
            shellSize = distanceToShell / SAFETY
            shellSize = max( shellSize, single.particle.species.radius )

        shellSize = min( shellSize, self.getMaxShellSize() )

        self.determineSingleEvent(single, self.t, shellSize)
        self.moveShell(single.shellID, singlepos, shellSize)

    def firePair( self, pair ):
        assert self.checkObj( pair )

        if __debug__:
            log.info( 'fire: %s eventType %s' % ( pair, pair.eventType ) )

        particle1 = pair.single1.particle
        particle2 = pair.single2.particle
        
        oldInterParticle = particle2.pos - particle1.pos
        oldCoM = self.calculatePairCoM(pair)
        self.applyBoundary( oldCoM )

        # Three cases:
        #  0. Reaction
        #  1. Escaping through a_r.
        #  2. Escaping through a_R.
        #  3. Single reaction 

        # First handle single reaction case.
        if pair.eventType == 3:

            reactingsingle = pair.reactingsingle

            if __debug__:
                log.info( 'pair: single reaction %s' % str( reactingsingle ) )

            if reactingsingle == pair.single1:
                theothersingle = pair.single2
            else:
                theothersingle = pair.single1

            self.burstPair( pair )

            self.addSingleEvent( theothersingle )

            try:
                self.removeFromShellMatrix( reactingsingle )
                self.fireSingleReaction( reactingsingle )
            except NoSpace:
                self.addToShellMatrix( reactingsingle )
                self.rejectedMoves += 1
                reactingsingle.dt = 0
                self.addSingleEvent( reactingsingle )

            return
        


        #
        # 0. Reaction
        #
        if pair.eventType == EventType.REACTION:

            if __debug__:
                log.info( 'reaction' )

            if len( pair.rt.products ) == 1:
                
                species3 = pair.rt.products[0]

                rnd = uniform( size=2 )

                # calculate new R
            
                sgf = FirstPassageGreensFunction(pair.D_R)
                sgf.seta(pair.a_R)

                r_R = pair.drawR_single( sgf, pair.dt )
            
                displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R
                
                if __debug__:
                    shellSize = self.getShell(pair.shellID)[1]
                    assert self.distance( oldCoM, newCoM ) + species3.radius <\
                        shellSize

                #FIXME: SURFACE
                self.applyBoundary( newCoM )

                self.removeParticle( particle1 )
                self.removeParticle( particle2 )

                particle = self.createParticle( species3, newCoM )
                newsingle = self.createSingle( particle )
                self.addSingleEvent( newsingle )

                self.reactionEvents += 1

                self.lastReaction = Reaction( pair.rt, [particle1, particle2],
                                              [particle] )

                if __debug__:
                    log.info( 'product; %s' % str( newsingle ) )

            else:
                raise NotImplementedError,\
                      'num products >= 2 not supported.'

            self.removeFromShellMatrix( pair )

            return


        #
        # Escape 
        #

        r0 = self.distance( particle1.pos, particle2.pos )

        # 1 Escaping through a_r.
        if pair.eventType == EventType.ESCAPE:

            rnd = uniform( size=4 )

            # calculate new R
            
            sgf = FirstPassageGreensFunction(pair.D_R)
            sgf.seta(pair.a_R)

            r_R = pair.drawR_single( sgf, pair.dt )
                
            displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R

            # calculate new r
            theta_r = pair.drawTheta_pair( rnd[2], pair.a_r, r0, pair.dt, 
                                           pair.a_r )
            phi_r = rnd[3] * 2 * Pi
            newInterParticleS = numpy.array( [ pair.a_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
                
            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
            self.applyBoundary( newpos1 )
            self.applyBoundary( newpos2 )


        # 2 escaping through a_R.
        elif pair.eventType == 2:

            rnd = uniform( size = 4 )

            # calculate new r
            r = pair.drawR_pair( r0, pair.dt, pair.a_r )
            if __debug__:
                log.debug( 'new r = %g' % r )
            #assert r >= pair.sigma
            
            theta_r = pair.drawTheta_pair( rnd[0], r, r0, pair.dt, pair.a_r )
            phi_r = rnd[1] * 2*Pi
            newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
                
            # calculate new R
            displacement_R_S = [ pair.a_R, rnd[2] * Pi, rnd[3] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            
            newCoM = oldCoM + displacement_R
                
            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
            self.applyBoundary( newpos1 )
            self.applyBoundary( newpos2 )

        else:
            raise SystemError, 'Bug: invalid eventType.'

        # this has to be done before the following clearVolume()

        assert self.checkPairPos(pair, newpos1, newpos2, oldCoM,\
                                     self.getShell(pair.shellID)[1] )

        self.removeFromShellMatrix( pair )

        assert self.checkOverlap( newpos1, particle1.species.radius,
                                  ignore = [ particle1, particle2 ] )
        assert self.checkOverlap( newpos2, particle2.species.radius,
                                  ignore = [ particle1, particle2 ] )

        single1, single2 = pair.single1, pair.single2

        single1.initialize(self.t)
        single2.initialize(self.t)

        self.moveSingle(single1, newpos1, single1.particle.species.radius)
        self.moveSingle(single2, newpos2, single2.particle.species.radius)
            
        self.addSingleEvent( single1 )
        self.addSingleEvent( single2 )

        assert self.checkObj( single1 )
        assert self.checkObj( single2 )

        return

    def fireMulti( self, multi ):
        sim = multi.sim

        sim.step()
        #sim.sync()

        if sim.lastReaction:
            if __debug__:
                log.info( 'bd reaction' )

            self.breakUpMulti( multi )
            self.reactionEvents += 1
            self.lastReaction = sim.lastReaction
            return

        if sim.escaped:
            if __debug__:
                log.info( 'multi particle escaped.' )

            self.breakUpMulti( multi )
            return

        #if __debug__:
        #   log.info( 'multi stepped %d steps, duration %g, dt = %g' %
        #          ( additionalSteps + 1, sim.t - startT + sim.dt, dt ) )

        self.addMultiEvent(multi)
        return


    def breakUpMulti( self, multi ):
        self.removeFromShellMatrix( multi )

        singles = []
        for particle in multi.sim.particleList:
            single = self.createSingle( particle )
            self.addSingleEvent( single )
            singles.append( single )

        return singles

    def burstMulti( self, multi ):
        #multi.sim.sync()
        singles = self.breakUpMulti( multi )

        return singles

    def burstSingle( self, single ):
        assert self.t >= single.lastTime
        assert self.t <= single.lastTime + single.dt

        dt = self.t - single.lastTime

        particleRadius = single.particle.species.radius

        oldpos, shellSize = self.getShell(single.shellID)
        r = single.drawR(dt, shellSize)
        displacement = randomVector( r )

        newpos = oldpos + displacement

        self.applyBoundary( newpos )

        assert self.distance( newpos, oldpos ) <=\
            self.calculateSingleMobilityRadius(single)
        assert self.distance( newpos, oldpos ) - r <= r * 1e-6
        assert self.checkOverlap( newpos, particleRadius,\
                                  ignore = [ single.particle, ] )

        single.initialize( self.t )

        self.moveSingle( single, newpos, single.particle.species.radius )

        self.updateEvent( self.t, single )

    def breakUpPair( self, pair ):
        assert self.t >= pair.lastTime
        assert self.t <= pair.lastTime + pair.dt

        dt = self.t - pair.lastTime 

        if dt > 0.0:

            single1 = pair.single1
            single2 = pair.single2
            particle1 = single1.particle
            particle2 = single2.particle

            pos1 = particle1.pos
            pos2 = particle2.pos

            oldInterParticle = pos2 - pos1
            oldCoM = self.calculatePairCoM(pair)
            r0 = self.distance(pos1, pos2)
            
            rnd = uniform( size = 4 )

            sgf = FirstPassageGreensFunction(pair.D_R)
            sgf.seta(pair.a_R)

            # calculate new CoM
            r_R = pair.drawR_single( sgf, dt )
            
            displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R
            
            # calculate new interparticle
            r_r = pair.drawR_pair( r0, dt, pair.a_r )
            theta_r = pair.drawTheta_pair( rnd[2], r_r, r0, dt, pair.a_r )
            phi_r = rnd[3] * 2 * Pi
            newInterParticleS = numpy.array( [ r_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
            self.applyBoundary( newpos1 )
            self.applyBoundary( newpos2 )
            assert self.checkOverlap( newpos1, particle1.species.radius,
                                      ignore = [ particle1, particle2 ] )
                                      
            assert self.checkOverlap( newpos2, particle2.species.radius,
                                      ignore = [ particle1, particle2 ] )
                                      
            if __debug__:
                shellSize = self.getShell(pair.shellID)[1]
                assert self.checkPairPos(pair, newpos1, newpos2, oldCoM,\
                                             shellSize)


        return pair.single1, pair.single2

    def burstPair( self, pair ):
        single1, single2 = self.breakUpPair( pair )
        single1.initialize( self.t )
        single2.initialize( self.t )
        
        self.removeFromShellMatrix( pair )
        self.moveSingle(single1, single1.particle.pos, single1.particle.species.radius)
        self.moveSingle(single2, single2.particle.pos, single2.particle.species.radius)

        return single1, single2

    def formPairOrMulti( self, single, singlepos, neighbors ):
        assert neighbors

        # Try forming a Pair.
        if isinstance( neighbors[0], Single ):
            obj = self.formPair( single, singlepos, neighbors[0], neighbors[1:] )
            if obj:
                return obj


        # Then, a Multi.
        minShell = single.particle.species.radius * ( 1.0 + self.MULTI_SHELL_FACTOR )
        neighborDists = self.objDistanceArray( singlepos, neighbors )
        neighbors = [ neighbors[i] for i in 
                      ( neighborDists <= minShell ).nonzero()[0] ]

        if not neighbors:
            return None

        closest = neighbors[0]

        if isinstance( closest, Single ):

            multi = self.createMulti()
            self.addToMulti( single, multi )
            self.removeFromShellMatrix( single )
            for neighbor in neighbors:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )
            
            #self.addToShellMatrix( multi )
            self.addMultiEvent( multi )

            return multi

        elif isinstance( closest, Multi ):

            multi = closest
            if __debug__:
                log.info( 'multi merge %s %s' % ( single, multi ) )

            self.removeFromShellMatrix( multi )

            self.addToMulti( single, multi )
            self.removeFromShellMatrix( single )
            for neighbor in neighbors[1:]:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )

            #self.addToShellMatrix( multi )
            self.updateEvent( self.t + multi.dt, multi )

            return multi


        assert False, 'do not reach here'

    def formPair( self, single1, pos1, single2, bursted ):
        if __debug__:
           log.debug( 'trying to form %s' %
                  'Pair( %s, %s )' % ( single1.particle, 
                                       single2.particle ) )

        assert single1.isReset()
        assert single2.isReset()

        species1 = single1.particle.species
        species2 = single2.particle.species

        radius1 = species1.radius
        radius2 = species2.radius
        sigma = radius1 + radius2

        D1, D2 = species1.D, species2.D
        D12 = D1 + D2

        assert (pos1 - self.getShell(single1.shellID)[0]).sum() == 0
        pos2 = self.getShell(single2.shellID)[0]
        pairDistance = self.distance(pos1, pos2)
        r0 = pairDistance - sigma
        assert r0 >= 0, 'r0 (pair gap) between %s and %s = %g < 0' \
            % ( single1, single2, r0 )

        shellSize1 = pairDistance * D1 / D12 + radius1
        shellSize2 = pairDistance * D2 / D12 + radius2
        shellSizeMargin1 = radius1 * 2 #* self.SINGLE_SHELL_FACTOR
        shellSizeMargin2 = radius2 * 2 #* self.SINGLE_SHELL_FACTOR
        shellSizeWithMargin1 = shellSize1 + shellSizeMargin1
        shellSizeWithMargin2 = shellSize2 + shellSizeMargin2
        if shellSizeWithMargin1  >= shellSizeWithMargin2:
            minShellSize = shellSize1
            shellSizeMargin = shellSizeMargin1
        else:
            minShellSize = shellSize2
            shellSizeMargin = shellSizeMargin2

        # 1. Shell cannot be larger than max shell size or sim cell size.
        com = calculate_pair_CoM( pos1, pos2, D1, D2, self.getWorldSize() )
        self.applyBoundary( com )
        minShellSizeWithMargin = minShellSize + shellSizeMargin
        maxShellSize = min( self.getMaxShellSize(),
                            r0 * 100 + sigma + shellSizeMargin )

        if minShellSizeWithMargin >= maxShellSize:
            if __debug__:
                log.debug( '%s not formed: minShellSize >= maxShellSize' %
                       ( 'Pair( %s, %s )' % ( single1.particle, 
                                              single2.particle ) ) )
            return None

        # Here, we have to take into account of the burst Singles in
        # this step.  The simple check for closest below could miss
        # some of them, because sizes of these Singles for this
        # distance check has to include SINGLE_SHELL_FACTOR, while
        # these burst objects have zero mobility radii.  This is not
        # beautiful, a cleaner framework may be possible.

        closest, closestShellDistance = DummySingle(), numpy.inf
        for b in bursted:
            if isinstance( b, Single ):
                bpos = self.getShell(b.shellID)[0]
                d = self.distance( com, bpos ) \
                    - b.particle.species.radius * ( 1.0 + self.SINGLE_SHELL_FACTOR )
                if d < closestShellDistance:
                    closest, closestShellDistance = b, d

        if closestShellDistance <= minShellSizeWithMargin:
            if __debug__:
                log.debug( '%s not formed: squeezed by burst neighbor %s' %
                       ( 'Pair( %s, %s )' % ( single1.particle, 
                                              single2.particle ), closest ) )
            return None

        assert closestShellDistance > 0
        c, d = self.getClosestObj( com, ignore=[ single1, single2 ] )
        if d < closestShellDistance:
            closest, closestShellDistance = c, d

        if __debug__:
            log.debug( 'Pair closest neighbor: %s %g, minShellWithMargin %g' %
                   ( closest, closestShellDistance, minShellSizeWithMargin ) )

        assert closestShellDistance > 0

        if isinstance( closest, Single ):

            D_closest = closest.particle.species.D
            D_tot = D_closest + D12
            closestDistance = self.distance( com, closest.particle.pos ) ##??

            closestMinRadius = closest.particle.species.radius
            closestMinShell = closestMinRadius * \
                ( self.SINGLE_SHELL_FACTOR + 1.0 )

            shellSize = min( ( D12 / D_tot ) *
                             ( closestDistance - minShellSize 
                               - closestMinRadius ) + minShellSize,
                             closestDistance - closestMinShell,
                             closestShellDistance )

            shellSize /= SAFETY
            assert shellSize < closestShellDistance

        else:
            assert isinstance( closest, ( Pair, Multi, DummySingle ) )

            shellSize = closestShellDistance / SAFETY

        if shellSize <= minShellSizeWithMargin:
            if __debug__:
                log.debug( '%s not formed: squeezed by %s' %
                       ( 'Pair( %s, %s )' % ( single1.particle, 
                                              single2.particle ), closest ) )
            return None


        d1 = self.distance(com, pos1)
        d2 = self.distance(com, pos2)

        if shellSize < max( d1 + single1.particle.species.radius *
                            ( 1.0 + self.SINGLE_SHELL_FACTOR ), \
                                d2 + single2.particle.species.radius * \
                                ( 1.0 + self.SINGLE_SHELL_FACTOR ) ) * 1.3:
            if __debug__:
                log.debug( '%s not formed: singles are better' %
                       'Pair( %s, %s )' % ( single1.particle, 
                                            single2.particle ) )
            return None

        # 3. Ok, Pair makes sense.  Create one.
        shellSize = min( shellSize, maxShellSize )

        pair = self.createPair(single1, single2, shellSize)

        self.determinePairEvent(pair, self.t, pos1, pos2, shellSize)

        self.removeFromShellMatrix( single1 )
        self.removeFromShellMatrix( single2 )

        self.addPairEvent( pair )
        # single1 will be removed by the scheduler.
        self.removeEvent( single2 )

        assert closestShellDistance == numpy.inf or shellSize < closestShellDistance
        assert shellSize >= minShellSizeWithMargin
        assert shellSize <= maxShellSize

        if __debug__:
            log.info( '%s, dt=%g, pairDistance=%g, shell=%g,' %
                  ( pair, pair.dt, pairDistance, shellSize ) + 
                  'closest=%s, closestShellDistance=%g' %
                  ( closest, closestShellDistance ) )

        assert self.checkObj( pair )

        return pair
    
    def addToMultiRecursive( self, obj, multi ):
        if isinstance( obj, Single ):
            if obj.particle in multi.sim.particleList:  # Already in the Multi.
                return
            assert obj.isReset()
            objpos = self.getShell(obj.shellID)[0]
            
            self.addToMulti( obj, multi )
            self.removeFromShellMatrix( obj )
            self.removeEvent( obj )

            radius = obj.particle.species.radius *\
                ( 1.0 + self.MULTI_SHELL_FACTOR )
            neighbors = self.getNeighborsWithinRadiusNoSort( objpos, radius,
                                                             ignore=[obj,] )
            bursted = self.burstNonMultis( neighbors )
            neighborDists = self.objDistanceArray( objpos, bursted )
            neighbors = [ bursted[i] for i in 
                          ( neighborDists <= radius ).nonzero()[0] ]

            for obj in neighbors:
                self.addToMultiRecursive( obj, multi )

        elif isinstance( obj, Multi ):
            if not obj.sim.particleList[0] in multi.sim.particleList:
                self.mergeMultis( obj, multi )
                self.removeFromShellMatrix( obj )
                self.removeEvent( obj )
            else:
                if __debug__:
                    log.debug( '%s already added. skipping.' % obj )
        else:
            assert False, 'do not reach here.'  # Pairs are bursted

    def addToMulti( self, single, multi ):
        if __debug__:
            log.info( 'adding %s to %s' % ( single, multi ) )

        shellSize = single.particle.species.radius * \
            ( 1.0 + self.MULTI_SHELL_FACTOR )
        multi.addParticle( single.particle )
        multi.addShell( single.particle.pos, shellSize )

    def mergeMultis( self, multi1, multi2 ):
        '''
        merge multi1 into multi2
        '''
        if __debug__:
            log.info( 'merging %s to %s' % ( multi1, multi2 ) )

        assert not multi1.sim.particleList[0] in multi2.sim.particleList

        for i, particle in enumerate( multi1.sim.particleList ):
            
            # FIXME: shells should be renewed

            multi2.addParticle( particle )
            shellID = multi1.shellIDList[i]
            shell = self.getShell(shellID)
            multi2.addShell( shell[0], shell[1] )

        multi2.initialize( self.t )

    def getNeighborShells( self, pos, n=None ):
        '''
        Find closest n shells.

        This method returns a tuple ( neighbors, distances ).
        '''
        neighbors, distances = self.shellMatrix.getNeighbors( pos, n )

        if len( neighbors ) == 0:
            return [( DummySingle(), 0 ),], [numpy.inf,]
        return neighbors, distances

    def getNeighborShellsNoSort( self, pos, n=None ):
        return self.shellMatrix.getNeighborsNoSort( pos, n )

    def getNeighborShellsWithinRadius( self, pos, radius ):
        return self.shellMatrix.getNeighborsWithinRadius( pos, radius )

    def getNeighborShellsWithinRadiusNoSort( self, pos, radius ):
        return self.shellMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

    def getNeighborsWithinRadius( self, pos, radius, ignore=[] ):
        shells, distances =\
            self.shellMatrix.getNeighborsWithinRadius( pos, radius )

        neighbors = [ s[0] for s in shells if s[0] not in ignore ]
        neighbors = uniq( neighbors )

        return neighbors

    def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        shells, distances =\
            self.shellMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        neighbors = uniq( [ s[0] for s in shells if s[0] not in ignore ] )

        return neighbors

    def getNeighbors( self, pos, radius=numpy.inf, ignore=[] ):
        shells, dists = self.shellMatrix.getNeighbors( pos )

        seen = dict.fromkeys( ignore )
        neighbors = []
        distances = []

        for i, shell in enumerate( shells ):
            if not shell[0] in seen:
                seen[ shell ] = None
                neighbors.append(shell[0])
                distances.append(dists[i])
                if dists[i] > radius:
                    return neighbors, distances

        return neighbors + [DummySingle()], numpy.concatenate( [ distances,
                                                                 [numpy.inf] ] )

    def getClosestObj( self, pos, ignore=[] ):
        shells, distances = self.getNeighborShells( pos )

        for i, shell in enumerate( shells ):
            neighbor = shell[0]
            if neighbor not in ignore:
                return neighbor, distances[i]

        return DummySingle(), numpy.inf

    def objDistance( self, pos, obj ):
        dists = numpy.zeros( len( obj.shellIDList ) )
        for i, shellID in enumerate(obj.shellIDList):
            shellpos, shellradius = self.getShell(shellID)
            dists[i] = self.distance( pos, shellpos ) - shellradius
        return min( dists )

    def objDistanceArray( self, pos, objs ):
        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        return dists
            

    #
    # consistency checkers
    #

    def checkObj( self, obj ):
        obj.check()

        for shellID in obj.shellIDList:

            shellpos, radius = self.getShell(shellID)

            closest, distance = self.getClosestObj( shellpos,
                                                    ignore = [obj] )

            assert radius <= self.getUserMaxShellSize(),\
                '%s shell size larger than user-set max shell size' % \
                str( shellID )

            assert radius <= self.getMaxShellSize(),\
                '%s shell size larger than simulator cell size / 2' % \
                str( shellID )

            assert distance - radius >= 0.0,\
                '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                % ( str( obj ), str( closest ), radius, distance,\
                        distance - radius )

        return True

    def checkObjForAll( self ):
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            self.checkObj( obj )

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

        shellPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            shellPopulation += len( obj.shellIDList )

        if shellPopulation != self.shellMatrix.size:
            raise RuntimeError,\
                'num shells != self.shellMatrix.size'
        
        self.shellMatrix.check()


    def checkPairPos( self, pair, pos1, pos2, com, radius ):
        species1 = pair.single1.particle.species
        species2 = pair.single2.particle.species

        oldCoM = com
        
        # debug: check if the new positions are valid:
        newDistance = distance_Simple( pos1, pos2 )
        particleRadius12 = species1.radius + species2.radius

        # check 1: particles don't overlap.
        if newDistance <= particleRadius12:
            if __debug__:
                log.info( 'rejected move: radii %g, particle distance %g',
                          ( species1.radius + species2.radius, newDistance ) )
            if __debug__:
                log.debug( 'DEBUG: pair.dt %g, pos1 %s, pos2 %s' %
                           ( pair.dt, str( pos1 ), str( pos2 ) ) )
            raise RuntimeError, 'New particles overlap'

        # check 2: particles within mobility radius.
        d1 = self.distance( oldCoM, pos1 ) + species1.radius
        d2 = self.distance( oldCoM, pos2 ) + species2.radius
        if d1 > radius or d2 > radius:
            raise RuntimeError, \
                'New particle(s) out of protective sphere. %s' % \
                'radius = %g, d1 = %g, d2 = %g ' % ( radius, d1, d2 )
                
        

        return True




    def check( self ):
        ParticleSimulatorBase.check( self )

        assert self.scheduler.check()

        assert self.t >= 0.0
        assert self.dt >= 0.0

        self.checkShellMatrix()

        self.checkEventStoichiometry()
        
        self.checkObjForAll()

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



