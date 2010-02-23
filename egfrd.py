#!/usr/env python


from weakref import ref
import math

import numpy

from _gfrd import (
    EventScheduler,
    FirstPassageGreensFunction,
    FirstPassagePairGreensFunction,
    FirstPassageNoCollisionPairGreensFunction,
    BasicPairGreensFunction,
    FreePairGreensFunction,
    EventType,
    Particle,
    SphericalShell,
    SphericalShellContainer,
    DomainIDGenerator,
    ShellIDGenerator,
    DomainID,
    ParticleContainer
    )

from surface import CuboidalRegion

from gfrdbase import *
from single import *
from pair import *
from multi import *
from utils import *
import myrandom

import logging
import os

log = logging.getLogger('ecell')

class Delegate( object ):
    def __init__( self, obj, method ):
        self.ref = ref( obj )
        self.method = method

    def __call__( self, *arg ):
        return self.method( self.ref(), *arg )


class EGFRDSimulator( ParticleSimulatorBase ):
    def __init__( self ):
        self.shellMatrix = None
        self.domainIDGenerator = DomainIDGenerator(0)
        self.shellIDGenerator = ShellIDGenerator(0)

        ParticleSimulatorBase.__init__( self )

        self.MULTI_SHELL_FACTOR = 0.05
        self.SINGLE_SHELL_FACTOR = 0.1

        self.isDirty = True
        self.scheduler = EventScheduler()

        self.smallT = 1e-8  # FIXME: is this ok?

        self.userMaxShellSize = numpy.inf

        self.domains = {}

        self.reset()

    def setWorldSize( self, size ):
        ParticleSimulatorBase.setWorldSize( self, size )
        self.shellMatrix = SphericalShellContainer( self.worldSize, self.matrixSize )

    def setMatrixSize( self, size ):
        ParticleSimulatorBase.setMatrixSize( self, size )
        self.shellMatrix = SphericalShellContainer( self.worldSize, self.matrixSize )

    def getMatrixCellSize( self ):
        return self.shellMatrix.cell_size

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
        self.single_steps = {EventType.SINGLE_ESCAPE:0,
                             EventType.SINGLE_REACTION:0}
        self.pair_steps = {EventType.SINGLE_REACTION:0,
                           EventType.PAIR_REACTION:0,
                           EventType.IV_ESCAPE:0,
                           EventType.COM_ESCAPE:0}
        self.multi_steps = {EventType.MULTI_ESCAPE:0,
                            EventType.MULTI_REACTION:0, 2:0}
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
        self.shellMatrix = SphericalShellContainer(self.worldSize, self.matrixSize)
        self.domains = {}

        singles = []
        for pid_particle_pair in self.particleMatrix:
            singles.append(self.createSingle(pid_particle_pair))
        assert len(singles) == len(self.particleMatrix)
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
            
        if __debug__:
            if int("0" + os.environ.get("ECELL_CHECK", ""), 10):
                self.check()
        
        self.stepCounter += 1

        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getArg()

        if __debug__:
            domain_counts = self.count_domains()
            log.info( '\n%d: t=%g dt=%g\tSingles: %d, Pairs: %d, Multis: %d'
                      % (( self.stepCounter, self.t, self.dt ) + domain_counts ))
            log.info( 'event=%s reactions=%d rejectedmoves=%d' 
                      % ( self.lastEvent, self.reactionEvents, 
                          self.rejectedMoves ) )
        
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

    def createSingle( self, pid_particle_pair ):
        rt = self.getReactionRule1(pid_particle_pair[1].sid)
        domain_id = self.domainIDGenerator()
        shell_id = self.shellIDGenerator()

        # Get surface.
        species = self.speciesList[pid_particle_pair[1].sid]
        surface = self.getSurface(species)

        # Create single. The type of the single that will be created depends 
        # on the surface this particle is on. Either SphericalSingle, 
        # PlanarSurfaceSingle, or CylindricalSurfaceSingle.
        TypeOfSingle = surface.DefaultSingle
        single = TypeOfSingle(domain_id, pid_particle_pair, shell_id, rt)

        single.initialize(self.t)
        self.moveShell(single.shell)
        self.domains[domain_id] = single
        return single

    def createPair(self, single1, single2, shellSize):
        assert single1.dt == 0.0
        assert single2.dt == 0.0

        rt = self.getReactionRule2(single1.pid_particle_pair[1].sid, single2.pid_particle_pair[1].sid)[ 0 ]

        domain_id = self.domainIDGenerator()
        shell_id = self.shellIDGenerator()

        CoM = calculate_pair_CoM(single1.pid_particle_pair[1].position,
                                 single2.pid_particle_pair[1].position,
                                 single1.pid_particle_pair[1].D,
                                 single2.pid_particle_pair[1].D,
                                 self.worldSize)

        pos1 = single1.shell[1].position
        pos2 = single2.shell[1].position
        r0 = self.distance(pos1, pos2)

        # Get surface.
        species = self.speciesList[single1.pid_particle_pair[1].sid]
        surface = self.getSurface(species)

        # Create pair. The type of the pair that will be created depends on 
        # the surface the particles are on. Either SphericalPair, 
        # PlanarSurfacePair, or CylindricalSurfacePair.
        TypeOfPair = surface.DefaultPair
        pair = TypeOfPair(domain_id, CoM, single1, single2, shell_id, 
                          r0, shellSize, rt)

        pair.initialize( self.t )

        self.moveShell(pair.shell)
        self.domains[domain_id] = pair
        return pair

    def createMulti( self ):
        domain_id = self.domainIDGenerator()
        multi = Multi( domain_id, self )
        self.domains[domain_id] = multi
        return multi

    def moveSingle(self, single, position, radius=None):
        self.moveSingleShell(single, position, radius)
        self.moveSingleParticle(single, position)

    def moveSingleShell(self, single, position, radius):
        if radius is None:
            # Shrink single to particle radius by default.
            radius = single.shell[1].radius

        # Reuse shell_id and domain_id.
        shell_id = single.shell[0]
        domain_id = single.domain_id

        # Replace shell.
        shell = single.createNewShell(position, radius, domain_id)
        new_sid_shell_pair = (shell_id, shell) 

        single.shell = new_sid_shell_pair
        self.moveShell(new_sid_shell_pair)

        # Rescale domains of single.
        single.rescaleCoordinates(radius)

    def moveSingleParticle(self, single, position):
        new_pid_particle_pair = (single.pid_particle_pair[0],
                          Particle(position,
                                   single.pid_particle_pair[1].radius,
                                   single.pid_particle_pair[1].D,
                                   single.pid_particle_pair[1].sid))
        single.pid_particle_pair = new_pid_particle_pair

        self.moveParticle(new_pid_particle_pair)

    def removeDomain( self, obj ):
        del self.domains[obj.domain_id]
        for shell_id, _ in obj.shell_list:
            del self.shellMatrix[shell_id]
    
    def moveShell(self, shell_id_shell_pair):
        self.shellMatrix.update(shell_id_shell_pair)

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
            # TODO. Compare with gfrd.
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
        reactantSpeciesRadius = single.pid_particle_pair[1].radius
        oldpos = single.pid_particle_pair[1].position
        
        rt = single.drawReactionRule()

        if len( rt.products ) == 0:
            
            self.removeParticle(single.pid_particle_pair)

            self.lastReaction = Reaction( rt, [single.pid_particle_pair[1]], [] )

            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            if reactantSpeciesRadius < productSpecies.radius:
                self.clearVolume( oldpos, productSpecies.radius )

            if self.checkOverlap( oldpos, productSpecies.radius,
                                  ignore = [ single.pid_particle_pair[0], ] ):
                if __debug__:
                    log.info( 'no space for product particle.' )
                raise NoSpace()

            self.removeParticle(single.pid_particle_pair)
            newparticle = self.createParticle( productSpecies.id, oldpos )
            newsingle = self.createSingle( newparticle )
            self.addSingleEvent( newsingle )

            self.lastReaction = Reaction( rt, [single.pid_particle_pair[1]], [newparticle] )

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

            for _ in range(self.dissociation_retry_moves):
                unitVector = randomUnitVector()
                vector = unitVector * particleRadius12 * ( 1.0 + 1e-7 )
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?

                while 1:
                    newpos1 = oldpos + vector * ( D1 / D12 )
                    newpos2 = oldpos - vector * ( D2 / D12 )
                    newpos1 = self.applyBoundary( newpos1 )
                    newpos2 = self.applyBoundary( newpos2 )

                    if self.distance( newpos1, newpos2 ) >= particleRadius12:
                        break

                    vector *= 1.0 + 1e-7


                # accept the new positions if there is enough space.
                if (not self.checkOverlap(newpos1, particleRadius1,
                                          ignore=[single.pid_particle_pair[0]])) and \
                   (not self.checkOverlap(newpos2, particleRadius2,
                                          ignore=[single.pid_particle_pair[0]])):
                    break
            else:
                if __debug__:
                    log.info( 'no space for product particles.' )
                raise NoSpace()

            self.removeParticle(single.pid_particle_pair)

            particle1 = self.createParticle(productSpecies1.id, newpos1)
            particle2 = self.createParticle(productSpecies2.id, newpos2)
            newsingle1 = self.createSingle(particle1)
            newsingle2 = self.createSingle(particle2)

            self.addSingleEvent( newsingle1 )
            self.addSingleEvent( newsingle2 )

            self.lastReaction = Reaction( rt, [single.pid_particle_pair[1]], 
                                          [particle1, particle2] )

            if __debug__:
                log.info( 'products; %s %s' % 
                      ( str( newsingle1 ), str( newsingle2 ) ) )

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1

    def propagateSingle(self, single, isEscape=False, isBurst=False):
        """The difference between a burst and a propagate is that a burst 
        always takes place before the actual scheduled event for the single, 
        while propagateSingle can be called for an escape event.

        Another subtle difference is that burstSingle always reschedules 
        (updateEvent) the single, while just calling propagate does not. 
        So whoever calls propagateSingle directly should reschedule the single 
        afterwards.

        """
        if __debug__:
            log.debug( "single.dt=%g, single.lastTime=%g, self.t=%g" % (
                single.dt, single.lastTime, self.t ) )

        if not isBurst:
            assert abs(single.dt + single.lastTime - self.t) <= 1e-18 * self.t
        
        newpos = single.drawNewPosition(self.t - single.lastTime, isEscape) 
        newpos = self.applyBoundary(newpos)

        if __debug__:
            log.debug( "propagate %s: %s => %s" % ( single, single.pid_particle_pair[1].position, newpos ) )

        if self.checkOverlap(newpos,
                             single.pid_particle_pair[1].radius,
                             ignore=[single.pid_particle_pair[0]]):
            raise RuntimeError('propagateSingle: checkOverlap failed.')

        if single.eventType == EventType.SINGLE_REACTION and isBurst == False:
            # SINGLE_REACTION, and not a burst. No need to update, single is 
            # removed anyway.
            pass
        else:
            # Todo. if isinstance(single, InteractionSingle):
            single.initialize(self.t)
            self.moveSingle(single, newpos, single.pid_particle_pair[1].radius)

    def fireSingle( self, single ):
        # Reaction.
        if single.eventType == EventType.SINGLE_REACTION:
            if __debug__:
                log.info('fireSingle: eventType %s' % single.eventType)

            self.single_steps[single.eventType] += 1

            if __debug__:
                log.info( 'single reaction %s' % str( single ) )

            self.propagateSingle(single, isEscape=False)

            try:
                self.removeDomain( single )
                self.fireSingleReaction( single )
                return
            except NoSpace:
                if __debug__:
                    log.info( 'single reaction; placing product failed.' )
                self.shellMatrix.update( single.shell )
                self.rejectedMoves += 1
                single.reset()
                self.addSingleEvent(single)
                return

        # Propagate, if not reaction.
        single.eventType = EventType.SINGLE_ESCAPE
        if __debug__:
            log.info('fireSingle: eventType %s' % single.eventType)
        self.single_steps[single.eventType] += 1

        # Handle immobile case first.
        if single.getD() == 0:
            # no propagation, just calculate next reaction time.
            single.dt, single.eventType, single.activeCoordinate = \
                single.determineNextEvent() 
            single.lastTime = self.t
            self.addSingleEvent(single)
            return
        
        if single.dt != 0.0:
            # Propagate this particle to the exit point on the shell.
            self.propagateSingle(single, isEscape=True)

        singlepos = single.shell[1].position

        # (2) Clear volume.

        minShell = single.pid_particle_pair[1].radius * ( 1.0 + self.SINGLE_SHELL_FACTOR )

        intruders = []   # intruders are domains within minShell
        closest = None   # closest is the closest domain, excluding intruders.
        closestDistance = numpy.inf # distance to the shell of the closet.
        
        neighbors = self.getNeighbors(singlepos)
        seen = set([single.domain_id,])
        for n in neighbors:
            did = n[0][1].did
            distance = n[1]
            if distance > minShell:
                closest = self.domains[did]
                closestDistance = distance
                break
            elif did not in seen:
                seen.add(did)
                intruders.append(self.domains[did])

        if __debug__:
            log.debug( "intruders: %s, closest: %s (dist=%g)" %\
                           (intruders, closest, closestDistance) )

        burst = []
        if intruders:
            burst = self.burstNonMultis(intruders)

            obj = self.formPairOrMulti(single, singlepos, burst)

            if obj:
                return

            # if nothing was formed, recheck closest and restore shells.
            closest, closestDistance = \
                self.getClosestObj( singlepos, ignore = [ single.domain_id, ] )

        self.updateSingle( single, closest, closestDistance )

        burst = uniq( burst )
        burstSingles = [ s for s in burst if isinstance( s, Single ) ]
        self.restoreSingleShells( burstSingles )
            
        if __debug__:
            log.info( 'single shell %s dt %g.' %\
                          ( single.shell, single.dt ) )

        self.addSingleEvent(single)
        return

    def restoreSingleShells( self, singles ):
        for single in singles:
            assert single.isReset()
            c, d = self.getClosestObj( single.shell[1].position, ignore = [single.domain_id,] )

            self.updateSingle( single, c, d )
            self.updateEvent( self.t + single.dt, single )
            if __debug__:
                log.debug( 'restore shell %s %g dt %g closest %s %g' %
                       ( single, single.shell[1].radius, single.dt, c, d ) )

    def calculateSingleShellSize( self, single, closest, 
                                  distance, shellDistance ):
        assert isinstance( closest, Single )

        minRadius1 = single.pid_particle_pair[1].radius
        D1 = single.getD()

        if D1 == 0:
            return minRadius1

        D2 = closest.getD()
        minRadius2 = closest.pid_particle_pair[1].radius
        minRadius12 = minRadius1 + minRadius2
        sqrtD1 = math.sqrt( D1 )
            
        shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                         * ( distance - minRadius12 ) + minRadius1,
                         shellDistance / SAFETY )
        if shellSize < minRadius1:
            shellSize = minRadius1

        return shellSize

    def updateSingle( self, single, closest, distanceToShell ): 
        # Todo. assert not isinstance(single, InteractionSingle)

        singlepos = single.shell[1].position
        if isinstance( closest, Single ):
            closestpos = closest.shell[1].position
            distanceToClosest = self.distance(singlepos, closestpos)
            new_shell_size = self.calculateSingleShellSize( single, closest, 
                                                       distanceToClosest,
                                                       distanceToShell )
        else:  # Pair or Multi
            new_shell_size = distanceToShell / SAFETY
            new_shell_size = max( new_shell_size, single.pid_particle_pair[1].radius )

        new_shell_size = min( new_shell_size, self.getMaxShellSize() )

        # Resize shell, don't change position.
        # Note: this should be done before determineNextEvent.
        self.moveSingleShell(single, singlepos, new_shell_size)        

        single.dt, single.eventType, single.activeCoordinate = \
            single.determineNextEvent()
        single.lastTime = self.t

    def firePair( self, pair ):
        assert self.checkObj( pair )

        particle1 = pair.single1.pid_particle_pair
        particle2 = pair.single2.pid_particle_pair
        
        oldInterParticle = particle2[1].position - particle1[1].position
        oldCoM = calculate_pair_CoM(
            particle1[1].position,
            particle2[1].position,
            particle1[1].D,
            particle2[1].D,
            self.worldSize)
        oldCoM = self.applyBoundary( oldCoM )

        # Two cases:
        #  1. Single reaction
        #  2. Not a single reaction

        #
        # 1. Single reaction
        #
        if pair.eventType == EventType.SINGLE_REACTION:
            self.pair_steps[pair.eventType] += 1
            if __debug__:
                log.info('firePair: eventType %s' % pair.eventType)

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
                self.removeDomain( reactingsingle )
                self.fireSingleReaction( reactingsingle )
            except NoSpace:
                self.shellMatrix.update( reactingsingle.shell )
                self.rejectedMoves += 1
                reactingsingle.dt = 0
                self.addSingleEvent( reactingsingle )

            return
        
        #
        # 2. Not a single reaction
        #

        # When this pair was initialized, pair.determineNextEvent() was called 
        # and the pair.activeCoordinate we use here was set.
        #
        # Possibilities if the IV coordinate is active:
        #   2a. Pair reaction
        #   2b. IV escape
        #
        # Possiblities if the CoM coordinate is active:
        #   2c. CoM escape
        pair.eventType = pair.activeCoordinate.drawEventType(pair.dt)
        self.pair_steps[pair.eventType] += 1
        if __debug__:
            log.info('firePair: eventType %s' % pair.eventType)

        #
        # 2a. Pair reaction
        #
        if pair.eventType == EventType.PAIR_REACTION:
            if __debug__:
                log.info( 'reaction' )

            if len( pair.rt.products ) == 1:
                
                species3 = pair.rt.products[0]

                # calculate new R
            
                sgf = FirstPassageGreensFunction(pair.D_R)
                sgf.seta(pair.a_R)

                r_R = pair.drawR_single( sgf, pair.dt )
            
                displacement_R_S = [ r_R, myrandom.uniform() * Pi,
                                     myrandom.uniform() * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R
                
                if __debug__:
                    shellSize = pair.shell[1].radius
                    assert self.distance( oldCoM, newCoM ) + species3.radius <\
                        shellSize

                #FIXME: SURFACE
                newCoM = self.applyBoundary( newCoM )

                self.removeParticle(pair.single1.pid_particle_pair)
                self.removeParticle(pair.single2.pid_particle_pair)

                particle = self.createParticle( species3.id, newCoM )
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

            self.removeDomain( pair )

            return


        r0 = self.distance( particle1[1].position, particle2[1].position )

        #
        # 2b. Escaping through a_r.
        #
        if pair.eventType == EventType.IV_ESCAPE:

            # calculate new R
            
            sgf = FirstPassageGreensFunction(pair.D_R)
            sgf.seta(pair.a_R)

            r_R = pair.drawR_single( sgf, pair.dt )
                
            displacement_R_S = [r_R, myrandom.uniform() * Pi, 
                                myrandom.uniform() * 2 * Pi]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R

            # calculate new r
            theta_r = pair.drawTheta_pair(myrandom.uniform(), pair.a_r, 
                                          r0, pair.dt, pair.a_r )
            phi_r = myrandom.uniform() * 2 * Pi
            newInterParticleS = numpy.array( [ pair.a_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
                
            newpos1, newpos2 = pair.calculatePairPos(newCoM,
                                                     newInterParticle,
                                                     oldInterParticle)
            newpos1 = self.applyBoundary( newpos1 )
            newpos2 = self.applyBoundary( newpos2 )
        
        #
        # 2c. escaping through a_R.
        #
        elif pair.eventType == EventType.COM_ESCAPE:

            # calculate new r
            r = pair.drawR_pair( r0, pair.dt, pair.a_r )
            if __debug__:
                log.debug( 'new r = %g' % r )
            #assert r >= pair.sigma
            
            theta_r = pair.drawTheta_pair(myrandom.uniform(), r, r0, pair.dt, 
                                          pair.a_r)
            phi_r = myrandom.uniform() * 2*Pi
            newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
                
            # calculate new R
            displacement_R_S = [ pair.a_R, myrandom.uniform() * Pi, 
                                 myrandom.uniform() * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            
            newCoM = oldCoM + displacement_R
                
            newpos1, newpos2 = pair.calculatePairPos(newCoM, 
                                                     newInterParticle,
                                                     oldInterParticle)
            newpos1 = self.applyBoundary( newpos1 )
            newpos2 = self.applyBoundary( newpos2 )

        else:
            raise SystemError, 'Bug: invalid eventType.'

        # this has to be done before the following clearVolume()

        assert self.checkPairPos(pair, newpos1, newpos2, oldCoM, pair.shell[1].radius)

        self.removeDomain( pair )

        assert not self.checkOverlap(newpos1, particle1[1].radius,
                                     ignore=[particle1[0], particle2[0]])
        assert not self.checkOverlap(newpos2, particle2[1].radius,
                                     ignore=[particle1[0], particle2[0]])

        single1, single2 = pair.single1, pair.single2

        single1.initialize(self.t)
        single2.initialize(self.t)

        if __debug__:
            log.debug("firePair: #1 { %s: %s => %s }" % (single1, particle1[1].position, newpos1))
            log.debug("firePair: #2 { %s: %s => %s }" % (single2, particle2[1].position, newpos2))

        self.domains[single1.domain_id] = single1
        self.domains[single2.domain_id] = single2

        self.moveSingle(single1, newpos1, single1.pid_particle_pair[1].radius)
        self.moveSingle(single2, newpos2, single2.pid_particle_pair[1].radius)
            
        self.addSingleEvent( single1 )
        self.addSingleEvent( single2 )

        assert self.checkObj( single1 )
        assert self.checkObj( single2 )

        return

    def fireMulti( self, multi ):
            
        sim = multi.sim

        self.multi_steps[2] += 1  # multi_steps[2]: total multi steps
        sim.step()

        if __debug__:
            event_type = ''
            if multi.sim.lastReaction:
                event_type = 'reaction'
            elif multi.sim.escaped:
                event_type = 'escaped'
            log.info('fireMulti: %s' % event_type)

        if sim.lastReaction:
            self.breakUpMulti( multi )
            self.reactionEvents += 1
            self.lastReaction = sim.lastReaction
            self.multi_steps[EventType.MULTI_REACTION] += 1
            return

        if sim.escaped:
            self.breakUpMulti( multi )
            self.multi_steps[EventType.MULTI_ESCAPE] += 1
            return

        self.addMultiEvent(multi)

        return

    def breakUpMulti( self, multi ):
        self.removeDomain( multi )

        singles = []
        for pid in multi.sim.particleList:
            single = self.createSingle((pid, multi.sim.particleMatrix[pid]))
            self.addSingleEvent( single )
            singles.append( single )

        return singles

    def burstMulti( self, multi ):
        #multi.sim.sync()
        assert isinstance(multi, Multi)
        singles = self.breakUpMulti( multi )

        return singles

    def burstSingle( self, single ):
        assert self.t >= single.lastTime
        assert self.t <= single.lastTime + single.dt

        oldpos, oldRadius, _ = single.shell[1]
        particleRadius = single.pid_particle_pair[1].radius

        self.propagateSingle(single, isEscape=False, isBurst=True)

        newpos = single.pid_particle_pair[1].position
        assert self.distance(newpos, oldpos) <= oldRadius - particleRadius
        # Displacement check is in NonInteractionSingle.drawNewPosition.

        # Todo. if isinstance(single, InteractionSingle):
        self.updateEvent( self.t, single )

        assert single.shell[1].radius == particleRadius

    def burstPair( self, pair ):
        if __debug__:
            log.debug('burstPair: %s', pair)

        assert self.t >= pair.lastTime
        assert self.t <= pair.lastTime + pair.dt

        single1 = pair.single1
        single2 = pair.single2

        dt = self.t - pair.lastTime 

        particle1 = single1.pid_particle_pair
        particle2 = single2.pid_particle_pair

        if dt > 0.0:

            pos1 = particle1[1].position
            pos2 = particle2[1].position
            D1 = particle1[1].D
            D2 = particle2[1].D

            oldInterParticle = pos2 - pos1
            oldCoM = calculate_pair_CoM(pos1, pos2, D1, D2, self.worldSize)
            r0 = self.distance(pos1, pos2)
            
            sgf = FirstPassageGreensFunction(pair.D_R)
            sgf.seta(pair.a_R)

            # calculate new CoM
            r_R = pair.drawR_single( sgf, dt )
            
            displacement_R_S = [r_R, myrandom.uniform() * Pi, 
                                 myrandom.uniform() * 2 * Pi]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R
            
            # calculate new interparticle
            r_r = pair.drawR_pair( r0, dt, pair.a_r )

            theta_r = pair.drawTheta_pair(myrandom.uniform(), r_r, r0, dt, 
                                          pair.a_r)
            phi_r = myrandom.uniform() * 2 * Pi
            newInterParticleS = numpy.array( [ r_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
            newpos1, newpos2 = pair.calculatePairPos(newCoM, 
                                                     newInterParticle,
                                                     oldInterParticle)

            newpos1 = self.applyBoundary(newpos1)
            newpos2 = self.applyBoundary(newpos2)
            assert not self.checkOverlap(newpos1, particle1[1].radius,
                                         ignore=[particle1[0], particle2[0]])
            assert not self.checkOverlap(newpos2, particle2[1].radius,
                                         ignore=[particle1[0], particle2[0]])
            assert self.checkPairPos(pair, newpos1, newpos2, oldCoM,\
                                         pair.shell[1].radius)
        else:
            newpos1 = particle1[1].position
            newpos2 = particle2[1].position

        single1.initialize( self.t )
        single2.initialize( self.t )
        
        self.removeDomain( pair )
        assert single1.domain_id not in self.domains
        assert single2.domain_id not in self.domains
        self.domains[single1.domain_id] = single1
        self.domains[single2.domain_id] = single2
        self.moveSingle(single1, newpos1, particle1[1].radius)
        self.moveSingle(single2, newpos2, particle2[1].radius)

        assert self.shellMatrix[single1.shell[0]].radius == single1.shell[1].radius
        assert self.shellMatrix[single2.shell[0]].radius == single2.shell[1].radius
        assert single1.shell[1].radius == particle1[1].radius
        assert single2.shell[1].radius == particle2[1].radius

        return single1, single2

    def formPairOrMulti( self, single, singlepos, neighbors ):
        assert neighbors

        # sort burst neighbors by distance
        dists = self.objDistanceArray(singlepos, neighbors)
        if len(dists) >= 2:
            n = dists.argsort()
            dists = dists.take(n)
            neighbors = numpy.take(neighbors, n)

        # First, try forming a Pair.
        if isinstance(neighbors[0], Single):
            obj = self.formPair(single, singlepos, neighbors[0], neighbors[1:])
            if obj:
                return obj

        # If a Pair is not formed, then try forming a Multi.
        obj = self.formMulti(single, neighbors, dists)
        if obj:
            return obj


    def formPair( self, single1, pos1, single2, burst ):
        if __debug__:
           log.debug( 'trying to form %s' %
                  'Pair( %s, %s )' % ( single1.pid_particle_pair, 
                                       single2.pid_particle_pair ) )

        assert single1.isReset()
        assert single2.isReset()

        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius

        sigma = radius1 + radius2

        D1, D2 = single1.pid_particle_pair[1].D, single2.pid_particle_pair[1].D
        D12 = D1 + D2

        assert (pos1 - single1.shell[1].position).sum() == 0
        pos2 = single2.shell[1].position
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
        com = self.applyBoundary( com )
        minShellSizeWithMargin = minShellSize + shellSizeMargin
        maxShellSize = min( self.getMaxShellSize(),
                            r0 * 100 + sigma + shellSizeMargin )

        if minShellSizeWithMargin >= maxShellSize:
            if __debug__:
                log.debug( '%s not formed: minShellSize >= maxShellSize' %
                       ( 'Pair( %s, %s )' % ( single1.pid_particle_pair[0], 
                                              single2.pid_particle_pair[0] ) ) )
            return None

        # Here, we have to take into account of the burst Singles in
        # this step.  The simple check for closest below could miss
        # some of them, because sizes of these Singles for this
        # distance check has to include SINGLE_SHELL_FACTOR, while
        # these burst objects have zero mobility radii.  This is not
        # beautiful, a cleaner framework may be possible.

        closest, closestShellDistance = None, numpy.inf
        for b in burst:
            if isinstance( b, Single ):
                bpos = b.shell[1].position
                d = self.distance( com, bpos ) \
                    - b.pid_particle_pair[1].radius * ( 1.0 + self.SINGLE_SHELL_FACTOR )
                if d < closestShellDistance:
                    closest, closestShellDistance = b, d

        if closestShellDistance <= minShellSizeWithMargin:
            if __debug__:
                log.debug( '%s not formed: squeezed by burst neighbor %s' %
                       ( 'Pair( %s, %s )' % ( single1.pid_particle_pair[0], 
                                              single2.pid_particle_pair[0]),
                         closest ) )
            return None

        assert closestShellDistance > 0
        c, d = self.getClosestObj( com, ignore=[ single1.domain_id, single2.domain_id ] )
        if d < closestShellDistance:
            closest, closestShellDistance = c, d

        if __debug__:
            log.debug( 'Pair closest neighbor: %s %g, minShellWithMargin %g' %
                   ( closest, closestShellDistance, minShellSizeWithMargin ) )

        assert closestShellDistance > 0

        if isinstance( closest, Single ):

            D_closest = closest.pid_particle_pair[1].D
            D_tot = D_closest + D12
            closestDistance = self.distance( com, closest.pid_particle_pair[1].position ) ##??

            closestMinRadius = closest.pid_particle_pair[1].radius
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
            assert isinstance( closest, ( Pair, Multi, None.__class__ ) )

            shellSize = closestShellDistance / SAFETY

        if shellSize <= minShellSizeWithMargin:
            if __debug__:
                log.debug( '%s not formed: squeezed by %s' %
                       ( 'Pair( %s, %s )' % ( single1.pid_particle_pair[0], 
                                              single2.pid_particle_pair[0]),
                         closest ) )
            return None


        d1 = self.distance(com, pos1)
        d2 = self.distance(com, pos2)

        if shellSize < max( d1 + single1.pid_particle_pair[1].radius *
                            ( 1.0 + self.SINGLE_SHELL_FACTOR ), \
                                d2 + single2.pid_particle_pair[1].radius * \
                                ( 1.0 + self.SINGLE_SHELL_FACTOR ) ) * 1.3:
            if __debug__:
                log.debug( '%s not formed: singles are better' %
                       'Pair( %s, %s )' % ( single1.pid_particle_pair[0], 
                                            single2.pid_particle_pair[0] ) )
            return None

        # 3. Ok, Pair makes sense.  Create one.
        shellSize = min( shellSize, maxShellSize )

        pair = self.createPair(single1, single2, shellSize)

        r0 = self.distance(pos1, pos2)

        pair.dt, pair.eventType, pair.reactingSingle, pair.activeCoordinate = \
            pair.determinePairEvent()
        assert pair.dt >= 0

        self.lastTime = self.t

        self.removeDomain( single1 )
        self.removeDomain( single2 )

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
    

    def formMulti(self, single, neighbors, dists):

        minShell = single.pid_particle_pair[1].radius * ( 1.0 + self.MULTI_SHELL_FACTOR )
        # Multis shells need to be contiguous.
        if dists[0] > minShell:
            return None

        neighbors = [neighbors[i] for i in (dists <= minShell).nonzero()[0]]

        closest = neighbors[0]

        # if the closest to this Single is a Single, create a new Multi
        if isinstance( closest, Single ):

            multi = self.createMulti()
            self.addToMulti( single, multi )
            self.removeDomain( single )
            for neighbor in neighbors:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )
            
            self.addMultiEvent( multi )

            return multi

        # if the closest to this Single is a Multi, reuse the Multi.
        elif isinstance( closest, Multi ):

            multi = closest
            if __debug__:
                log.info( 'multi merge %s %s' % ( single, multi ) )

            self.addToMulti( single, multi )
            self.removeDomain( single )
            for neighbor in neighbors[1:]:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )

            self.updateEvent( self.t + multi.dt, multi )

            return multi


        assert False, 'do not reach here'


    def addToMultiRecursive( self, obj, multi ):
        if isinstance( obj, Single ):
            if obj.pid_particle_pair[0] in multi.sim.particleList:  # Already in the Multi.
                return
            assert obj.isReset()
            objpos = obj.shell[1].position
            
            self.addToMulti( obj, multi )
            self.removeDomain( obj )
            self.removeEvent( obj )

            radius = obj.pid_particle_pair[1].radius *\
                ( 1.0 + self.MULTI_SHELL_FACTOR )
            neighbors = self.getNeighborsWithinRadiusNoSort( objpos, radius,
                                                             ignore=[obj.domain_id] )

            burst = self.burstNonMultis( neighbors )
            neighborDists = self.objDistanceArray( objpos, burst )
            neighbors = [ burst[i] for i in 
                          ( neighborDists <= radius ).nonzero()[0] ]

            for obj in neighbors:
                self.addToMultiRecursive( obj, multi )

        elif isinstance( obj, Multi ):
            if obj.sim.particleList.isdisjoint(multi.sim.particleList):
                self.mergeMultis(obj, multi)
                self.removeDomain(obj)
                self.removeEvent(obj)
            else:
                if __debug__:
                    log.debug( '%s already added. skipping.' % obj )
        else:
            assert False, 'do not reach here.'  # Pairs are burst

    def addToMulti( self, single, multi ):
        if __debug__:
            log.info( 'adding %s to %s' % ( single, multi ) )
        shellSize = single.pid_particle_pair[1].radius * \
            ( 1.0 + self.MULTI_SHELL_FACTOR )
        multi.addParticleAndShell(single.pid_particle_pair, shellSize)

    def mergeMultis( self, multi1, multi2 ):
        '''
        merge multi1 into multi2
        '''
        if __debug__:
            log.info( 'merging %s to %s' % ( multi1, multi2 ) )

        assert not multi1.sim.particleList[0] in multi2.sim.particleList

        for pid, sid in multi1.pid_shell_id_map.iteritems():
            # FIXME: shells should be renewed
            multi2.addParticleAndShell(
                multi1.sim.particleMatrix[pid],
                multi1.sim.shellMatrix[sid].radius)

    def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        '''
        Get neighbor domains within given radius.

        ignore: domain ids.
        '''

        result = self.shellMatrix.get_neighbors_within_radius(pos, radius)
        return [self.domains[did] for did in uniq(s[0][1].did for s in result) if did not in ignore]

    def getNeighbors(self, pos):
        return self.shellMatrix.get_neighbors(pos)

    def getNeighborsWithinRadius( self, pos, radius=numpy.inf, ignore=[] ):
        '''
        ignore: domain ids
        '''
        result = self.shellMatrix.get_neighbors_within_radius(pos, radius)

        seen = set(ignore)
        neighbors = []
        distances = []

        for item in result:
            did = item[0][1].did
            if not did in seen:
                seen.add(did)
                neighbors.append(self.domains[did])
                distances.append(item[1])
        return neighbors, distances

    def getClosestObj( self, pos, ignore=[] ):
        '''
        ignore: domain ids.
        '''

        result = self.shellMatrix.get_neighbors(pos)

        for item in result:
            if item[0][1].did not in ignore:
                return self.domains[item[0][1].did], item[1]

        return None, numpy.inf

    def objDistance( self, pos, obj ):
        dists = numpy.zeros( len( obj.shell_list ) )
        for i, shell_id_shell_pair in enumerate(obj.shell_list):
            dists[i] = self.distance( pos, shell_id_shell_pair[1].position ) - shell_id_shell_pair[1].radius
        return min( dists )

    def objDistanceArray( self, pos, objs ):
        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        return dists
            

    #
    # statistics reporter
    #

    def print_report(self, out=None):
        report = '''
t = %g
steps = %d 
\tSingle:\t%d\t(escape: %d, reaction: %d)
\tPair:\t%d\t(escape r: %d, R: %d, reaction pair: %d, single: %d)
\tMulti:\t%d\t(escape: %d, reaction: %d)
total reactions = %d
rejected moves = %d
'''\
            % (self.t, self.stepCounter,
               numpy.array(self.single_steps.values()).sum(),
               self.single_steps[EventType.SINGLE_ESCAPE],
               self.single_steps[EventType.SINGLE_REACTION],
               numpy.array(self.pair_steps.values()).sum(),
               self.pair_steps[EventType.IV_ESCAPE],
               self.pair_steps[EventType.COM_ESCAPE],
               self.pair_steps[EventType.PAIR_REACTION],
               self.pair_steps[EventType.SINGLE_REACTION],
               self.multi_steps[2], # total multi steps
               self.multi_steps[EventType.MULTI_ESCAPE],
               self.multi_steps[EventType.MULTI_REACTION],
               self.reactionEvents,
               self.rejectedMoves
               )

        print >> out, report

    #
    # consistency checkers
    #

    def checkObj( self, obj ):
        obj.check()

        for shell_id, shell in obj.shell_list:
            closest, distance = self.getClosestObj( shell.position,
                                                    ignore = [obj.domain_id] )

            assert shell.radius <= self.getUserMaxShellSize(),\
                '%s shell size larger than user-set max shell size' % \
                str( shell_id )

            assert shell.radius <= self.getMaxShellSize(),\
                '%s shell size larger than simulator cell size / 2' % \
                str( shell_id )

            assert distance - shell.radius >= 0.0,\
                '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                % ( str( obj ), str( closest ), shell.radius, distance,\
                        distance - shell.radius )

        return True

    def checkObjForAll( self ):
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            self.checkObj( obj )

    def checkEventStoichiometry( self ):
        population = 0
        for pool in self.particlePool.itervalues():
            population += len(pool)

        eventPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            eventPopulation += obj.multiplicity

        if population != eventPopulation:
            raise RuntimeError, 'population %d != eventPopulation %d' %\
                  ( population, eventPopulation )

    def checkShellMatrix( self ):
        if self.worldSize != self.shellMatrix.world_size:
            raise RuntimeError,\
                'self.worldSize != self.shellMatrix.worldSize'

        shellPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            shellPopulation += len(obj.shell_list)

        if shellPopulation != len(self.shellMatrix):
            raise RuntimeError,\
                'num shells (%d) != self.shellMatrix.size (%d)' % (shellPopulation, len(self.shellMatrix))
        
        self.shellMatrix.check()

    def check_domains(self):
        domains = set(self.domains.itervalues())
        for i in range(self.scheduler.getSize()):
            obj = self.scheduler.getEventByIndex(i).getArg()
            if obj not in domains:
                raise RuntimeError,\
                    '%s in EventScheduler not in self.domains' % obj
            domains.remove(obj)

        # self.domains always include a None  --> this can change in future
        if len(domains):
            raise RuntimeError,\
                'following domains in self.domains not in Event Scheduler: %s' \
                % str(tuple(domains))

    def checkPairPos( self, pair, pos1, pos2, com, radius ):
        particle1 = pair.single1.pid_particle_pair[1]
        particle2 = pair.single2.pid_particle_pair[1]

        oldCoM = com
        
        # debug: check if the new positions are valid:
        newDistance = distance_Simple( pos1, pos2 )
        particleRadius12 = particle1.radius + particle2.radius

        # check 1: particles don't overlap.
        if newDistance <= particleRadius12:
            if __debug__:
                log.info( 'rejected move: radii %g, particle distance %g',
                          ( particle1.radius + particle2.radius, newDistance ) )
            if __debug__:
                log.debug( 'DEBUG: pair.dt %g, pos1 %s, pos2 %s' %
                           ( pair.dt, str( pos1 ), str( pos2 ) ) )
            raise RuntimeError, 'New particles overlap'

        # check 2: particles within mobility radius.
        d1 = self.distance( oldCoM, pos1 ) + particle1.radius
        d2 = self.distance( oldCoM, pos2 ) + particle2.radius
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
        self.check_domains()
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

    def count_domains(self):
        '''
        Returns a tuple (# Singles, # Pairs, # Multis).
        '''

        numSingles = 0
        numPairs = 0
        numMultis = 0
        for d in self.domains.itervalues():
            if isinstance(d, Single):
                numSingles += 1
            elif isinstance(d, Pair):
                numPairs += 1
            elif isinstance(d, Multi):
                numMultis += 1
            else:
                raise RuntimeError, 'NO NOT GET HERE'

        return (numSingles, numPairs, numMultis)
