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
    CylindricalShellContainer,
    DomainIDGenerator,
    ShellIDGenerator,
    DomainID,
    ParticleContainer
    )

from surface import *

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
    def __init__(self, world):
        ParticleSimulatorBase.__init__(self, world)

        self.domainIDGenerator = DomainIDGenerator(0)
        self.shellIDGenerator = ShellIDGenerator(0)

        self.MULTI_SHELL_FACTOR = 0.05
        self.SINGLE_SHELL_FACTOR = 0.1

        self.isDirty = True
        self.scheduler = EventScheduler()

        self.userMaxShellSize = numpy.inf

        self.domains = {}

        self.reset()

    def getMatrixCellSize( self ):
        return self.containers[0].cell_size

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
                           EventType.IV_REACTION:0,
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

    def initialize( self ):
        ParticleSimulatorBase.initialize( self )

        self.scheduler.clear()
        self.containers = [SphericalShellContainer(self.world.world_size, 
                                                   self.world.matrix_size),
                           CylindricalShellContainer(self.world.world_size, 
                                                     self.world.matrix_size)]
        self.domains = {}

        singles = []
        for pid_particle_pair in self.world:
            single = self.createSingle(pid_particle_pair)
            if __debug__:
                log.debug("%s as single %s", pid_particle_pair[0], single.domain_id)
            singles.append(single)
        assert len(singles) == self.world.num_particles
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

        if __debug__:
            if self.scheduler.getSize() == 0:
                raise RuntimeError('No particles in scheduler.')

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

        if __debug__:
            if self.scheduler.getSize() == 0:
                raise RuntimeError('Zero particles left.')

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

    def createSingle( self, pid_particle_pair ):
        rt = self.getReactionRule1(pid_particle_pair[1].sid)
        domain_id = self.domainIDGenerator()
        shell_id = self.shellIDGenerator()

        # Get surface.
        species = self.world.get_species(pid_particle_pair[1].sid)
        surface = self.model.getSurface(species)

        # Create single. The type of the single that will be created depends 
        # on the surface this particle is on. Either SphericalSingle, 
        # PlanarSurfaceSingle, or CylindricalSurfaceSingle.
        TypeOfSingle = surface.DefaultSingle
        single = TypeOfSingle(domain_id, pid_particle_pair, shell_id, rt, 
                              surface)

        single.initialize(self.t)
        self.moveShell(single.shell_id_shell_pair)
        self.domains[domain_id] = single
        return single

    def createPair(self, single1, single2, CoM, r0, shellSize):
        assert single1.dt == 0.0
        assert single2.dt == 0.0

        rt = self.getReactionRule2(single1.pid_particle_pair[1].sid, single2.pid_particle_pair[1].sid)[ 0 ]

        domain_id = self.domainIDGenerator()
        shell_id = self.shellIDGenerator()

        pos1 = single1.shell.position
        pos2 = single2.shell.position

        # Get surface.
        species = self.world.get_species(single1.pid_particle_pair[1].sid)
        surface = self.model.getSurface(species)

        # Create pair. The type of the pair that will be created depends on 
        # the surface the particles are on. Either SphericalPair, 
        # PlanarSurfacePair, or CylindricalSurfacePair.
        TypeOfPair = surface.DefaultPair
        pair = TypeOfPair(domain_id, CoM, single1, single2, shell_id, 
                          r0, shellSize, rt, surface)

        pair.initialize( self.t )

        self.moveShell(pair.shell_id_shell_pair)
        self.domains[domain_id] = pair
        return pair

    def createMulti( self ):
        domain_id = self.domainIDGenerator()
        multi = Multi( domain_id, self )
        self.domains[domain_id] = multi
        if __debug__:
            try:
                # Option to make multis run faster for nicer visualization.
                multi.sim.dtFactor = DEFAULT_DT_FACTOR * self.bd_dt_factor
            except AttributeError:
                multi.sim.dtFactor = DEFAULT_DT_FACTOR 
        return multi

    def moveSingle(self, single, position, radius=None):
        self.moveSingleShell(single, position, radius)
        self.moveSingleParticle(single, position)

    def moveSingleShell(self, single, position, radius=None):
        if radius == None:
            # By default, don't change radius.
            radius = single.shell.radius

        # Reuse shell_id and domain_id.
        shell_id = single.shell_id
        domain_id = single.domain_id

        # Replace shell.
        shell = single.createNewShell(position, radius, domain_id)
        shell_id_shell_pair = (shell_id, shell) 

        single.shell_id_shell_pair = shell_id_shell_pair
        self.moveShell(shell_id_shell_pair)

    def moveSingleParticle(self, single, position):
        new_pid_particle_pair = (single.pid_particle_pair[0],
                          Particle(position,
                                   single.pid_particle_pair[1].radius,
                                   single.pid_particle_pair[1].D,
                                   single.pid_particle_pair[1].sid))
        single.pid_particle_pair = new_pid_particle_pair

        self.moveParticle(new_pid_particle_pair)

    def get_container(self, shell):
        if type(shell) is SphericalShell:
            return self.containers[0]
        elif type(shell) is CylindricalShell:
            return self.containers[1]

    def removeDomain( self, obj ):
        del self.domains[obj.domain_id]
        for shell_id, shell in obj.shell_list:
            container = self.get_container(shell)
            del container[shell_id]

    def moveShell(self, shell_id_shell_pair):
        shell = shell_id_shell_pair[1]
        container = self.get_container(shell)
        container.update(shell_id_shell_pair)

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
            # Don't schedule events in burst/propagatePair, because scheduling 
            # is different after a single reaction in firePair.
            self.addSingleEvent(single1)
            self.addSingleEvent(single2)
            self.removeEvent( obj )
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
        currentSurface = single.surface
        
        rt = single.drawReactionRule()

        if len( rt.products ) == 0:
            
            self.removeParticle(single.pid_particle_pair)

            self.lastReaction = Reaction( rt, [single.pid_particle_pair[1]], [] )

            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            if reactantSpeciesRadius < productSpecies.radius:
                self.clearVolume( oldpos, productSpecies.radius )

            if self.getParticlesWithinRadius( oldpos, productSpecies.radius,
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
                vector = \
                    currentSurface.randomVector(particleRadius12 *
                                                MINIMAL_SEPERATION_FACTOR)
            
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
                if (not self.getParticlesWithinRadius(newpos1, particleRadius1,
                                          ignore=[single.pid_particle_pair[0]])) and \
                   (not self.getParticlesWithinRadius(newpos2, particleRadius2,
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

    def propagateSingle(self, single):
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

        newpos = single.drawNewPosition(single.dt, single.eventType) 
        newpos = self.applyBoundary(newpos)

        if __debug__:
            log.debug( "propagate %s: %s => %s" % ( single, single.pid_particle_pair[1].position, newpos ) )

            if self.getParticlesWithinRadius(newpos,
                                 single.pid_particle_pair[1].radius,
                                 ignore=[single.pid_particle_pair[0]]):
                raise RuntimeError('propagateSingle: checkOverlap failed.')

        if(single.eventType == EventType.SINGLE_REACTION and
           single.eventType != EventType.BURST):
            # SINGLE_REACTION, and not a burst. No need to update, single is 
            # removed anyway.
            self.moveSingleParticle(single, newpos)
        else:
            # Todo. if isinstance(single, InteractionSingle):
            single.initialize(self.t)
            self.moveSingle(single, newpos, single.pid_particle_pair[1].radius)

    def fireSingle( self, single ):
        assert abs(single.dt + single.lastTime - self.t) <= 1e-18 * self.t

        # Reaction.
        if single.eventType == EventType.SINGLE_REACTION:
            if __debug__:
                log.info('fireSingle: eventType %s' % single.eventType)

            self.single_steps[single.eventType] += 1

            if __debug__:
                log.info( 'single reaction %s' % str( single ) )

            self.propagateSingle(single)

            try:
                self.removeDomain( single )
                self.fireSingleReaction( single )
            except NoSpace:
                self.reject_single_reaction(single)

            return

        # Propagate, if not reaction.
        single.eventType = EventType.SINGLE_ESCAPE
        if __debug__:
            log.info('fireSingle: eventType %s' % single.eventType)
        self.single_steps[single.eventType] += 1

        # Handle immobile case first.
        if single.getD() == 0:
            # no propagation, just calculate next reaction time.
            single.dt, single.eventType = single.determineNextEvent() 
            single.lastTime = self.t
            self.addSingleEvent(single)
            return
        
        if single.dt != 0.0:
            # Propagate this particle to the exit point on the shell.
            self.propagateSingle(single)

        singlepos = single.shell.position

        # (2) Clear volume.

        minShell = single.pid_particle_pair[1].radius * ( 1.0 + self.SINGLE_SHELL_FACTOR )

        intruders, closest, closestDistance = \
            self.get_intruders(singlepos, minShell, ignore=[single.domain_id, ])

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
                          ( single.shell_id_shell_pair, single.dt ) )

        self.addSingleEvent(single)
        return

    def reject_single_reaction(self, single):
        if __debug__:
            log.info( 'single reaction; placing product failed.' )
        self.domains[single.domain_id] = single
        self.moveShell(single.shell_id_shell_pair)
        self.rejectedMoves += 1
        single.initialize(self.t)
        self.addSingleEvent(single)

    def restoreSingleShells( self, singles ):
        for single in singles:
            assert single.isReset()
            c, d = self.getClosestObj( single.shell.position, ignore = [single.domain_id,] )

            self.updateSingle( single, c, d )
            self.updateEvent( self.t + single.dt, single )
            if __debug__:
                log.debug( 'restore shell %s %g dt %g closest %s %g' %
                       ( single, single.shell.radius, single.dt, c, d ) )

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

        singlepos = single.shell.position
        if isinstance( closest, Single ):
            closestpos = closest.shell.position
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

        single.dt, single.eventType = single.determineNextEvent()
        single.lastTime = self.t

    def firePair( self, pair ):
        assert self.checkObj( pair )

        single1 = pair.single1
        single2 = pair.single2
        particle1 = single1.pid_particle_pair
        particle2 = single2.pid_particle_pair
        pos1 = particle1[1].position
        pos2 = particle2[1].position
        
        if pair.eventType == EventType.IV_EVENT:
            # Draw actual pair event for iv at very last minute.
            r0 = self.distance(pos1, pos2)
            pair.eventType = pair.draw_iv_event_type(r0)

        self.pair_steps[pair.eventType] += 1

        if __debug__:
            log.info('firePair: eventType %s' % pair.eventType)


        oldCoM = pair.CoM

        # Four cases:
        #  1. Single reaction
        #  2. Pair reaction
        #  3a. IV escape
        #  3b. CoM escape

        #
        # 1. Single reaction
        #
        if pair.eventType == EventType.SINGLE_REACTION:
            reactingsingle = pair.reactingsingle

            if __debug__:
                log.info( 'pair: single reaction %s' % str( reactingsingle ) )

            if reactingsingle == single1:
                theothersingle = single2
            else:
                theothersingle = single1

            self.burstPair( pair )

            self.addSingleEvent( theothersingle )

            try:
                self.removeDomain( reactingsingle )
                self.fireSingleReaction( reactingsingle )
            except NoSpace:
                self.reject_single_reaction(reactingsingle)

            return
        
        #
        # 2. Pair reaction
        #
        if pair.eventType == EventType.IV_REACTION:
            if __debug__:
                log.info( 'reaction' )

            if len( pair.rt.products ) == 1:
                
                species3 = pair.rt.products[0]

                # calculate new R
                eventType = pair.eventType
                newCoM = pair.drawNewCoM(pair.dt, eventType)
                
                if __debug__:
                    shellSize = pair.get_shell_size()
                    assert self.distance( oldCoM, newCoM ) + species3.radius <\
                        shellSize

                newCoM = self.applyBoundary( newCoM )

                self.removeParticle(single1.pid_particle_pair)
                self.removeParticle(single2.pid_particle_pair)

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

        #
        # 3a. Escaping through a_r.
        # 3b. Escaping through a_R.
        #
        elif(pair.eventType == EventType.IV_ESCAPE or
             pair.eventType == EventType.COM_ESCAPE):
            dt = pair.dt
            eventType = pair.eventType
            single1, single2 = self.propagatePair(pair, dt, eventType)
            self.addSingleEvent(single1)
            self.addSingleEvent(single2)
        else:
            raise SystemError, 'Bug: invalid eventType.'

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

        oldpos = single.shell.position
        old_shell_size = single.get_shell_size()

        particleRadius = single.pid_particle_pair[1].radius

        # Override dt, burst happens before single's scheduled event.
        single.dt = self.t - single.lastTime
        # Override eventType. Always call gf.drawR on BURST.
        single.eventType = EventType.BURST
        self.propagateSingle(single)

        newpos = single.pid_particle_pair[1].position
        assert self.distance(newpos, oldpos) <= old_shell_size - particleRadius
        # Displacement check is in NonInteractionSingle.drawNewPosition.

        # Todo. if isinstance(single, InteractionSingle):
        self.updateEvent( self.t, single )

        assert single.shell.radius == particleRadius

    def burstPair( self, pair ):
        if __debug__:
            log.debug('burstPair: %s', pair)

        assert self.t >= pair.lastTime
        assert self.t <= pair.lastTime + pair.dt

        dt = self.t - pair.lastTime 
        # Override eventType. Always call sgf.drawR and pgf.drawR on BURST.
        eventType = EventType.BURST
        single1, single2 = self.propagatePair(pair, dt, eventType)

        return single1, single2

    def propagatePair(self, pair, dt, eventType):
        single1 = pair.single1
        single2 = pair.single2

        particle1 = single1.pid_particle_pair
        particle2 = single2.pid_particle_pair

        pos1 = particle1[1].position
        pos2 = particle2[1].position

        if dt > 0.0:
            D1 = particle1[1].D
            D2 = particle2[1].D

            pos2t = self.world.cyclic_transpose(pos2, pos1)
            oldInterParticle = pos2t - pos1
            r0 = self.distance(pos1, pos2)
            assert feq(r0, length(oldInterParticle))

            oldCoM = pair.CoM

            newpos1, newpos2 = pair.drawNewPositions(dt, r0, 
                                                     oldInterParticle, 
                                                     eventType)

            newpos1 = self.applyBoundary(newpos1)
            newpos2 = self.applyBoundary(newpos2)
            assert not self.getParticlesWithinRadius(newpos1, particle1[1].radius,
                                         ignore=[particle1[0], particle2[0]])
            assert not self.getParticlesWithinRadius(newpos2, particle2[1].radius,
                                         ignore=[particle1[0], particle2[0]])
            assert self.checkPairPos(pair, newpos1, newpos2, oldCoM,\
                                         pair.get_shell_size())
        else:
            newpos1 = particle1[1].position
            newpos2 = particle2[1].position

        if __debug__:
            log.debug("firePair: #1 { %s: %s => %s }" % (single1, pos1, newpos1))
            log.debug("firePair: #2 { %s: %s => %s }" % (single2, pos2, newpos2))

        single1.initialize( self.t )
        single2.initialize( self.t )
        
        self.removeDomain( pair )
        assert single1.domain_id not in self.domains
        assert single2.domain_id not in self.domains
        self.domains[single1.domain_id] = single1
        self.domains[single2.domain_id] = single2
        self.moveSingle(single1, newpos1, particle1[1].radius)
        self.moveSingle(single2, newpos2, particle2[1].radius)

        if __debug__:
            container = self.get_container(single1.shell)
            assert container[single1.shell_id].radius == single1.shell.radius
            assert container[single2.shell_id].radius == single2.shell.radius

            if type(single1.shell) is CylindricalShell:
                assert container[single1.shell_id].size == single1.shell.size
                assert container[single2.shell_id].size == single2.shell.size

        assert single1.shell.radius == particle1[1].radius
        assert single2.shell.radius == particle2[1].radius

        assert self.checkObj(single1)
        assert self.checkObj(single2)

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

        assert (pos1 - single1.shell.position).sum() == 0
        pos2 = single2.shell.position
        r0 = self.distance(pos1, pos2)
        distanceFromSigma = r0 - sigma
        assert distanceFromSigma >= 0, \
            ('distanceFromSigma (pair gap) between %s and %s = %g < 0' %
             (single1, single2, distanceFromSigma))

        shellSize1 = r0 * D1 / D12 + radius1
        shellSize2 = r0 * D2 / D12 + radius2
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
        com = self.world.calculate_pair_CoM( pos1, pos2, D1, D2 )
        com = self.applyBoundary( com )
        minShellSizeWithMargin = minShellSize + shellSizeMargin
        maxShellSize = min( self.getMaxShellSize(),
                            distanceFromSigma * 100 + sigma + shellSizeMargin )

        if minShellSizeWithMargin >= maxShellSize:
            if __debug__:
                log.debug('%s not formed: minShellSize %g >= maxShellSize %g' %
                          ('Pair( %s, %s )' % (single1.pid_particle_pair[0], 
                                               single2.pid_particle_pair[0]),
                           minShellSizeWithMargin, maxShellSize))
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
                bpos = b.shell.position
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

        pair = self.createPair(single1, single2, com, r0, shellSize)

        pair.dt, pair.eventType, pair.reactingsingle = \
            pair.determineNextEvent(r0)

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
            log.info( '%s, dt=%g, r0=%g, shell=%g,' %
                  ( pair, pair.dt, r0, shellSize ) + 
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
            objpos = obj.shell.position
            
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

            some_particle_of_multi1 = None
            try:
                some_particle_of_multi1 = iter(multi1.sim.particleList).next()
            except:
                pass
            assert some_particle_of_multi1 not in multi2.sim.particleList

        for pid, sid in multi1.pid_shell_id_map.iteritems():
            # FIXME: shells should be renewed
            multi2.addParticleAndShell(
                (pid, multi1.sim.particleMatrix[pid]),
                multi1.sim.sphere_container[sid].radius)

    def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        """Get neighbor domains within given radius.

        ignore: domain ids.

        Only returns neighbors, not the distances towards their shells. Can 
        for example be used to try to clear all objects from a certain volume.

        """
        neighbors = []
        for container in self.containers:
            result = container.get_neighbors_within_radius(pos, radius)
            # result = [((shell_id_shell_pair), distance), ]
            # Since a domain can have more than 1 shell (multis for example), 
            # and for each shell there is an entry in the shell container, we 
            # make sure each domain occurs only once in the returned list 
            # here.
            neighbors.extend(self.domains[did]
                             for did in uniq(s[0][1].did for s in result)
                                     if did not in ignore)
        return neighbors

    def get_intruders(self, position, radius, ignore):
        intruders = []   # intruders are domains within radius
        closest_domain = None   # closest domain, excluding intruders.
        closest_distance = numpy.inf # distance to the shell of the closest.

        seen = set(ignore)
        for container in self.containers:
            neighbors = container.get_neighbors(position)
            for n in neighbors:
                domain_id = n[0][1].did
                distance = n[1]
                if distance > radius:
                    if distance < closest_distance:
                        # This is domain (the first one for this container) 
                        # that has a shell that is more than radius away from 
                        # pos.  If it is closer than the closest such one we 
                        # found so far: store it. Always break out of the 
                        # inner for loop and check the other containers.
                        closest_domain = self.domains[domain_id]
                        closest_distance = distance
                        break
                    else:
                        break
                elif domain_id not in seen:
                    # Since a domain can have more than 1 shell (multis for 
                    # example), and for each shell there is an entry in the 
                    # shell container, we make sure each domain occurs only 
                    # once in the returned list here.
                    seen.add(domain_id)
                    intruders.append(self.domains[domain_id])

        return intruders, closest_domain, closest_distance

    def getClosestObj( self, pos, ignore=[] ):
        '''
        ignore: domain ids.

        '''
        closest_domain = None
        closest_distance = numpy.inf

        for container in self.containers:
            result = container.get_neighbors(pos)

            for shell_id_shell_pair, distance in result:
                domain_id = shell_id_shell_pair[1].did 

                if domain_id not in ignore and distance < closest_distance:
                    domain = self.domains[domain_id]
                    closest_domain, closest_distance = domain, distance
                    # Found yet a closer domain. Break out of inner for loop 
                    # and check other containers.
                    break   

        return closest_domain, closest_distance

    def objDistance( self, pos, obj ):
        return min(self.world.distance(shell, pos) for i, (_, shell) in enumerate(obj.shell_list))

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
               self.pair_steps[EventType.IV_REACTION],
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
            if(type(obj) is CylindricalSurfaceSingle or
               type(obj) is CylindricalSurfacePair):
                shell_size = shell.size
            else:
                shell_size = shell.radius

            assert shell_size <= self.getUserMaxShellSize(),\
                '%s shell size larger than user-set max shell size' % \
                str( shell_id )

            assert shell_size <= self.getMaxShellSize(),\
                '%s shell size larger than simulator cell size / 2' % \
                str( shell_id )

            assert distance - shell_size >= 0.0,\
                '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                % ( str( obj ), str( closest ), shell_size, distance,\
                        distance - shell_size )

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
        for container in self.containers:
            if self.world.world_size != container.world_size:
                raise RuntimeError,\
                    'self.world.world_size != container.world_size'

        shellPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            shellPopulation += len(obj.shell_list)
  
        matrix_population = sum(len(container) for container in self.containers)
        if shellPopulation != matrix_population:
            raise RuntimeError(
                'num shells (%d) != matrix population (%d)' % 
                (shellPopulation, matrix_population))
        
        for container in self.containers:
            container.check()

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
        newDistance = distance( pos1, pos2 )
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
