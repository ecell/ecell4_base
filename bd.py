#!/usr/env python

import weakref

import math

import random
import numpy
from numpy.random import uniform

from utils import *
from surface import *

from gfrdbase import *
import _gfrd

import logging

log = logging.getLogger( 'epdp' )

DEFAULT_DT_FACTOR = 1e-5

class BDSimulatorCoreBase( object ):
    '''
    BDSimulatorCore borrows the following from the main simulator:
    - speciesList
    - reactionTypes list (both 1 and 2)
    
    '''

    def __init__( self, main ):
        self.main = weakref.proxy( main )

        self.particleList = set()

        self.t = 0.0
        self.dt = 0.0

        self.dtFactor = DEFAULT_DT_FACTOR

        self.stepCounter = 0

        self.reactionEvents = 0

        self.lastReaction = None

        self.P_acct = {}

    def initialize( self ):
        self.determineDt()

    #@staticmethod  # requires python 2.4 or later.
    def calculateBDDt( self, speciesList, factor ):
        D_list = []
        radius_list = []
        for species in speciesList:
            if self.main.particlePool[species.serial]:
                D_list.append( species.D )
                radius_list.append( species.radius )
        D_max = max( D_list ) * 2  # max relative diffusion speed
        radius_min = min( radius_list )
        sigma_min = radius_min * 2

        dt = factor * sigma_min ** 2 / D_max  
        if __debug__:
            log.debug( 'bd dt = %g' % dt )

        return dt


    def clearParticleList( self ):
        self.particleList = set()

    def addToParticleList( self, pid ):
        self.particleList.add( pid )

    def removeFromParticleList( self, pid ):
        self.particleList.remove( pid )

    def getNextTime( self ):
        return self.t + self.dt

    def stop( self, t ):
        # dummy
        self.t = t

    def determineDt( self ):
        self.dt = self.calculateBDDt( self.main.speciesList.itervalues(), self.dtFactor )

    def getP_acct( self, rt, D, sigma ):
        try:
            return self.P_acct[ rt ]
        except KeyError:
            I = _gfrd.I_bd( sigma, self.dt, D )
            p = rt.k * self.dt / ( I * 4.0 * numpy.pi )
            if not 0.0 <= p < 1.0:
                raise RuntimeError,\
                    'Invalid acceptance ratio (%s) for reaction %s.' \
                    % ( p, rt )
            self.P_acct[ rt ] = p
            return p

    def step( self ):
        self.stepCounter += 1
        self.lastReaction = None

        self.propagate()

        self.t += self.dt

    def propagate( self ):
        self.particlesToStep = list(self.particleList)
        random.shuffle(self.particlesToStep)
        while self.particlesToStep:
            pid = self.particlesToStep.pop() # take the last one
            pid_particle_pair = (pid, self.main.particleMatrix[pid])
            sid = pid_particle_pair[1].sid

            rt1 = self.attemptSingleReactions(sid)
            if rt1:
                try:
                    self.fireReaction1( particle, rt1 )
                except NoSpace:
                    if __debug__:
                        log.info( 'fireReaction1 rejected.' )
                continue

            D = pid_particle_pair[1].D
            if D == 0.0:
                continue

            displacement = drawR_free( self.dt, D )

            newpos = pid_particle_pair[1].position + displacement
            newpos = self.main.applyBoundary(newpos)

            neighbors = self.getParticlesWithinRadiusNoSort(
                newpos, pid_particle_pair[1].radius,
                ignore=[pid_particle_pair[0]] )
            if neighbors:

                if len( neighbors ) >= 2:
                    if __debug__:
                        log.info( 'collision two or more particles; move rejected' )
                    continue

                closest = neighbors[0]

                rt = self.main.getReactionRule2(sid, closest[1].sid)[0]

                if rt.k != 0.0:
                    radius12 = pid_particle_pair[1].radius + closest[1].radius
                    D12 = D + closest[1].D

                    p = self.getP_acct( rt, D12, radius12 )

                    rnd = uniform()

                    if p > rnd:
                        if __debug__:
                            log.info( 'fire reaction2' )
                        try:
                            self.fireReaction2( pid_particle_pair, closest, rt )
                        except NoSpace:
                            if __debug__:
                                log.info( 'fireReaction2 move rejected' )
                        continue

                else:
                    if __debug__:
                        log.info( 'collision move rejected' )

                continue

            try:
                self.clearVolume(newpos, pid_particle_pair[1].radius, ignore=[pid_particle_pair[0]])
                self.moveParticle((pid_particle_pair[0],
                                   _gfrd.Particle(newpos,
                                            pid_particle_pair[1].radius,
                                            pid_particle_pair[1].D,
                                            pid_particle_pair[1].sid)))
            except NoSpace:
                if __debug__:
                    log.info( 'propagation move rejected.' )

    def attemptSingleReactions(self, sid):
        reactionTypes = self.main.getReactionRule1(sid)
        if not reactionTypes:
            return None  # no reaction

        rnd = uniform() / self.dt

        # handle the most common case efficiently.
        if len( reactionTypes ) == 1:  
            if reactionTypes[0].k >= rnd:
                return reactionTypes[0]
            else:
                return None

        # if there are more than one possible reaction types..
        k_array = numpy.add.accumulate( [ rt.k for rt in reactionTypes ] )

        if k_array[-1] < rnd:
            return None

        i = numpy.searchsorted( k_array, rnd )

        return reactionTypes[i]

    def fireReaction1( self, particle, rt ):
        oldpos = particle.pos.copy()

        if len( rt.products ) == 0:
            
            self.removeParticle( particle )

            self.lastReaction = Reaction( rt, [particle], [] )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]
            radius = productSpecies.radius

            if not self.checkOverlap( oldpos, radius,
                                      ignore = [ particle, ] ):
                if __debug__:
                    log.info( 'no space for product particle.' )
                raise NoSpace()

            self.clearVolume( oldpos, radius, ignore = [ particle ] )
                
            self.removeParticle( particle )
            newparticle = self.createParticle( productSpecies, oldpos )

            self.lastReaction = Reaction( rt, [particle], [newparticle] )

            
        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]
            
            D1 = productSpecies1.D
            D2 = productSpecies2.D
            D12 = D1 + D2
            
            radius1 = productSpecies1.radius
            radius2 = productSpecies2.radius
            radius12 = radius1 + radius2

            for i in range( 100 ):

                rnd = uniform()
                pairDistance = drawR_gbd( rnd, radius12, self.dt, D12 )

                unitVector = randomUnitVector()
                vector = unitVector * pairDistance # * (1.0 + 1e-10) # safety
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                newpos1 = oldpos + vector * ( D1 / D12 )
                newpos2 = oldpos - vector * ( D2 / D12 )
                #FIXME: check surfaces here
            
                newpos1 = self.main.applyBoundary( newpos1 )
                newpos2 = self.main.applyBoundary( newpos2 )

                # accept the new positions if there is enough space.
                if self.checkOverlap( newpos1, radius1,
                                      ignore = [ particle, ]) and \
                                      self.checkOverlap( newpos2, radius2,
                                                         ignore = 
                                                         [ particle, ]):
                    break
            else:
                if __debug__:
                    log.info( 'no space for product particles.' )
                raise NoSpace()

            self.clearVolume( newpos1, radius1, ignore = [ particle ] )
            self.clearVolume( newpos2, radius2, ignore = [ particle ] )

            # move accepted
            self.removeParticle( particle )

            newparticle1 = self.createParticle( productSpecies1, newpos1 )
            newparticle2 = self.createParticle( productSpecies2, newpos2 )

            self.lastReaction = Reaction( rt, [particle], 
                                          [newparticle1, newparticle2] )

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1

    def fireReaction2( self, pid_particle_pair1, pid_particle_pair2, rt ):
        if len(rt.products) == 1:
            productSpecies = rt.products[0]

            D1 = pid_particle_pair1[1].D
            D2 = pid_particle_pair2[1].D

            pos2t = cyclic_transpose(pid_particle_pair1[1].position,
                                     pid_particle_pair2[1].position,
                                     self.main.worldSize )
            newPos = (D2 * pid_particle_pair1[1].position + D1 * pos2t) / (D1 + D2)
            newPos = self.main.applyBoundary(newPos)

            if not self.checkOverlap(newPos, productSpecies.radius,
                                     ignore=[pid_particle_pair1[0],
                                             pid_particle_pair2[0]]):
                raise NoSpace()
            self.clearVolume(newPos, productSpecies.radius,
                             ignore=[pid_particle_pair1[0],
                                     pid_particle_pair2[0]])

            self.removeParticle(pid_particle_pair1)
            self.removeParticle(pid_particle_pair2)
            newparticle = self.createParticle(productSpecies.serial, newPos)

            try:
                self.particlesToStep.remove(pid_particle_pair2[0])
            except ValueError:  
                pass     # particle2 already stepped, which is fine.

            self.reactionEvents += 1

            self.lastReaction = Reaction(
                rt, [pid_particle_pair1, pid_particle_pair2], 
                [newparticle])

            return
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'

    def check( self ):
        # particles don't overlap

        for pid in self.particleList:
            particle = self.main.particleMatrix[pid]
            assert self.checkOverlap( particle.position, particle.radius,
                                      ignore=[pid,] )


class BDSimulatorCore( BDSimulatorCoreBase ):
    def __init__( self, main ):
        BDSimulatorCoreBase.__init__( self, main )


    def getParticlesWithinRadius( self, pos, radius, ignore=[] ):
        return self.main.getParticlesWithinRadius( pos, radius, ignore )

    def getParticlesWithinRadiusNoSort( self, pos, radius, ignore=[] ): 
        return self.main.getParticlesWithinRadiusNoSort( pos, radius, ignore )

    def initialize( self ):
        BDSimulatorCoreBase.initialize( self )

        self.updateParticleList()

    def updateParticleList( self ):
        self.clearParticleList()
        for s in self.main.speciesList.itervalues():
            for pid in self.main.particlePool[s.serial]:
                self.addToParticleList(pid) 

    def addParticle(self, pid_particle_pair):
        self.main.addParticle(pid_particle_pair)
        self.addToParticleList(pid_particle_pair[0])

    def moveParticle(self, pid_particle_pair):
        return self.main.moveParticle(pid_particle_pair)

    def checkOverlap( self, pos, radius, ignore=() ):
        return self.main.checkOverlap( pos, radius, ignore )

    def removeParticle(self, pid_particle_pair):
        self.main.removeParticle(pid_particle_pair)
        self.removeFromParticleList(pid_particle_pair[0])

    def createParticle(self, species, pos):
        particle = self.main.createParticle( species, pos )
        self.addToParticleList( particle )

    def clearVolume( self, pos, radius, ignore=[] ):
        '''
        This method is a customization point for implementing
        BD in protective domains.
        '''
        pass


class BDSimulator( ParticleSimulatorBase ):
    def __init__( self ):
        ParticleSimulatorBase.__init__( self )
        self.core = BDSimulatorCore( self )
        self.isDirty = True

    def t( self ):
        return self.core.t

    def sett( self, t ):
        self.core.t = t

    t = property(t, sett)

    def getDt( self ):
        return self.core.dt

    def getStepCounter( self ):
        return self.core.stepCounter

    dt = property( getDt )
    stepCounter = property( getStepCounter )


    def initialize( self ):
        self.core.initialize()
        self.isDirty = False

    def getNextTime( self ):
        return self.core.t + self.core.dt

    def reset( self ):
        # DUMMY
        self.core.t=0

    def stop( self, t ):
        # dummy
        self.core.stop( t )

    def step( self ):
        self.reactionType = None

        if self.isDirty:
            self.initialize()

        self.core.step()

        if __debug__:
            log.info( '%d: t=%g dt=%g, reactions=%d, rejectedMoves=%d' %
                  ( self.stepCounter, self.t, self.dt, self.reactionEvents,
                    self.rejectedMoves ) )

    def check( self ):
        pass
