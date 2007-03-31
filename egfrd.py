#!/usr/env python


import math
import random

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

        #FIXME: si and i are fragile
        self.particle = particle
        self.sim = sim
        self.lastTime = 0.0
        self.dt = 0.0
        self.closest = (-1, -1)
        self.eventID = None

    def fire( self ):
        dt = self.sim.fireSingle( self )
        return dt

    def update( self, t ):

        self.lastTime = t
        
        neighbors, drs = self.sim.getNeighbors( self.particle.getPos() )
        closest = neighbors[1]
        dr = drs[1] - self.particle.species.radius
        closestParticle = Particle( closest[0], index=closest[1] )
        closestSingle = self.sim.findSingle( closestParticle )
        
        self.closest = closestSingle
        
        #print 'dr', dr
        #FIXME: take different D into account
        self.setDr( dr * .5 )
        
        self.dt = self.calculateFirstPassageTime()


    def isDependentOn( self, event ):
        #print event
        return False

    def setPosition( self, pos ):
        self.particle.species.pool.positions[self.particle.i] = pos

    def getPosition( self ):
        return self.particle.species.pool.positions[self.particle.i]

    def setDr( self, dr ):
        pool = self.particle.species.pool
        pool.drs[ pool.getIndex( self.particle.serial ) ] = dr

    def getDr( self ):
        pool = self.particle.species.pool
        return pool.drs[ pool.getIndex( self.particle.serial ) ]

    def calculateFirstPassageTime( self ):
        
        species = self.particle.species
        fpgf = FirstPassageGreensFunction( species.D )
        rnd = random.random()
        dr = self.getDr()
        if dr <= 0.0:
            raise 'dr <= 0.0: %s' % str(dr)
        dt = fpgf.drawTime( rnd, dr )
        print dt
        if dt <= 0.0:
            raise 'dt <= 0.0: %s' % str(dt)
        return dt


    def __str__( self ):
        return str(self.particle) + str(self.getDr())



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
        self.pgf = FirstPassagePairGreensFunction( D12, rt.k,
                                                   self.sigma )

        self.eventID = None



    def __del__( self ):
        print 'del', str( self )



    '''
    Calculate and return Center of Mass (== CoM) of this pair.
    '''

    def getCoM( self ):
        particle1 = self.single1.particle
        particle2 = self.single2.particle
        
        D1 = particle1.species.D
        D2 = particle2.species.D


        pos1 = particle1.getPos()
        pos2 = particle2.getPos()
        
        #com = sqrtD2D1 * pos1 + self.sqrtD1D2 * pos2
        com = ( pos1 + self.sqrtD1D2 * pos2 ) * .5

        return com


    def newPositions( self, newCoM, newInterParticle, oldInterParticle ):

        # Now I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotationAxis2 = crossproduct( [ 0,0,1 ], interParticle )
        
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
        species2 = particle1.species

        pos1 = particle1.getPos()
        pos2 = particle2.getPos()

        D1 = species1.D
        D2 = species2.D

        oldInterParticle = pos2 - pos1

        # 1. now we handle the reaction case first.
        if self.eventType == EventType.REACTION:

            if len( self.rt.products ) == 1:
                
                species3 = self.rt.products[0]

                if D1 == 0.0:
                    newR = pos1
                elif D2 == 0.0:
                    newR = pos2
                else:
                    R0 = self.getCoM()
                    dR = p_free( D1/4, D2/4, self.dt )
                    newR = R0 + dR
                
                
                #FIXME: SURFACE
                newR = self.sim.applyBoundary( newR )

                species1.removeParticleBySerial( particle1.serial )
                species2.removeParticleBySerial( particle2.serial )

                particle = self.sim.createParticle( species3, newR )


                single = self.sim.createSingle( particle )
                self.sim.createSingleEvent( single )

                # self.sim.scheduler.removeEvent( self.eventID )
                # returning -1 will make the scheduler removing this event.
                return -1

            else:
                raise NotImplementedError,\
                      'num products >= 2 not supported yet.'

        # 2. escape cases.

        # 2.1 escape r
        if self.eventType == EventType.ESCAPE:

            print 'escape r'

            rnd = numpy.random.uniform( size=2 )

            # calculate new R
            
            r_R = self.sgf.drawR( rnd[0], self.dt, self.a_R )
            
            displacement_R_S = numpy.array( [ r_R,
                                              random.uniform( 0.0, Pi ),
                                              random.uniform( 0.0, 2*Pi ) ] )
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = self.getCoM() + displacement_R


            # calculate new r
            print ( rnd[1], self.a_r, self.r0, self.dt )
            theta_r = self.pgf.drawTheta( rnd[1], self.a_r*.09, self.r0, self.dt )
            phi_r = random.uniform( 0.0, 2*Pi )
            newInterParticleS = numpy.array( [ self.a_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = self.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )



            #return 0.0
            #raise NotImplementedError,'ESCAPE'


        # 2.2 escape R
        elif self.eventType == 2:

            print 'escape R'

            # calculate new R
            displacement_R_S = numpy.array( [ self.a_R,
                                              random.uniform( 0.0, Pi ),
                                              random.uniform( 0.0, 2*Pi ) ] )
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = self.getCoM() + displacement_R


            # calculate new r
            rnd = numpy.random.uniform( size = 2 )
            r = self.pgf.drawR( rnd[0], self.a_r, self.dt )
            print ( rnd[1], r, self.r0, self.dt )
            theta_r = self.pgf.drawTheta( rnd[1], r, self.r0, self.dt )
            phi_r = random.uniform( 0.0, 2*Pi )
            newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = self.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
                
            # raise NotImplementedError,'ESCAPE2'  # escape R
        else:
            raise SystemError, 'Bug: invalid eventType.'

        particle1.setPos( newpos1 )
        particle2.setPos( newpos2 )

        # here decide whether this pair still continues or breaks up

        dt, eventType = self.nextEvent( self.dr )
        return dt


    def update( self, t ):
        print 'update'

    def isDependentOn( self, event ):
        #print event
        return False

    def __str__( self ):
        return str(self.single1.particle) + ' ' + str(self.single2.particle)


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

        self.formPairs()

        #debug
        #self.checkShellForAll()

        #self.scheduler.updateAllEventDependency()

        self.isDirty = False



    def step( self ):

        if self.isDirty:
            self.initialize()

        self.lastEvent = self.scheduler.getTopEvent()[1]

        self.t = self.scheduler.getTime()

        self.scheduler.step()

        self.dt = self.scheduler.getTime() - self.t
        
        # if the same single stepped in the last n steps,
        # reinitialize everything.
        # FIXME: don't need to initialize everything.
        #        (1) recalculate dr to the closest with
        #            its dr shrunken.
        #        (2) new dr is
        #            min( dr to the second closest,
        #                 dr to the closest, with its dr
        #                 shrunken )

        nextEvent = self.scheduler.getTopEvent()[1]
        #print self.hoggerCounter, nextEvent, self.lastEvent
        if self.lastEvent is nextEvent:
            self.hoggerCounter += 1
        else:
            self.hoggerCounter = 0

        if self.hoggerCounter >= 10: # or self.dt < 1e-15:
                print 'reinitialize'
                self.hoggerCounter = 0
                self.initialize()
                #nextEvent = self.scheduler.getTopEvent()[1]
                self.dt = self.scheduler.getTime() - self.t

        #if self.dt == 0.0:
        #    raise 'dt=0'


        print 'dt', self.dt,\
              'reactions', self.reactionEvents,\
              'rejected moves', self.rejectedMoves,\
              'hogger counter', self.hoggerCounter
        print ''
        


    def fireSingle( self, single ):

        species = single.particle.species

        displacementS = numpy.array( ( single.getDr(),
                                       numpy.random.uniform( 0.0, Pi ),
                                       numpy.random.uniform( 0.0, 2*Pi ) ) )
        displacement = sphericalToCartesian( displacementS )

        #self.checkShellForAll()
        
        pos = single.particle.getPos()
        pos += displacement

        # BOUNDARY
        self.applyBoundary( pos )
        
        #print 'displacement', length(displacement), single.getDr()

        neighbors, distances = self.getNeighborShells( pos )
        newdr = distances[1] - single.particle.species.radius
        
        single.closest = self.findSingle( Particle( neighbors[1][0],
                                                    neighbors[1][1] ) )

        #print 'newdr', newdr
        if newdr <= 0:
            print newdr, single.closest
            raise 'Fatal newdr <= 0'

        single.setDr( newdr )

        #print single, single.closest
        #debug
        #self.checkShellForAll()

        return single.calculateFirstPassageTime()


    def initializeSingleMap( self ):

        self.singleMap = {}

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                single = self.createSingle( particle )
                self.singleMap[ ( species, particle.serial ) ] = single        

    def findSingle( self, particle ):
        return self.singleMap.get( ( particle.species, particle.serial ) )

    def createSingle( self, particle ):
        single = Single( self, particle )
        single.setDr( 0.0 )
        return single

    def initializeSingles( self ):

        for single in self.singleMap.values():
            self.createSingleEvent( single )
        
        #FIXME: here perhaps a second pass to get drs maximally large.


    def createSingleEvent( self, single ):
        single.update( self.t )
        nextt = single.lastTime + single.dt
        self.scheduler.addEvent( nextt, single )


    def formPairs( self ):
        #pass
        self.formPairsModestly()

    def formPairsModestly( self ):

        for single in self.singleMap.values():

            #neighbors, drs = self.getNeighbors( single.particle.getPos() )
            #closest = Particle( neighbors[1][0], neighbors[1][1] )
            closest = single.closest

            if single.particle < closest.particle:
                continue

            partnerNeighbors, partnerDrs = self.getNeighbors( closest.particle.getPos())
            partnerClosest = Particle( partnerNeighbors[1][0],
                                       index=partnerNeighbors[1][1] )

            if single.particle == partnerClosest:
                
                pos1 = single.particle.getPos()
                pos2 = closest.particle.getPos()
                r0 = self.distance( pos1, pos2 )

                pair = self.createPair( single, closest )
                com = pair.getCoM()

                neighbors, drs = self.getNeighbors( com, 3 )
                pairDr = drs[2]

                if pairDr < r0:   # this happens with a small probability
                    print 'pairDr < r0', pairDr, r0, pairDr - r0
                    #raise ''
                    break
                
                if pairDr > r0 * 5:
                    pairDr = r0 * 5

                nextEvent = pair.nextEvent( pairDr )
                dt = nextEvent[0]
                #eventType = nextEvent[1]

                if pair.single1.eventID != None:
                    self.scheduler.removeEvent( pair.single1.eventID )
                if pair.single2.eventID != None:
                    self.scheduler.removeEvent( pair.single2.eventID )

                pair.eventID = self.scheduler.addEvent( self.t + dt, pair ) 
                

    def formPairsGreedily( self ):

        for single in self.singleMap.values():
            print single, single.closest, single.getDr()

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

            
        
    def checkShell( self, single ):
        neighbors, drs = self.getNeighborShells( single.particle.getPos() )
        closest, distance = neighbors[1], drs[1]
        distance -= single.particle.species.radius
        if single.getDr() - distance >= 1e-18:
            dr = single.getDr()
            print single.particle, closest, dr, distance, dr - distance
            raise RuntimeError, 'Fatal: shells overlap.'

    def checkShellForAll( self ):
        for single in self.singleMap.values():
            self.checkShell( single )



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

            indices, distances = self.getNeighborsInSpecies( pos,
                                                             positions,
                                                             n )
            distances -= species.radius

            topDistances.extend( distances )
            topNeighbors.extend( [ ( species, i ) for i in indices ] )

        topargs = numpy.argsort( topDistances )[:n]
        topDistances = numpy.take( topDistances, topargs )
        topNeighbors = [ topNeighbors[arg] for arg in topargs ]

        return topNeighbors, topDistances


    def getNeighborsInSpecies( self, position1, positions, n=2 ):

        distances = self.distanceSqArray( position1, positions )
        
        indices = distances.argsort()[:n]
        distances = distances.take( indices )
        distances = numpy.sqrt( distances )
        
        return indices, distances


    def getNeighborShells( self, pos, n=2, speciesList=None ):

        topNeighbors = []
        topDistances = []

        if speciesList == None:
            speciesList = self.speciesList.values()

        for species in speciesList:

            # empty
            if species.pool.size == 0:
                continue

            positions = species.pool.positions
            indices, distances = self.getNeighborShellsInSpecies( pos,
                                                                  species,
                                                                  positions,
                                                                  n )
            topDistances.extend( distances )
            topNeighbors.extend( [ ( species, i ) for i in indices ] )

        topargs = numpy.argsort( topDistances )[:n]
        topDistances = numpy.take( topDistances, topargs )
        topNeighbors = [ topNeighbors[arg] for arg in topargs ]

        return topNeighbors, topDistances


    def getNeighborShellsInSpecies( self, position1, species, positions, n=2 ):

        distances = self.distanceSqArray( position1, positions )

        distances = numpy.sqrt( distances )
        drs = species.pool.drs
        distances -= drs

        indices = distances.argsort()[:n]
        distances = distances.take( indices )
        return indices, distances

