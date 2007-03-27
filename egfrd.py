#!/usr/env python


import math
import random

import numpy
#import scipy
#import scipy.optimize


from utils import *
from surface import *
import gfrdfunctions
import _gfrd

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
        self.dt = 0.0
        self.closest = (-1, -1)
        self.eventID = None


    def fire( self ):
        dt = self.sim.fireSingle( self )
        return dt

    def update( self, t ):
        print 'update'

    def isDependentOn( self, event ):
        #print event
        return False

    def setPosition( self, pos ):
        self.particle.species.pool.positions[self.particle.i] = pos

    def getPosition( self ):
        return self.particle.species.pool.positions[self.particle.i]

    def setDr( self, dr ):
        pool = self.particle.species.pool
        pool.drs[ pool.getSerialByIndex( self.particle.serial ) ] = dr

    def getDr( self ):
        pool = self.particle.species.pool
        return pool.drs[ pool.getSerialByIndex( self.particle.serial ) ]

    def calculateFirstPassageTime( self ):
        
        species = self.particle.species
        fpgf = _gfrd.FirstPassageGreensFunction( species.D )
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
        self.dt = 0.0
        self.closest = None

        particle1 = self.single1.particle
        particle2 = self.single2.particle

        D12 = particle1.species.D + particle2.species.D
        self.sigma = particle1.species.radius + particle2.species.radius

        self.sgf = _gfrd.FirstPassageGreensFunction( D12 / 4.0 )
        self.pgf = _gfrd.FirstPassagePairGreensFunction( D12, rt.k,
                                                         self.sigma )



    def __del__( self ):
        print 'del', str( self )



    def getCoM( self ):
        particle1 = self.single1.particle
        particle2 = self.single2.particle
        
        D1 = particle1.species.D
        D2 = particle2.species.D

        sqrtD2D1 = math.sqrt( D2 / D1 ) 
        sqrtD1D2 = math.sqrt( D1 / D2 )

        pos1 = particle1.getPos()
        pos2 = particle2.getPos()
        
        #com = sqrtD2D1 * pos1 + sqrtD1D2 * pos2
        com = ( pos1 + sqrtD1D2 * pos2 ) * .5
        return com

    def nextEvent( self, dr ):
        rnd = numpy.random.uniform( size=3 )

        pos1 = self.single1.particle.getPos()
        pos2 = self.single2.particle.getPos()

        r0 = self.sim.distance( pos1, pos2 )

        a_r = ( dr + r0 ) * .5
        a_R = a_r - r0

        print 'dr', dr, 'r0', r0, 'a_r', a_r, 'a_R', a_R, dr - a_r - a_R

        self.pgf.seta( a_r )

        self.t_R = self.sgf.drawTime( rnd[0], a_R )
        self.t_r = self.pgf.drawTime( rnd[1], r0 )

        if self.t_R < self.t_r:
            t = self.t_R
            self.eventType = 2
        else:
            t = self.t_r
            self.eventType = self.pgf.drawEventType( rnd[2], r0, self.t_r )

        return t, self.eventType


    def fire( self ):
        return 0.0


    def update( self, t ):
        print 'update'

    def isDependentOn( self, event ):
        #print event
        return False

    def __str__( self ):
        return str(self.single1.particle) + str(self.single2.particle)


class EGFRDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):
        GFRDSimulatorBase.__init__( self )

        self.isDirty = True

        self.scheduler = _gfrd.EventScheduler()

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

        for single in self.singleMap.values():
            dt = single.calculateFirstPassageTime()
            single.eventID = self.scheduler.addEvent( self.t + dt, single )

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
        self.checkBoundary( pos )
        
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
                single = Single( self, particle )
                single.setDr( 0.0 )
                self.singleMap[ ( species, particle.serial ) ] = single

    def findSingle( self, particle ):
        return self.singleMap.get( ( particle.species, particle.serial ) )

    def initializeSingles( self ):

        for single in self.singleMap.values():
            neighbors, drs = self.getNeighbors( single.particle.getPos() )
            closest = neighbors[1]
            dr = drs[1] - single.particle.species.radius
            closestParticle = Particle( closest[0], index=closest[1] )
            closestSingle = self.findSingle( closestParticle )
            
            single.closest = closestSingle
            #print 'dr', dr
            #FIXME: take different D into account
            single.setDr( dr * .5 )
        
        #FIXME: here perhaps a second pass to get drs maximally large.


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
                pairPartner = neighbors[2]
                pairDr = drs[2]

                if pairDr < r0:   # this happens with a small probability
                    print 'pairDr < r0', pairDr, r0, pairDr - r0
                    #raise ''
                    break
                
                if pairDr > r0 * 5:
                    pairDr = r0 * 5

                nextEvent = pair.nextEvent( pairDr )
                dt = nextEvent[0]
                eventType = nextEvent[1]

                if pair.single1.eventID != None:
                    self.scheduler.removeEvent( pair.single1.eventID )
                if pair.single2.eventID != None:
                    self.scheduler.removeEvent( pair.single2.eventID )

                pair.eventID = self.scheduler.addEvent( self.t + dt, pair ) 
                

    def formPairsGreedily( self ):

        for single in self.singleMap.values():
            print single, closest, dr

    def createPair( self, single1, single2 ):
        rt = self.reactionTypeMap2.get( ( single1.particle.species,\
                                          single2.particle.species ) )
        return Pair( self, single1, single2, rt )

            
        
    def checkShell( self, single ):
        neighbors, drs = self.getNeighborShells( single.particle.getPos() )
        closest, distance = neighbors[1], drs[1]
        distance -= single.particle.species.radius
        if single.getDr() - distance >= 1e-18:
            dr = single.getDr()
            print single.particle, closest, dr, distance, dr - distance
            raise 'Fatal: shells overlap.'

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

