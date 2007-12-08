#!/usr/env python

import weakref

import math

import numpy
#import scipy
#import scipy.optimize


from utils import *
from surface import *

from gfrdbase import *
import _gfrd


class BDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        GFRDSimulatorBase.__init__( self )

        self.isDirty = True

        self.t = 0.0
        self.dt = 1e-9

        self.stepCounter = 0

        self.smallT = 1e-8  # FIXME: is this ok?

        self.lastEvent = None

        self.clearPopulationChanged()


    def initialize( self ):

        self.setAllRepulsive()

        self.isDirty = False


    def step( self ):

        self.clearPopulationChanged()

        if self.isDirty:
            self.initialize()

        if self.stepCounter % 10000 == 0:
            self.checkInvariants()

        self.stepCounter += 1

        self.propagate()

        self.t += self.dt
        print self.stepCounter, ': t = ', self.t, 'dt = ', self.dt, 
        'reactions', self.reactionEvents, 'rejected moves', self.rejectedMoves
        print ''


    def populationChanged( self ):

        return self.isPopulationChanged

    def clearPopulationChanged( self ):

        self.isPopulationChanged = False

    def setPopulationChanged( self ):

        self.isPopulationChanged = True
        

    def propagate( self ):
        
        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                self.propagateParticle( species, i )


    def propagateParticle( self, species, i ):

        D = species.D
        if D == 0.0:
            return

        pos = species.pool.positions[i]

        displacement = drawR_free( self.dt, D )
        newpos = pos + displacement
        
        n, d = self.getNeighborParticles( newpos )

        if len( n ) <= 1:  # no other particles; free diffusion.
            species.pool.positions[i] = self.applyBoundary( newpos )
            return

        closest = n[1]
        dist = d[1]
        
        if dist <= species.radius:
            species2 = closest.species

            rt = self.reactionTypeMap2.get( ( species, species2 ) )
            if rt == None:
                k = 0.0
            else:
                k = rt.k

            if k != 0.0:
                radius12 = species.radius + species2.radius
                D12 = species.D + species2.D

                I = _gfrd.I_bd( radius12, self.dt, D12 )
                print 'I', I
                kdt = k * self.dt
                p = k * self.dt / I
                print 'p', p
                assert p <= 1.0 and p >= 0.0

                rnd = numpy.random.uniform()
                if p < rnd:
                    print 'fire reaction2'
                    self.fireReaction2( Particle( species, i ), closest, rt )
                    return
                else:
                    print 'reaction reject'
                    return

            else:
                print 'reject'
                return


        species.pool.positions[i] = self.applyBoundary( newpos )
        print pos

    def fireReaction2( self, particle1, particle2, rt ):

        if len( rt.products ) == 1:
                
            species3 = rt.products[0]

            D1 = particle1.species.D
            D2 = particle2.species.D

            newPos = ( D2 * particle1.pos + D1 * particle2.pos ) / ( D1 + D2 )

            self.removeParticle( particle1 )
            self.removeParticle( particle2 )

            particle = self.createParticle( species3, newPos )

            self.reactionEvents += 1
            self.setPopulationChanged()
            return
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'




    def checkInvariants( self ):
        pass
