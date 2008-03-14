#!/usr/env python

import weakref

import math

import numpy
#import scipy
import scipy.optimize


from utils import *
from surface import *

from gfrdbase import *
import _gfrd


class BDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        GFRDSimulatorBase.__init__( self )

        self.isDirty = True

        self.t = 0.0
        self.dt = 0.0

        self.stepCounter = 0

        self.clearPopulationChanged()


    def initialize( self ):

        self.setAllRepulsive()

        self.determineDt()

        self.isDirty = False


    def getNextTime( self ):
        return self.t + self.dt

    def stop( self, t ):
        # dummy
        self.t = t

    def determineDt( self ):

        D_list = []
        radius_list = []
        for species in self.speciesList.values():
            if species.pool.size != 0:
                D_list.append( species.D )
                radius_list.append( species.radius )
        D_max = max( D_list )
        radius_min = min( radius_list )
        sigma_min = radius_min * 2

        DT_FACTOR = 1e-6

        self.dt = DT_FACTOR * sigma_min ** 2 / D_max  
        print 'dt = ', self.dt


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
        
        particleList = []
        for species in self.speciesList.values():
            if species.D != 0.0:
                particleList.extend( [ ( species, s )\
                                           for s in species.pool.serials ] )

        #print particleList

        random.shuffle( particleList )
        for p in particleList:
            self.propagateParticle( p[0], p[1] )


    def propagateParticle( self, species, serial ):

        try:
            i = species.pool.indexMap[ serial ]
        except KeyError:  # already deleted by reaction
            print 'already deleted: ', species.id, serial
            return

        rt1 = self.attemptSingleReactions( species )

        if rt1:
            print rt1
            try:
                self.fireReaction1( Particle( species, serial ), rt1 )
                return
            except NoSpace:
                pass  #FIXME:


        D = species.D
        #if D == 0.0:  # already checked in propagate()
        #    return

        #displacement = drawR_free( self.dt, D )
        displacement = numpy.random.normal( 0.0, 
                                            math.sqrt( 2 * D * self.dt ), 3 )

        newpos = species.pool.positions[i] + displacement
        
        n, d = self.getNeighborParticles( newpos )

        if len( n ) <= 1:  # no other particles; free diffusion.
            species.pool.positions[i] = self.applyBoundary( newpos )
            return

        closest = n[1]
        dist = d[1]
        
        if dist <= species.radius:  # collision
            species2 = closest.species

            rt = self.reactionTypeMap2.get( ( species, species2 ) )
            k = rt.k

            if k != 0.0:
                radius12 = species.radius + species2.radius
                D12 = species.D + species2.D

                print ( radius12, self.dt, D12 )

                I = _gfrd.I_bd( radius12, self.dt, D12 )
                p = k * self.dt / ( I * 4.0 * numpy.pi )
                #print 'p', p
                assert p <= 1.0 and p >= 0.0

                rnd = numpy.random.uniform()
                if p > rnd:
                    print 'fire reaction2'
                    self.fireReaction2( Particle( species, serial ), 
                                        closest, rt )
                    return
                else:
                    print 'reaction reject'
                    return

            else:
                print 'reject'
                return


        species.pool.positions[i] = self.applyBoundary( newpos )
        #print species.pool.positions[i]


    def attemptSingleReactions( self, species ):

        reactionTypes = self.getReactionType1( species )
        if not reactionTypes:
            return None  # no reaction

        k_array = [ rt.k * self.dt for rt in reactionTypes ]
        k_array = numpy.add.accumulate( k_array )
        k_max = k_array[-1]

        rnd = numpy.random.uniform()
        if k_max < rnd:
            return None

        i = numpy.searchsorted( k_array, rnd )

        return reactionTypes[i]


    def fireReaction1( self, particle, rt ):
        
        reactantSpecies = particle.species
        oldpos = particle.pos.copy()

        if len( rt.products ) == 0:
            
            self.removeParticle( particle )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            particle.pos = NOWHERE

            if not self.checkOverlap( oldpos, productSpecies.radius ):
                print 'no space for product particle.'
                particle.pos = oldpos
                raise NoSpace()
                
            self.removeParticle( particle )
            self.placeParticle( productSpecies, oldpos )
            
        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]
            
            D1 = productSpecies1.D
            D2 = productSpecies2.D
            D12 = D1 + D2
            
            particle.pos = NOWHERE

            radius1 = productSpecies1.radius
            radius2 = productSpecies2.radius
            radius12 = radius1 + radius2

            for i in range( 100 ):

                pairDistance = self.drawR_gbd( radius12, self.dt, D12 )
                print pairDistance

                unitVector = randomUnitVector()
                vector = unitVector * pairDistance # * (1.0 + 1e-10) # safety
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                newpos1 = oldpos + vector * ( D1 / D12 )
                newpos2 = oldpos - vector * ( D2 / D12 )
                print oldpos, vector
                print newpos1, newpos2
                #FIXME: check surfaces here
            
                newpos1 = self.applyBoundary( newpos1 )
                newpos2 = self.applyBoundary( newpos2 )

                # accept the new positions if there is enough space.
                if self.checkOverlap( newpos1, radius1 ) and \
                       self.checkOverlap( newpos2, radius2 ):
                    break
            else:
                print 'no space for product particles.'
                particle.pos = oldpos
                raise NoSpace()

            # move accepted
            self.removeParticle( particle )

            particle1 = self.createParticle( productSpecies1, newpos1 )
            particle2 = self.createParticle( productSpecies2, newpos2 )

            print newpos1, newpos2

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1
        self.setPopulationChanged()



    def fireReaction2( self, particle1, particle2, rt ):

        oldPos1 = particle1.pos.copy()
        oldPos2 = particle2.pos.copy()

        particle1.pos = particle2.pos = NOWHERE

        if len( rt.products ) == 1:
                
            productSpecies = rt.products[0]

            D1 = particle1.species.D
            D2 = particle2.species.D

            newPos = ( D2 * oldPos1 + D1 * oldPos2 ) / ( D1 + D2 )

            if not self.checkOverlap( newPos, productSpecies.radius ):
                print 'no space for product particle.'
                particle1.pos = oldPos1
                particle2.pos = oldPos2
                raise NoSpace()
                
            # move accepted
            self.removeParticle( particle1 )
            self.removeParticle( particle2 )
            particle = self.createParticle( productSpecies, newPos )

            self.reactionEvents += 1
            self.setPopulationChanged()
            return
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'



    def drawR_gbd( self, sigma, t, D ):

        def f( r, sigma, t, D, I ):
            return _gfrd.g_bd( r, sigma, t, D ) / I

        I = _gfrd.I_bd( sigma, t, D )

        result = scipy.optimize.brent( f, ( sigma, t, D, I ), 
                                       ( sigma, 
                                         sigma + 6 * math.sqrt( 6 * D * t ) ) )
        r = result

        return r



    def checkInvariants( self ):
        pass
