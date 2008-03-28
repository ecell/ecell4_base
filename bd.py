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

DEFAULT_DT_FACTOR = 1e-4

def drawR_gbd( sigma, t, D ):
    
    def f( r, sigma, t, D, I ):
        return _gfrd.g_bd( r, sigma, t, D ) / I

    I = _gfrd.I_bd( sigma, t, D )
    
    result = scipy.optimize.brent( f, ( sigma, t, D, I ), 
                                   ( sigma, 
                                     sigma + 6 * math.sqrt( 6 * D * t ) ) )
    return result

def calculateBDDt( speciesList, factor = DEFAULT_DT_FACTOR ):

        D_list = []
        radius_list = []
        for species in speciesList:
            if species.pool.size != 0:
                D_list.append( species.D )
                radius_list.append( species.radius )
        D_max = max( D_list )
        radius_min = min( radius_list )
        sigma_min = radius_min * 2

        dt = factor * sigma_min ** 2 / D_max  
        print 'dt = ', dt

        return dt

class BDSimulatorCoreBase( object ):
    
    '''
    BDSimulatorCore borrows the following from the main simulator:
    - speciesList
    - reactionTypes list (both 1 and 2)
    
    '''

    def __init__( self, main, particleMatrix=None ):

        self.main = weakref.proxy( main )

        self.particleList = []

        self.t = 0.0
        self.dt = 0.0

        self.stepCounter = 0

        self.reactionEvents = 0

        self.speciesList = self.main.speciesList
        self.getReactionType1 = self.main.getReactionType1
        self.getReactionType2 = self.main.getReactionType2
        self.applyBoundary = self.main.applyBoundary


        
    def initialize( self ):
        self.determineDt()

    def clearParticleList( self ):
        self.particleList = []

    def addToParticleList( self, particle ):
        self.particleList.append( particle )

    def removeFromParticleList( self, particle ):
        self.particleList.remove( particle )


    def getNextTime( self ):
        return self.t + self.dt

    def stop( self, t ):
        # dummy
        self.t = t

    def determineDt( self ):

        self.dt = calculateBDDt( self.speciesList.values() )


    def step( self ):

        self.stepCounter += 1

        self.propagate()

        self.t += self.dt


    def propagate( self ):
        
        particles = self.particleList[:]

        #print particleList

        random.shuffle( particles )
        for p in particles:
            self.propagateParticle( p )


    def propagateParticle( self, particle ):

        species = particle.species
        rt1 = self.attemptSingleReactions( species )

        if rt1:
            print rt1
            try:
                self.fireReaction1( particle, rt1 )
                return
            except NoSpace:
                pass  #FIXME:


        D = species.D
        if D == 0.0:
            return

        #displacement = drawR_free( self.dt, D )
        displacement = numpy.random.normal( 0.0, 
                                            math.sqrt( 2 * D * self.dt ), 3 )

        newpos = particle.pos + displacement
        
        closest, dist =\
            self.getClosestParticle( newpos, 
                                     ignore = [ particle, ] )
        
        if closest:  # no other particles; free diffusion.
            newpos = self.applyBoundary( newpos )
            self.moveParticle( particle, newpos )
            return

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


        newpos = self.applyBoundary( newpos )
        self.moveParticle( particle, newpos )



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

            if not self.checkOverlap( oldpos, productSpecies.radius,
                                           ignore = [ particle, ] ):
                print 'no space for product particle.'
                raise NoSpace()
                
            self.removeParticle( particle )
            self.placeParticle( productSpecies, oldpos )
            
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

                pairDistance = drawR_gbd( radius12, self.dt, D12 )
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
                if self.checkOverlap( newpos1, radius1,
                                           ignore = [ particle, ]) and \
                       self.checkOverlap( newpos2, radius2,
                                               ignore = [ particle, ]):
                    break
            else:
                print 'no space for product particles.'
                raise NoSpace()

            # move accepted
            self.removeParticle( particle )

            self.createParticle( productSpecies1, newpos1 )
            self.createParticle( productSpecies2, newpos2 )

            print newpos1, newpos2

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1
        self.populationChanged = True



    def fireReaction2( self, particle1, particle2, rt ):

        pos1 = particle1.pos
        pos2 = particle2.pos

        if len( rt.products ) == 1:
                
            productSpecies = rt.products[0]

            D1 = particle1.species.D
            D2 = particle2.species.D

            newPos = ( D2 * pos1 + D1 * pos2 ) / ( D1 + D2 )

            if not self.checkOverlap( newPos, productSpecies.radius,
                                      ignore=[ particle1, particle2 ]):
                print 'no space for product particle.'
                raise NoSpace()
                
            # move accepted
            self.removeParticle( particle1 )
            self.removeParticle( particle2 )
            self.createParticle( productSpecies, newPos )

            self.reactionEvents += 1
            self.populationChanged = True
            return
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'


    def check( self ):
        pass


'''
IndependentBDSimulatorCore holds its own ObjectMatrix to
calculate distances between particles.   
'''

class IndependentBDSimulatorCore( BDSimulatorCoreBase ):
    
    def __init__( self, main ):

        BDSimulatorCoreBase.__init__( self, main )

        self.checkOverlap = self.main.checkOverlap

        self.particleMatrix = SimpleObjectMatrix()
        self.particleMatrix.setWorldSize( self.main.worldSize )

        self.getNeighborParticles = self._getNeighborParticles
        self.getClosestParticle = self._getClosestParticle

        

    def initialize( self ):
        BDSimulatorCoreBase.initialize( self )

    def addParticle( self, particle ):
        self._addParticle( particle )

    def _addParticle( self, particle ):
        self.addToParticleList( particle )
        self.particleMatrix.add( ( particle[0], particle[1] ),
                                 particle.pos, particle.radius )

    def _updateParticle( self, particle, pos ):
        self.particleMatrix.update( ( particle[0], particle[1] ), pos,
                                    particle.radius )

    def removeParticle( self, particle ):
        self.main.removeParticle( particle )
        self.removeFromParticleList( particle )
        self.particleMatrix.remove( ( particle[0], particle[1] ) )

    def createParticle( self, species, pos ):
        particle = self.main.createParticle( species, pos )
        self._addParticle( particle )

    def placeaParticle( self, species, pos ):
        self.main.placeParticle( species, pos )
        self._addParticle( Particle( species, pos ) )

    def moveParticle( self, particle, pos ):
        self.main.moveParticle( particle, pos )
        self._updateParticle( particle, pos )


    def _getNeighborParticles( self, pos, n=None, dummy=None ):

        n, d = self.particleMatrix.getNeighbors( pos, n, dummy )
        neighbors = [ Particle( i[0], i[1] ) for i in n ]
        return neighbors, d


    def _getClosestParticle( self, pos, ignore=[] ):

        neighbors, distances =\
            self.getNeighborParticles( pos, 
                                       len( ignore ) + 1,
                                       dummy = ( None, -1 ) )

        for i in range( len( neighbors ) ): 
            if neighbors[i] not in ignore:
                closest, distance = neighbors[i], distances[i]

                #assert not closest in ignore
                return closest, distance

        # default case: none left.
        return None, numpy.inf
        #return DummyParticle(), numpy.inf




    def check( self ):
        pass


'''
DependentBDSimulatorCore uses its main simulator's ObjectMatrix
to hold and calculate distances between particles.
'''

class DependentBDSimulatorCore( BDSimulatorCoreBase ):
    

    def __init__( self, main ):

        BDSimulatorCoreBase.__init__( self, main )

        self.checkOverlap = self.main.checkOverlap

        self.moveParticle = self.main.moveParticle

        self.getNeighborParticles = main.getNeighborParticles
        self.getClosestParticle = main.getClosestParticle

        
    def initialize( self ):
        BDSimulatorCoreBase.initialize( self )

        self.updateParticleList()

    def updateParticleList( self ):

        self.clearParticleList()

        for s in self.speciesList.values():
            for i in range( s.pool.size ):
                self.addToParticleList( Particle( s, index = i ) )

    def removeParticle( self, particle ):
        self.main.removeParticle( particle )
        self.removeFromParticleList( particle )

    def createParticle( self, species, pos ):
        particle = self.main.createParticle( species, pos )
        self.addToParticleList( particle )

    def placeaParticle( self, species, pos ):
        self.main.placeParticle( species, pos )
        self.addToParticleList( particle )





class BDSimulator( GFRDSimulatorBase ):
    
    def __init__( self ):

        GFRDSimulatorBase.__init__( self )

        self.core = DependentBDSimulatorCore( self )
        self.isDirty = True

    def gett( self ):
        return self.core.t

    def sett( self, t ):
        self.core.t = t

    def getDt( self ):
        return self.core.dt

    def getStepCounter( self ):
        return self.core.stepCounter

    t = property( gett, sett )
    dt = property( getDt )
    stepCounter = property( getStepCounter )

    def initialize( self ):

        self.setAllRepulsive()

        self.core.initialize()

        self.isDirty = False



    def getNextTime( self ):
        return self.core.t + self.core.dt

    def stop( self, t ):
        # dummy
        self.core.stop( t )

    def step( self ):

        self.populationChanged = False

        if self.isDirty:
            self.initialize()

        if self.stepCounter % 10000 == 0:
            self.check()

        self.core.step()

        print self.stepCounter, ': t = ', self.t,\
            'dt = ', self.dt,\
            'reactions:', self.reactionEvents,\
            'rejected moves:', self.rejectedMoves
        #print ''



    def check( self ):
        pass
