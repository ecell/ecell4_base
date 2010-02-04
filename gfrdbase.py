 #!/usr/env python


import math
import sys

import numpy
import scipy


#from surface import *
import _gfrd
from utils import *

import os
import logging
import logging.handlers

import myrandom

__all__ = [
    'log',
    'setupLogging',
    'p_free',
    'drawR_free',
    'NoSpace',
    'Species',
    'ReactionRuleCache',
    'ParticleModel',
    'createUnimolecularReactionRule',
    'createDecayReactionRule',
    'createBindingReactionRule',
    'createUnbindingReactionRule',
    'Reaction',
    'ParticleSimulatorBase',
    ]

log = logging.getLogger('ecell')

def setupLogging():
    global log 

    if 'LOGFILE' in os.environ:
        if 'LOGSIZE' in os.environ and int( os.environ[ 'LOGSIZE' ] ) != 0:
            handler = logging.handlers.\
                RotatingFileHandler( os.environ[ 'LOGFILE' ], mode='w',
                                     maxBytes=int( os.environ[ 'LOGSIZE' ] ) )
        else:
            handler = logging.FileHandler( os.environ[ 'LOGFILE' ], 'w', )
            
    else:
        handler = logging.StreamHandler( sys.stdout )

    formatter = logging.Formatter( '%(message)s' )
    handler.setFormatter( formatter )
    if __debug__:
        log.addHandler( handler )
        
    LOGLEVELS = { 'CRITICAL': logging.CRITICAL,
                  'ERROR': logging.ERROR,
                  'WARNING': logging.WARNING,
                  'INFO': logging.INFO, 
                  'DEBUG': logging.DEBUG, 
                  'NOTSET': logging.NOTSET }

    if 'LOGLEVEL' in os.environ:
        if __debug__:
            log.setLevel( LOGLEVELS[ os.environ[ 'LOGLEVEL' ] ] )
    else:
        if __debug__:
            log.setLevel( logging.INFO )


setupLogging()


def p_free( r, t, D ):
    Dt4 = D * t * 4.0
    Pi4Dt = numpy.pi * Dt4
    rsq = r * r
    
    p = math.exp( - rsq / Dt4 ) / math.sqrt( Pi4Dt * Pi4Dt * Pi4Dt )

    jacobian = 4.0 * numpy.pi * rsq

    return p * jacobian
    

def drawR_free( t, D ):
    ro = math.sqrt( 2.0 * D * t )
    return myrandom.normal( 0.0, ro, 3 )


class NoSpace( Exception ):
    pass


class Species( object ):
    def __init__( self, type, D, radius ):
        self.type = type
        self.D = float(D)
        self.radius = float(radius)

    @property
    def serial( self ):
        return self.type.id

    @property
    def id( self ):
        return self.type["id"]

    def __str__(self):
        return str(self.type)


class ReactionRuleCache(object):
    def __init__(self, rt, products, k):
        self.rt = rt
        self.products = products
        self.k = k


class ParticleModel( _gfrd.Model ):
    def new_species_type( self, id, D, radius ):
        st = _gfrd.Model.new_species_type( self )
        st[ "id" ] = str( id )
        st[ "D" ] = str( D )
        st[ "radius" ] = str( radius )
        return st

    def set_all_repulsive( self ):
        nr = self.network_rules
        for species1 in self.species_types:
            for species2 in self.species_types:
                if nr.query_reaction_rule(species1, species2) is None:
                    nr.add_reaction_rule(
                        _gfrd.ReactionRule([species1, species2], [], 0.))


def createUnimolecularReactionRule( s1, p1, k ):
    return _gfrd.ReactionRule( [ s1, ], [ p1, ], k )


def createDecayReactionRule( s1, k ):
    return _gfrd.ReactionRule( [ s1, ], [], k )


def createBindingReactionRule( s1, s2, p1, k ):
    return _gfrd.ReactionRule( [ s1, s2 ], [ p1, ], k )


def createUnbindingReactionRule( s1, p1, p2, k ):
    return _gfrd.ReactionRule( [ s1, ], [ p1, p2 ], k )


class Reaction:
    def __init__( self, type, reactants, products ):
        self.type = type
        self.reactants = reactants
        self.products = products

    def __str__( self ):
        return 'Reaction( ' + str( self.type ) + ', ' + str( self.reactants )\
            + ', ' + str( self.products ) + ' )'


class ParticleSimulatorBase( object ):
    def __init__( self ):
        self.speciesList = {}
        self.particlePool = {}
        self.reactionRuleCache = {}

        self.surfaceList = []

        #self.dt = 1e-7
        #self.t = 0.0

        self.H = 3.0
        
        self.dtLimit = 1e-3
        self.dtMax = self.dtLimit

        # counters
        self.rejectedMoves = 0
        self.reactionEvents = 0

        self.particleMatrix = None

        self.maxMatrixSize = 0

        self.worldSize = 1.0
        self.matrixSize = 10

        self._distanceSq = None
        self._disatnceSqArray = None

        self.lastReaction = None

        self.model = None

        self.particleIDGenerator = _gfrd.ParticleIDGenerator(0)

    def setModel( self, model ):
        model.set_all_repulsive()
        self.speciesList.clear()
        self.particlePool.clear()
        self.reactionRuleCache.clear()
        for st in model.species_types:
            self.speciesList[st.id] = Species( st, st["D"], st["radius"] )
            self.particlePool[st.id] = _gfrd.ParticleIDSet()
        self.model = model

    def initialize( self ):
        pass

    def setWorldSize( self, size ):
        if isinstance( size, list ) or isinstance( size, tuple ):
            size = numpy.array( size )

        self.worldSize = size

        self.particleMatrix = _gfrd.ParticleContainer(
                self.worldSize, self.matrixSize )

        if isinstance( size, float ) and size == INF:
            self._distance = distance_Simple
            #self._distanceArray = distanceArray_Simple
            self._distanceSq = distanceSq_Simple
            self._distanceSqArray = distanceSqArray_Simple
        else:
            self._distance = distance_Cyclic
            #self._distanceArray = distanceSqArray_Cyclic
            self._distanceSq = distanceSq_Cyclic
            self._distanceSqArray = distanceSqArray_Cyclic

    def getWorldSize( self ):
        return self.worldSize

    def setMatrixSize( self, size ):
        if self.maxMatrixSize == 0:
            self.matrixSize = size
        else:
            self.matrixSize = max( size, self.maxMatrixSize )
        self.particleMatrix = _gfrd.ParticleContainer(
                self.worldSize, self.matrixSize )

    def applyBoundary( self, pos ):
        return apply_boundary(pos, self.worldSize)

    def getReactionRule1( self, species ):
        k = (species,)
        try:
            retval = self.reactionRuleCache[k]
        except:
            gen = self.model.network_rules.query_reaction_rule( species )
            if gen is None:
                retval = None
            else:
                retval = []
                for rt in gen:
                    products = [ self.speciesList[st] for st in rt.products ]
                    species1 = self.speciesList[rt.reactants[0]]
                    if len( products ) == 1:
                        if species1.radius * 2 < products[0].radius:
                            raise RuntimeError,\
                                'radius of product must be smaller ' \
                                + 'than radius of reactant.'
                    elif len( products ) == 2:
                        if species1.radius < products[0].radius or\
                                species1.radius < products[1].radius:
                            raise RuntimeError,\
                                'radii of both products must be smaller than ' \
                                + 'reactant.radius.'
                    retval.append(ReactionRuleCache(rt, products, rt.k))
            self.reactionRuleCache[k] = retval
        return retval

    def getReactionRule2( self, species1, species2 ):
        if species2 < species1:
            k = (species2, species1)
        else:
            k = (species1, species2)
        try:
            retval = self.reactionRuleCache[k]
        except:
            gen = self.model.network_rules.query_reaction_rule( species1, species2 )
            if gen is None:
                retval = None
            else:
                retval = []
                for rt in gen:
                    retval.append(
                        ReactionRuleCache(
                            rt,
                            [ self.speciesList[st] for st in rt.products ],
                            rt.k))
            self.reactionRuleCache[k] = retval
        return retval

    def getSpecies(self):
        return self.speciesList.itervalues()

    def getParticlePool(self, sid):
        return self.particlePool[sid]

    def distanceSq( self, position1, position2 ):
        return self._distanceSq( position1, position2, self.worldSize )

    def distance( self, position1, position2 ):
        return self._distance( position1, position2, self.worldSize )
        
    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1,\
                                                 positions ) )

    def addSurface( self, surface ):
        self.surfaceList.append( surface )

    def setAllRepulsive( self ):
        # TODO
        pass
        
    def throwInParticles( self, st, n, surface=[] ):
        species = self.speciesList[st.id]
        if __debug__:
            log.info( 'throwing in %s %s particles' % ( n, species.id ) )

        for i in range( int( n ) ):
            while 1:
                position = surface.randomPosition()
                if not self.checkOverlap( position, species.radius ):
                    break
                else:
                    if __debug__:
                        log.info( '%d-th particle rejected.' %i )
            
            p = self.createParticle(st.id, position)
            if __debug__:
                log.info(p)

    def placeParticle( self, st, pos ):
        species = self.speciesList[st.id]
        pos = numpy.array( pos )
        radius = species.radius

        if self.checkOverlap( pos, radius ):
            raise NoSpace, 'overlap check failed'
            
        particle = self.createParticle( st.id, pos )
        return particle

    def createParticle(self, sid, pos):
        pid = self.particleIDGenerator()
        species = self.speciesList[sid]
        newParticle = _gfrd.Particle(pos, species.radius, species.D, sid)
        return self.addToParticleMatrix(pid, newParticle)

    def removeParticle(self, pid_particle_pair):
        self.particlePool[pid_particle_pair[1].sid].remove(pid_particle_pair[0])
        del self.particleMatrix[pid_particle_pair[0]]

    def moveParticle( self, pid_particle_pair ):
        self.particleMatrix.update(pid_particle_pair)

    def addToParticleMatrix(self, pid, particle):
        assert isinstance(pid, _gfrd.ParticleID)
        assert isinstance(particle, _gfrd.Particle)
        self.particleMatrix[pid] = particle
        self.particlePool[particle.sid].add(pid)
        return (pid, particle)

    def checkOverlap( self, pos, radius, ignore=[] ):
        result = self.particleMatrix.get_neighbors_within_radius( pos, radius )
        for i in result:
            if i[0][0] not in ignore:
                return i
        return None

    def getParticlesWithinRadius(self, pos, radius, ignore=[]):
        result = self.particleMatrix.get_neighbors_within_radius(pos, radius)
        return [ p[0] for p in result if p[0][0] not in ignore ]

    def getParticlesWithinRadiusNoSort(self, pos, radius, ignore=[]):
        return self.getParticlesWithinRadius(pos, radius, ignore)

    def clear( self ):
        self.dtMax = self.dtLimit
        self.dt = self.dtLimit

    def checkParticleMatrix( self ):
        if self.worldSize != self.particleMatrix.world_size:
            raise RuntimeError,\
                'self.worldSize != self.particleMatrix.worldSize'

        total = sum(len(pool) for pool in self.particlePool.itervalues())

        if total != len(self.particleMatrix):
            raise RuntimeError,\
                'total number of particles %d != self.particleMatrix.size %d'\
                % ( total, self.particleMatrix.size )

    def checkParticles(self):
        for i in self.particleMatrix:
            pid = i[0]
            pos = i[1].position
            if (pos >= self.worldSize).any() or (pos < 0.0).any():
                raise RuntimeError,\
                    '%s at position %s out of the world (worldSize=%g).' %\
                    (pid, pos, self.worldSize)


    def check( self ):
        self.checkParticleMatrix()
        self.checkParticles()

    def dumpPopulation( self ):
        buf = ''
        for sid, pool in self.particlePool.iteritems():
            buf += str(sid) + ':' + str(pool) + '\t'

        return buf
