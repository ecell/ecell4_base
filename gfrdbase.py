 #!/usr/env python


import math
import sys

import numpy
import scipy


from surface import *
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
    'NoSpace',
    'ReactionRuleCache',
    'ParticleModel',
    'World',
    'createUnimolecularReactionRule',
    'createDecayReactionRule',
    'createBindingReactionRule',
    'createUnbindingReactionRule',
    'Reaction',
    'ParticleSimulatorBase',
    'Species',
    ]

World = _gfrd.World

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
    

class NoSpace( Exception ):
    pass


class ReactionRuleCache(object):
    def __init__(self, rt, products, k):
        self.rt = rt
        self.products = products
        self.k = k


class DummySpecies( object ):
    """This is needed internally during initialization for the virtual product 
    of a decay or surface absorption reaction.

    For alternative user interface.
    """
    def __init__(self, surface):
        self.radius = 0
        self.surface = surface


class Species(object):
    """For alternative user interface.

    """
    def __init__(self, name, D=None, radius=None, surface=None):
        self.name = name
        self.D = D
        self.radius = radius
        self.surface = surface


class ParticleModel( _gfrd.Model ):
    def new_species_type( self, id, D, radius ):
        st = _gfrd.Model.new_species_type( self )
        st[ "id" ] = str( id )
        st[ "D" ] = str( D )
        st[ "radius" ] = str( radius )
        return st

    def set_all_repulsive( self ):
        nr = self.network_rules
        # Maybe the user has defined a reaction rule for any 2 species since a 
        # previous call to this method, so remove *all* repulsive reaction 
        # rules first.
        for species1 in self.species_types:
            for species2 in self.species_types:
                gen = nr.query_reaction_rule(species1, species2)
                if gen is not None:
                    for reaction_rule in gen:
                        if float(reaction_rule['k']) == 0.0:
                            nr.remove_reaction_rule(reaction_rule)

        for species1 in self.species_types:
            for species2 in self.species_types:
                gen = nr.query_reaction_rule(species1, species2)
                if gen is None or len(set(gen)) == 0:
                    rr = _gfrd.ReactionRule([species1, species2], [])
                    rr['k'] = '0.0'
                    nr.add_reaction_rule(rr)


    '''
    Methods for alternative user interface.
    '''
    def addSpecies(self, species, surface=None, D=None, radius=None):
        """Add a new species.

        A species is a type of particles. By default the species is added to 
        the 'world'. If a surface is specified, it is added to that surface. Per 
        surface a different diffusion constant D and radius can be specified. 
        By default the ones for the 'world' are used.

        species -- a species created with Species().
        surface -- the surface this species can exist on.
        D       -- the diffusion constant of the particles.
        radius  -- the radii of the particles.

        """
        if surface == None:
            name = species.name
            surface = self.defaultSurface
        else:
            assert any(surface == s for s in self.surfaceList.itervalues()), \
                   '%s not in surfaceList.' % (surface)

            # Construct new name if species lives on a surface.
            name = '(' + species.name + ',' + str(surface) + ')'

        if D == None:
            assert species.D != None, \
                   'Diffusion constant of species %s not specified.' % species
            D = species.D

        if radius == None:
            assert species.radius != None, \
                   'Radius of species %s not specified.' % species
            radius = species.radius

        # Create a species type for internal use. Don't use the user defined 
        # species at all. The new name is a concatenation of the user defined 
        # name and the surface this species exists on.
        species_type = self.new_species_type(name, D, radius)
        species_type['surface'] = surface.name

        return species_type

    def get_species_type(self, key):
        """Return species type.

        """
        if isinstance(key, SpeciesType):
            return key
        elif isinstance(key, Species):
            name = key.name
        else:
            # Unpack (species, surface)-key.
            species = key[0]
            surface = key[1]

            if species == 0:
                # This is the virtual product of a decay or surface absorption 
                # reaction.
                species = DummySpecies(surface)

            # Note: see addSpecies for how name is constructed. 
            name = '(' + species.name + ',' + str(surface) + ')'

        return self.get_species_type_by_name(name)

    def get_species_type_by_name(self, name):
        for species_type in self.species_types:
            if species_type['id'] == name:
                return species_type

        raise RuntimeError('Species type %s does not exist.' % (name))

    def addReaction(self, reactants, products, k): 
        reactants = map(self.get_species_type, reactants)
        products  = map(self.get_species_type, products)

        rr = _gfrd.ReactionRule(reactants, products)
        rr['k'] = str(k)
        self.network_rules.add_reaction_rule(rr)
        return rr


def createUnimolecularReactionRule( s1, p1, k ):
    rr = _gfrd.ReactionRule( [ s1, ], [ p1, ] )
    rr['k'] = str(k)
    return rr


def createDecayReactionRule( s1, k ):
    rr = _gfrd.ReactionRule( [ s1, ], [] )
    rr['k'] = str(k)
    return rr


def createBindingReactionRule( s1, s2, p1, k ):
    rr = _gfrd.ReactionRule( [ s1, s2 ], [ p1, ] )
    rr['k'] = str(k)
    return rr


def createUnbindingReactionRule( s1, p1, p2, k ):
    rr = _gfrd.ReactionRule( [ s1, ], [ p1, p2 ] )
    rr['k'] = str(k)
    return rr


class Reaction:
    def __init__( self, type, reactants, products ):
        self.type = type
        self.reactants = reactants
        self.products = products

    def __str__( self ):
        return 'Reaction( ' + str( self.type ) + ', ' + str( self.reactants )\
            + ', ' + str( self.products ) + ' )'


class ParticleSimulatorBase( object ):
    def __init__(self, world):
        self.particlePool = {}
        self.reactionRuleCache = {}

        self.surfaceList = {}

        #self.dt = 1e-7
        #self.t = 0.0

        self.H = 3.0

        self.dissociation_retry_moves = 1
        
        self.dtLimit = 1e-3
        self.dtMax = self.dtLimit

        # counters
        self.rejectedMoves = 0
        self.reactionEvents = 0

        self.world = world

        self.maxMatrixSize = 0

        self._distanceSq = None
        self._disatnceSqArray = None

        self.lastReaction = None

        self.model = None

        self.network_rules = None
    
        # Particles of a Species whose surface is not specified will be added 
        # to the world. Dimensions don't matter, except for visualization.
        self.defaultSurface = \
            CuboidalRegion([0, 0, 0], [self.world.world_size] * 3, 'world')

    def setModel( self, model ):
        model.set_all_repulsive()
        self.reactionRuleCache.clear()

        for st in model.species_types:
            if st["surface"] == "":
                st["surface"] = self.defaultSurface.name

            try:
                self.world.get_species(st.id)
            except:
                self.world.add_species(_gfrd.SpeciesInfo(st.id, 
                                                   float(st["D"]), 
                                                   float(st["radius"]), 
                                                   st["surface"]))
            # else: keep species info for this species.

            if not self.particlePool.has_key(st.id):
                self.particlePool[st.id] = _gfrd.ParticleIDSet()
            # else: keep particle data for this species.

        self.network_rules = _gfrd.NetworkRulesWrapper(model.network_rules)
        self.model = model

    def initialize( self ):
        pass

    def getClosestSurface(self, pos, ignore):
        """Return sorted list of tuples with:
            - distance to surface
            - surface itself

        We can not use matrixSpace, it would miss a surface if the origin of the 
        surface would not be in the same or neighboring cells as pos.

        """
        surfaces = [None]
        distances = [numpy.inf]

        ignoreSurfaces = []
        for obj in ignore:
            if isinstance(obj.surface, Surface):
                # For example ignore surface that particle is currently on.
                ignoreSurfaces.append(obj.surface)

        for surface in self.surfaceList.itervalues():
            if surface not in ignoreSurfaces:
                posTransposed = \
                    self.world.cyclic_transpose(pos, surface.shape.position)
                distanceToSurface = self.distance(surface.shape, posTransposed)
                distances.append(distanceToSurface)
                surfaces.append(surface)

        return min(zip(distances, surfaces))

    def getClosestSurfaceWithinRadius(self, pos, radius, ignore):
        """Return sorted list of tuples with:
            - distance to surface
            - surface itself

        """
        distanceToSurface, closestSurface = self.getClosestSurface(pos, 
                                                                   ignore ) 
        if distanceToSurface < radius:
            return distanceToSurface, closestSurface
        else:
            return numpy.inf, None

        if isinstance( size, float ) and size == INF:
            self._distance = distance
        else:
            self._distance = distance_cyclic

        # Particles of a Species whose surface is not specified will be added 
        # to the world. Dimensions don't matter, except for visualization.

    def applyBoundary( self, pos ):
        return self.world.apply_boundary(pos)

    def getReactionRule1( self, species ):
        k = (species,)
        try:
            retval = self.reactionRuleCache[k]
        except:
            gen = self.network_rules.query_reaction_rule( species )
            if gen is None:
                retval = None
            else:
                retval = []
                for rt in gen:
                    products = [self.world.get_species(sid) for sid in rt.products]
                    species1 = self.world.get_species(rt.reactants[0])
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
            gen = self.network_rules.query_reaction_rule( species1, species2 )
            if gen is None:
                retval = None
            else:
                retval = []
                for rt in gen:
                    retval.append(
                        ReactionRuleCache(
                            rt,
                            [self.world.get_species(sid) for sid in rt.products],
                            rt.k))
            self.reactionRuleCache[k] = retval
        if __debug__:
            if len(retval) > 1:
                name1 = self.model.get_species_type_by_id(species1)['id']
                name2 = self.model.get_species_type_by_id(species2)['id']
                raise RuntimeError('More than 1 bimolecular reaction rule '
                                   'defined for %s + %s: %s' %
                                   (name1, name2, retval))
        return retval

    def getSpecies(self):
        return self.world.species

    def getParticlePool(self, sid):
        return self.particlePool[sid]

    def get_first_pid(self, sid):
        return iter(self.particlePool[sid]).next()

    def get_position(self, object):
        if type(object) is tuple and type(object[0]) is _gfrd.ParticleID:
            pid = object[0]
        elif type(object) is _gfrd.ParticleID:
            pid = object
        elif type(object) is _gfrd.Particle:
            pid = self.get_first_pid(object.sid)
        elif type(object) is _gfrd.SpeciesTypeID:
            pid = self.get_first_pid(object)

        return self.world.get_particle(pid)[1].position

    def distance_between_particles(self, object1, object2):
        pos1 = self.get_position(object1) 
        pos2 = self.get_position(object2)
        return self.distance(pos1, pos2)

    def distance( self, position1, position2 ):
        return self.world.distance( position1, position2 )
        
    def addPlanarSurface(self, name, origin, vectorX, vectorY, Lx, Ly, Lz=0):
        """Add a planar surface.

        name -- a descriptive name, should not be omitted.

        origin -- [x0, y0, z0] is the *center* of the planar surface.
        vectorX -- [x1, y1, z1] and
        vectorY -- [x2, y2, z2] are 2 perpendicular vectors that don't have 
        to be normalized that span the plane. For example [1,0,0] and [0,1,0]
        for a plane at z=0.

        Lx -- lx and 
        Ly -- ly are the distances from the origin of the plane along vectorX 
            or vectorY *to an edge* of the plane. PlanarSurfaces are finite.
        Lz -- dz, the thickness of the planar surface, can be omitted for Lz=0.

        """
        return self.addSurface(PlanarSurface(name, origin, vectorX, vectorY,
                                             Lx, Ly, Lz))

    def addCylindricalSurface(self, name, origin, radius, orientation, size):
        """Add a cylindrical surface.

        name -- a descriptive name, should not be omitted.

        origin -- [x0, y0, z0] is the *center* of the cylindrical surface.
        radius -- r is the radis of the cylinder.
        orientation -- [x1, y1, z1] is a vector that doesn't have to
            normalized that defines the orienation of the cylinder. For 
            example [0,0,1] for a for a cylinder along the z-axis.
        size -- lz is the distances from the origin of the cylinder along 
            the oriention vector to the end of the cylinder. So effectively
            the *half-length*. CylindricalSurfaces are finite.

        """
        return self.addSurface(CylindricalSurface(name, origin, radius, 
                                                  orientation, size))

    def addSurface( self, surface ):
        if(not isinstance(surface, Surface) or
           isinstance(surface, CuboidalRegion)):
            raise RuntimeError(str(surface) + ' is not a surface.')

        self.surfaceList[surface.name] = surface
        return surface

    def getSurface(self, species): 
        nameOfSurface = species.surface
        if nameOfSurface != self.defaultSurface.name:
            surface = self.surfaceList[nameOfSurface]
        else:
            # Default surface is not stored in surfaceList, because it's not 
            # really a surface, and should not be found in getClosestSurface 
            # searches.
            surface = self.defaultSurface
        return surface

    def setAllRepulsive( self ):
        # TODO
        pass
        
    def throwInParticles(self, st, n, surface=None):
        if surface == None or isinstance(surface, CuboidalRegion):
            # Look up species_type.
            # For alternative user interface.
            st = self.model.get_species_type(st)
        else:
            # Look up species_type, given a species and a surface.
            # For alternative user interface.
            st = self.model.get_species_type((st, surface))
        species = self.world.get_species(st.id)
        if surface == None:
            surface = self.getSurface(species)

        if __debug__:
            log.info('\tthrowing in %s %s particles to %s' % (n, species.id,
                                                              surface))

        # This is a bit messy, but it works.
        i = 0
        while i < int(n):
            position = surface.randomPosition()

            # Check overlap.
            if not self.getParticlesWithinRadius(position, species.radius):
                create = True
                # Check if not too close to a neighbouring surfaces for 
                # particles added to the world, or added to a self-defined 
                # box.
                if surface == self.defaultSurface or \
                   (surface != self.defaultSurface and 
                    isinstance( surface, CuboidalRegion)):
                    distance, closestSurface = self.getClosestSurface(position,
                                                                      [])
                    if (closestSurface and
                        distance < closestSurface.minimalDistanceFromSurface( 
                                    species.radius)):
                        if __debug__:
                            log.info('\t%d-th particle rejected. Too close to '
                                     'surface. I will keep trying.' % i)
                        create = False
                if create:
                    # All checks passed. Create particle.
                    p = self.createParticle(st.id, position)
                    i += 1
                    if __debug__:
                        log.info(p)
            elif __debug__:
                log.info('\t%d-th particle rejected. I will keep trying.' % i)

    def placeParticle( self, st, pos ):
        species = self.world.get_species(st.id)
        pos = numpy.array( pos )
        radius = species.radius

        if self.getParticlesWithinRadius( pos, radius ):
            raise NoSpace, 'overlap check failed'
            
        particle = self.createParticle( st.id, pos )
        return particle

    def createParticle(self, sid, pos):
        pid_particle_pair = self.world.new_particle(sid, pos)
        self.particlePool[sid].add(pid_particle_pair[0])
        return pid_particle_pair

    def removeParticle(self, pid_particle_pair):
        self.particlePool[pid_particle_pair[1].sid].remove(pid_particle_pair[0])
        self.world.remove_particle(pid_particle_pair[0])

    def moveParticle( self, pid_particle_pair ):
        self.world.update_particle(pid_particle_pair)

    def getParticlesWithinRadius(self, pos, radius, ignore=[]):
        return self.world.check_overlap((pos, radius), ignore)

    def clear( self ):
        self.dtMax = self.dtLimit
        self.dt = self.dtLimit

    def checkParticleMatrix( self ):
        total = sum(len(pool) for pool in self.particlePool.itervalues())

        if total != self.world.num_particles:
            raise RuntimeError,\
                'total number of particles %d != self.world.num_particles %d'\
                % (total, self.world.num_particles)

    def checkParticles(self):
        for pid_particle_pair in self.world:
            pid = pid_particle_pair[0]
            pos = pid_particle_pair[1].position
            if (pos >= self.world.world_size).any() or (pos < 0.0).any():
                raise RuntimeError,\
                    '%s at position %s out of the world (world size=%g).' %\
                    (pid, pos, self.world.world_size)


    def check( self ):
        self.checkParticleMatrix()
        self.checkParticles()

    def dumpPopulation( self ):
        buf = ''
        for sid, pool in self.particlePool.iteritems():
            buf += str(sid) + ':' + str(pool) + '\t'

        return buf

    def dump_reaction_rule(self, reaction_rule):
        '''Pretty print reaction rule.

        ReactionRule.__str__ would be good, but we are actually getting a 
        ReactionRuleInfo or ReactionRuleCache object.

        This method needs access to self.model.

        '''
        buf = ('k=%.3g' % reaction_rule.k + ': ').ljust(15)
        for index, sid in enumerate(reaction_rule.rt.reactants):
            if index != 0:
                buf += ' + '
            reactant = self.model.get_species_type_by_id(sid)
            buf += reactant['id'].ljust(15)
        if len(reaction_rule.products) == 0:
            if reaction_rule.k != 0:
                buf += '..decays'
        else:
            buf += '-> '

        for index, sid in enumerate(reaction_rule.products):
            if index != 0:
                buf += ' + '
            product = self.model.get_species_type_by_id(sid)
            buf += product['id'].ljust(15)

        return buf + '\n'

    def dump_reaction_rules(self):
        reaction_rules_1 = []
        reaction_rules_2 = []
        reflective_reaction_rules = []
        for si1 in self.world.species:
            for reaction_rule_cache in self.getReactionRule1(si1.id):
                string = self.dump_reaction_rule(reaction_rule_cache)
                reaction_rules_1.append(string)
            for si2 in self.world.species:
                for reaction_rule_cache in self.getReactionRule2(si1.id, si2.id):
                    string = self.dump_reaction_rule(reaction_rule_cache)
                    if reaction_rule_cache.k > 0:
                        reaction_rules_2.append(string)
                    else:
                        reflective_reaction_rules.append(string)

        reaction_rules_1 = uniq(reaction_rules_1)
        reaction_rules_1.sort()
        reaction_rules_2 = uniq(reaction_rules_2)
        reaction_rules_2.sort()
        reflective_reaction_rules = uniq(reflective_reaction_rules)
        reflective_reaction_rules.sort()

        return('\nMonomolecular reaction rules:\n' + ''.join(reaction_rules_1) +
               '\nBimolecular reaction rules:\n' + ''.join(reaction_rules_2) +
               '\nReflective bimolecular reaction rules:\n' +
               ''.join(reflective_reaction_rules))

