 #!/usr/env python


import math
import sys

import numpy
import scipy


import _gfrd
from utils import *

import os
import logging
import logging.handlers

import myrandom

import model

__all__ = [
    'log',
    'setup_logging',
    'p_free',
    'throw_in_particles',
    'place_particle',
    'NoSpace',
    'create_world',
    'ParticleSimulatorBase',
    ]

World = _gfrd.World

log = logging.getLogger('ecell')

def setup_logging():
    global log 

    if 'LOGFILE' in os.environ:
        if 'LOGSIZE' in os.environ and int(os.environ['LOGSIZE']) != 0:
            handler = logging.handlers.\
                RotatingFileHandler(os.environ['LOGFILE'], mode='w',
                                    maxBytes=int(os.environ['LOGSIZE']))
        else:
            handler = logging.FileHandler(os.environ['LOGFILE'], 'w', )
            
    else:
        handler = logging.StreamHandler(sys.stdout)

    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    if __debug__:
        log.addHandler(handler)
        
    LOGLEVELS = { 'CRITICAL': logging.CRITICAL,
                  'ERROR': logging.ERROR,
                  'WARNING': logging.WARNING,
                  'INFO': logging.INFO, 
                  'DEBUG': logging.DEBUG, 
                  'NOTSET': logging.NOTSET }

    if 'LOGLEVEL' in os.environ:
        if __debug__:
            log.setLevel(LOGLEVELS[os.environ['LOGLEVEL']])
    else:
        if __debug__:
            log.setLevel(logging.INFO)


setup_logging()


def p_free(r, t, D):
    Dt4 = D * t * 4.0
    Pi4Dt = numpy.pi * Dt4
    rsq = r * r
    
    p = math.exp(- rsq / Dt4) / math.sqrt(Pi4Dt * Pi4Dt * Pi4Dt)

    jacobian = 4.0 * numpy.pi * rsq

    return p * jacobian
    
class NoSpace(Exception):
    pass

def get_closest_surface(world, pos, ignore):
    """Return sorted list of tuples with:
        - distance to surface
        - surface itself

    We can not use matrix_space, it would miss a surface if the origin of the 
    surface would not be in the same or neighboring cells as pos.

    """
    structures = []
    distances = []

    for surface in world.structures:
        if isinstance(surface, _gfrd.Surface) and surface.id not in ignore:
            pos_transposed = \
                world.cyclic_transpose(pos, surface.shape.position)
            distance_to_surface = world.distance(surface.shape, pos_transposed)
            distances.append(distance_to_surface)
            structures.append(surface)

    if distances:
        return min(zip(distances, structures))
    else:
        return None

def create_world(m, matrix_size=10):
    world_region = m.get_structure("world")
    if not isinstance(world_region, _gfrd.CuboidalRegion):
        raise TypeError("the world should be a CuboidalRegion")

    if not numpy.all(world_region.shape.extent == world_region.shape.extent[0]):
        raise NotImplementedError("non-cuboidal world is not supported")

    world_size = world_region.shape.extent[0] * 2

    world = _gfrd.World(world_size, matrix_size)

    for st in m.species_types:
        try:
            region = st["surface"]
        except _gfrd.NotFound:
            region = "world"
        world.add_species(
            _gfrd.SpeciesInfo(st.id, 
                              float(st["D"]), 
                              float(st["radius"]), 
                              region))

    for r in m.regions.itervalues():
        world.add_structure(r)

    return world
   
def create_network_rules_wrapper(model):
    return _gfrd.NetworkRulesWrapper(model.network_rules)

def throw_in_particles(world, sid, n):
    species = world.get_species(sid)
    surface = world.get_structure(species.structure_id)

    if __debug__:
        log.info('\tthrowing in %s %s particles to %s' % (n, species.id,
                                                          surface.id))

    # This is a bit messy, but it works.
    i = 0
    while i < int(n):
        position = _gfrd.random_position(surface, myrandom.rng)

        # Check overlap.
        if not world.check_overlap((position, species.radius)):
            create = True
            # Check if not too close to a neighbouring structures for 
            # particles added to the world, or added to a self-defined 
            # box.
            if isinstance(surface, _gfrd.CuboidalRegion):
                tmp = get_closest_surface(world, position, [])
                if tmp is not None:
                    distance, closest_surface = tmp
                    if distance < closest_surface.minimal_distance(species.radius):
                        if __debug__:
                            log.info('\t%d-th particle rejected. Too close to '
                                     'surface. I will keep trying.' % i)
                        create = False
            if create:
                # All checks passed. Create particle.
                p = world.new_particle(sid, position)
                i += 1
                if __debug__:
                    log.info(p)
        elif __debug__:
            log.info('\t%d-th particle rejected. I will keep trying.' % i)

def place_particle(world, sid, pos):
    species = world.get_species(sid)
    radius = species.radius

    if world.check_overlap((pos, radius)):
        raise NoSpace, 'overlap check failed'

    particle = world.new_particle(sid, pos)
    return particle

class ParticleSimulatorBase(object):
    def __init__(self, world, rng, network_rules):
        self.world = world
        self.rng = rng
        self.network_rules = network_rules

        #self.dt = 1e-7
        #self.t = 0.0

        self.H = 3.0

        self.dissociation_retry_moves = 1
        
        self.dt_limit = 1e-3
        self.dt_max = self.dt_limit

        # counters
        self.rejected_moves = 0
        self.reaction_events = 0

        self.max_matrix_size = 0

        self._distance_sq = None
        self._disatnce_sq_array = None

        self.last_reaction = None

    def initialize(self):
        pass

    def get_closest_surface_within_radius(self, pos, radius, ignore):
        """Return sorted list of tuples with:
            - distance to surface
            - surface itself

        """
        distance_to_surface, closest_surface = self.get_closest_surface(pos, 
                                                                   ignore) 
        if distance_to_surface < radius:
            return distance_to_surface, closest_surface
        else:
            return numpy.inf, None

        if isinstance(size, float) and size == INF:
            self._distance = distance
        else:
            self._distance = distance_cyclic

    def get_species(self):
        return self.world.species

    def get_first_pid(self, sid):
        return iter(self.world.get_particle_ids(sid)).next()

    def get_position(self, object):
        if type(object) is tuple and type(object[0]) is _gfrd.ParticleID:
            pid = object[0]
        elif type(object) is _gfrd.ParticleID:
            pid = object
        elif type(object) is _gfrd.Particle:
            pid = self.get_first_pid(object.sid)
        elif type(object) is _gfrd.SpeciesID:
            pid = self.get_first_pid(object)

        return self.world.get_particle(pid)[1].position


    def clear(self):
        self.dt_max = self.dt_limit
        self.dt = self.dt_limit

    def check_particle_matrix(self):
        total = sum(len(self.world.get_particle_ids(s.id)) for s in self.world.species)

        if total != self.world.num_particles:
            raise RuntimeError,\
                'total number of particles %d != self.world.num_particles %d'\
                % (total, self.world.num_particles)

    def check_particles(self):
        for pid_particle_pair in self.world:
            pid = pid_particle_pair[0]
            pos = pid_particle_pair[1].position
            if (pos >= self.world.world_size).any() or (pos < 0.0).any():
                raise RuntimeError,\
                    '%s at position %s out of the world (world size=%g).' %\
                    (pid, pos, self.world.world_size)


    def check(self):
        self.check_particle_matrix()
        self.check_particles()

    def dump_population(self):
        buf = ''
        for species in self.world.species:
            pool = self.world.get_particle_ids(species.id)
            buf += str(sid) + ':' + str(pool) + '\t'

        return buf

    def dump_reaction_rules(self):
        reaction_rules_1 = []
        reaction_rules_2 = []
        reflective_reaction_rules = []
        for si1 in self.world.species:
            for reaction_rule_cache in self.get_reaction_rule1(si1.id):
                string = self.model.dump_reaction_rule(reaction_rule_cache)
                reaction_rules_1.append(string)
            for si2 in self.world.species:
                for reaction_rule_cache in self.get_reaction_rule2(si1.id, si2.id):
                    string = self.model.dump_reaction_rule(reaction_rule_cache)
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


