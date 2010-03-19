#!/usr/env python

import weakref

import math

import numpy

from utils import *
from surface import *

from gfrdbase import *
import _gfrd

import logging

import myrandom

log = logging.getLogger('ecell')

DEFAULT_DT_FACTOR = 1e-5

class BDSimulatorCoreBase(object):
    '''
    BDSimulatorCore borrows the following from the main simulator:
    - species_list
    - reaction_types list (both 1 and 2)
    
    '''

    def __init__(self, main):
        self.main = weakref.proxy(main)

        self.particle_list = set()

        self.t = 0.0
        self.dt = 0.0

        self.dt_factor = DEFAULT_DT_FACTOR

        self.step_counter = 0

        self.reaction_events = 0

        self.last_reaction = None

        self.P_acct = {}

    def initialize(self):
        self.determine_dt()

    #@staticmethod  # requires python 2.4 or later.
    def calculate_bd_dt(self, species_list, factor):
        D_list = []
        radius_list = []
        for species in species_list:
            if self.main.particle_pool[species.id]:
                D_list.append(species.D)
                radius_list.append(species.radius)
        D_max = max(D_list) * 2  # max relative diffusion speed
        radius_min = min(radius_list)
        sigma_min = radius_min * 2

        dt = factor * sigma_min ** 2 / D_max  
        if __debug__:
            log.debug('bd dt = %g' % dt)

        return dt


    def clear_particle_list(self):
        self.particle_list = set()

    def add_to_particle_list(self, pid):
        assert type(pid) is ParticleID
        self.particle_list.add(pid)

    def remove_from_particle_list(self, pid):
        self.particle_list.remove(pid)

    def get_next_time(self):
        return self.t + self.dt

    def stop(self, t):
        # dummy
        self.t = t

    def determine_dt(self):
        self.dt = self.calculate_bd_dt(self.main.world.species, self.dt_factor)

    def getP_acct(self, rt, D, sigma):
        try:
            return self.P_acct[rt]
        except KeyError:
            I = _gfrd.I_bd(sigma, self.dt, D)
            p = rt.k * self.dt / (I * 4.0 * numpy.pi)
            if not 0.0 <= p < 1.0:
                raise RuntimeError,\
                    'Invalid acceptance ratio (%s) for reaction %s.' \
                    % (p, rt)
            self.P_acct[rt] = p
            return p

    def step(self):
        self.step_counter += 1
        self.last_reaction = None

        self.propagate()

        self.t += self.dt

    def propagate(self):
        self.particles_to_step = list(self.particle_list)
        myrandom.shuffle(self.particles_to_step)
        while self.particles_to_step:
            pid = self.particles_to_step.pop() # take the last one
            pid_particle_pair = self.main.world.get_particle(pid)
            sid = pid_particle_pair[1].sid

            rt1 = self.attempt_single_reactions(sid)
            if rt1:
                try:
                    self.fire_reaction1(particle, rt1)
                except NoSpace:
                    if __debug__:
                        log.info('fire_reaction1 rejected.')
                continue

            D = pid_particle_pair[1].D
            if D == 0.0:
                continue

            species = self.main.world.get_species(sid)
            surface = self.main.model.get_surface(species.surface_id)
            displacement = surface.draw_bd_displacement(self.dt, D)

            newpos = pid_particle_pair[1].position + displacement
            newpos = self.main.apply_boundary(newpos)

            neighbors = self.get_particles_within_radius(
                newpos, pid_particle_pair[1].radius,
                ignore=[pid_particle_pair[0]])
            if neighbors:

                if len(neighbors) >= 2:
                    if __debug__:
                        log.info('collision two or more particles; move rejected')
                    continue

                closest = neighbors[0][0]
                rt = self.main.get_reaction_rule2(sid, closest[1].sid)[0]

                if rt.k != 0.0:
                    radius12 = pid_particle_pair[1].radius + closest[1].radius
                    D12 = D + closest[1].D

                    p = self.getP_acct(rt, D12, radius12)

                    rnd = myrandom.uniform()

                    if p > rnd:
                        if __debug__:
                            log.info('fire reaction2')
                        try:
                            self.fire_reaction2(pid_particle_pair, closest, rt)
                        except NoSpace:
                            if __debug__:
                                log.info('fire_reaction2 move rejected')
                        continue

                else:
                    if __debug__:
                        log.info('collision move rejected')

                continue

            try:
                self.clear_volume(newpos, pid_particle_pair[1].radius, ignore=[pid_particle_pair[0]])
                self.move_particle((pid_particle_pair[0],
                                   _gfrd.Particle(newpos,
                                            pid_particle_pair[1].radius,
                                            pid_particle_pair[1].D,
                                            pid_particle_pair[1].sid)))
            except NoSpace:
                if __debug__:
                    log.info('propagation move rejected.')

    def attempt_single_reactions(self, sid):
        reaction_types = self.main.get_reaction_rule1(sid)
        if not reaction_types:
            return None  # no reaction

        rnd = myrandom.uniform() / self.dt

        # handle the most common case efficiently.
        if len(reaction_types) == 1:  
            if reaction_types[0].k >= rnd:
                return reaction_types[0]
            else:
                return None

        # if there are more than one possible reaction types..
        k_array = numpy.add.accumulate([rt.k for rt in reaction_types])

        if k_array[-1] < rnd:
            return None

        i = numpy.searchsorted(k_array, rnd)

        return reaction_types[i]

    def fire_reaction1(self, particle, rt):
        oldpos = particle.position

        if len(rt.products) == 0:
            
            self.remove_particle(particle)

            self.last_reaction = Reaction(rt, [particle], [])
            
        elif len(rt.products) == 1:
            
            product_species = rt.products[0]
            radius = product_species.radius

            if self.check_overlap(oldpos, radius,
                                 ignore = [particle, ]):
                if __debug__:
                    log.info('no space for product particle.')
                raise NoSpace()

            self.clear_volume(oldpos, radius, ignore = [particle])
                
            self.remove_particle(particle)
            newparticle = self.create_particle(product_species.id, oldpos)

            self.last_reaction = Reaction(rt, [particle], [newparticle])

            
        elif len(rt.products) == 2:
            
            product_species1 = rt.products[0]
            product_species2 = rt.products[1]
            
            D1 = product_species1.D
            D2 = product_species2.D
            D12 = D1 + D2
            
            radius1 = product_species1.radius
            radius2 = product_species2.radius
            radius12 = radius1 + radius2

            for i in range(self.main.dissociation_retry_moves):

                rnd = myrandom.uniform()
                pair_distance = _gfrd.drawR_gbd(rnd, radius12, self.dt, D12)

                unit_vector = random_unit_vector()
                vector = unit_vector * pair_distance # * (1.0 + 1e-10) # safety
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                newpos1 = oldpos + vector * (D1 / D12)
                newpos2 = oldpos - vector * (D2 / D12)
            
                newpos1 = self.main.apply_boundary(newpos1)
                newpos2 = self.main.apply_boundary(newpos2)

                # accept the new positions if there is enough space.
                if (not self.check_overlap(newpos1, radius1,
                                         ignore = [particle, ])) and \
                   (not self.check_overlap(newpos2, radius2,
                                         ignore = [particle, ])):
                    break
            else:
                if __debug__:
                    log.info('no space for product particles.')
                raise NoSpace()

            self.clear_volume(newpos1, radius1, ignore = [particle])
            self.clear_volume(newpos2, radius2, ignore = [particle])

            # move accepted
            self.remove_particle(particle)

            newparticle1 = self.create_particle(product_species1.id, newpos1)
            newparticle2 = self.create_particle(product_species2.id, newpos2)

            self.last_reaction = Reaction(rt, [particle], 
                                         [newparticle1, newparticle2])

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reaction_events += 1

    def fire_reaction2(self, pid_particle_pair1, pid_particle_pair2, rt):
        if len(rt.products) == 1:
            product_species = rt.products[0]

            D1 = pid_particle_pair1[1].D
            D2 = pid_particle_pair2[1].D

            pos2t = cyclic_transpose(pid_particle_pair1[1].position,
                                     pid_particle_pair2[1].position,
                                     self.main.world.world_size)
            new_pos = (D2 * pid_particle_pair1[1].position + D1 * pos2t) / (D1 + D2)
            new_pos = self.main.apply_boundary(new_pos)

            if self.check_overlap(new_pos, product_species.radius,
                                 ignore=[pid_particle_pair1[0],
                                         pid_particle_pair2[0]]):
                raise NoSpace()
            self.clear_volume(new_pos, product_species.radius,
                             ignore=[pid_particle_pair1[0],
                                     pid_particle_pair2[0]])

            self.remove_particle(pid_particle_pair1)
            self.remove_particle(pid_particle_pair2)
            newparticle = self.create_particle(product_species.id, new_pos)

            try:
                self.particles_to_step.remove(pid_particle_pair2[0])
            except ValueError:  
                pass     # particle2 already stepped, which is fine.

            self.reaction_events += 1

            self.last_reaction = Reaction(
                rt, [pid_particle_pair1, pid_particle_pair2], 
                [newparticle])

            return
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'

    def check(self):
        # particles don't overlap

        for pid, particle in self.main.world:
            assert not self.check_overlap(particle.position, particle.radius,
                                         ignore=[pid, ])


class BDSimulatorCore(BDSimulatorCoreBase):
    def __init__(self, main):
        BDSimulatorCoreBase.__init__(self, main)

    def get_particles_within_radius(self, pos, radius, ignore=[]): 
        return self.main.get_particles_within_radius(pos, radius, ignore)

    def initialize(self):
        BDSimulatorCoreBase.initialize(self)

        self.update_particle_list()

    def update_particle_list(self):
        self.clear_particle_list()
        for s in self.main.world.species:
            for pid in self.main.particle_pool[s.id]:
                self.add_to_particle_list(pid) 

    def add_particle(self, pid_particle_pair):
        self.main.add_particle(pid_particle_pair)
        self.add_to_particle_list(pid_particle_pair[0])

    def move_particle(self, pid_particle_pair):
        return self.main.move_particle(pid_particle_pair)

    def check_overlap(self, pos, radius, ignore=()):
        return self.main.get_particles_within_radius(pos, radius, ignore)

    def remove_particle(self, pid_particle_pair):
        self.main.remove_particle(pid_particle_pair)
        self.remove_from_particle_list(pid_particle_pair[0])

    def create_particle(self, sid, pos):
        particle = self.main.create_particle(sid, pos)
        self.add_to_particle_list(particle[0])

    def clear_volume(self, pos, radius, ignore=[]):
        '''
        This method is a customization point for implementing
        BD in protective domains.
        '''
        pass


class BDSimulator(ParticleSimulatorBase):
    def __init__(self, world):
        ParticleSimulatorBase.__init__(self, world)
        self.core = BDSimulatorCore(self)
        self.is_dirty = True

    def t(self):
        return self.core.t

    def sett(self, t):
        self.core.t = t

    t = property(t, sett)

    def get_dt(self):
        return self.core.dt

    def get_step_counter(self):
        return self.core.step_counter

    dt = property(get_dt)
    step_counter = property(get_step_counter)


    def initialize(self):
        self.core.initialize()
        self.is_dirty = False

    def get_next_time(self):
        return self.core.t + self.core.dt

    def reset(self):
        # DUMMY
        self.core.t=0

    def stop(self, t):
        # dummy
        self.core.stop(t)

    def step(self):
        self.reaction_type = None

        if self.is_dirty:
            self.initialize()

        self.core.step()

        if __debug__:
            log.info('%d: t=%g dt=%g, reactions=%d, rejected_moves=%d' %
                 (self.step_counter, self.t, self.dt, self.reaction_events,
                  self.rejected_moves))

    def check(self):
        pass
