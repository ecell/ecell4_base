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

import itertools

log = logging.getLogger('ecell')

DEFAULT_DT_FACTOR = 1e-5

class BDPropagator(object):
    def __init__(self, world, tx, network_rules, rng, dt, dissociation_retry_moves, particle_ids):
        self.world = world
        self.tx = tx
        self.nr = network_rules
        self.rng = rng
        self.dt = dt
        self.dissociation_retry_moves = dissociation_retry_moves
        self.reactions = []
        particles_to_step = list(particle_ids)
        for i in reversed(range(0, len(particles_to_step))):
            j = rng.uniform_int(0, i)
            particles_to_step[i], particles_to_step[j] = \
                particles_to_step[j], particles_to_step[i]
        self.particles_to_step = particles_to_step

    def getP_acct(self, rt, D, sigma):
        I = _gfrd.I_bd(sigma, self.dt, D)
        p = rt.k * self.dt / (I * 4.0 * numpy.pi)
        if not 0.0 <= p < 1.0:
            raise RuntimeError,\
                'Invalid acceptance ratio (%s) for reaction %s.' \
                % (p, rt)
        return p

    def __call__(self):
        if not self.particles_to_step:
            return False

        pid = self.particles_to_step.pop() # take the last one
        pid_particle_pair = self.world.get_particle(pid)
        sid = pid_particle_pair[1].sid

        rt1 = self.attempt_single_reactions(sid)
        if rt1:
            try:
                self.fire_reaction1(pid_particle_pair, rt1)
            except NoSpace:
                if __debug__:
                    log.info('fire_reaction1 rejected.')
            return

        D = pid_particle_pair[1].D
        if D == 0.0:
            return

        species = self.world.get_species(sid)
        surface = self.world.get_surface(species.surface_id)
        displacement = _gfrd.draw_bd_displacement(surface, math.sqrt(2.0 * D * self.dt), myrandom.rng)

        newpos = pid_particle_pair[1].position + displacement
        newpos = self.world.apply_boundary(newpos)

        neighbors = self.world.check_overlap(
            (newpos, pid_particle_pair[1].radius), pid_particle_pair[0])
        if neighbors:

            if len(neighbors) >= 2:
                if __debug__:
                    log.info('collision two or more particles; move rejected')
                return

            closest = neighbors[0][0]
            reactions = list(self.nr.query_reaction_rule(sid, closest[1].sid))
            assert len(reactions) == 1
            rt = reactions[0]

            if rt.k != 0.0:
                radius12 = pid_particle_pair[1].radius + closest[1].radius
                D12 = D + closest[1].D

                p = self.getP_acct(rt, D12, radius12)

                rnd = self.rng.uniform()

                if p > rnd:
                    if __debug__:
                        log.info('fire reaction2')
                    try:
                        self.fire_reaction2(pid_particle_pair, closest, rt)
                    except NoSpace:
                        if __debug__:
                            log.info('fire_reaction2 move rejected')
                    return 

            else:
                if __debug__:
                    log.info('collision move rejected')

            return 
    
        try:
            # self.clear_volume(newpos, pid_particle_pair[1].radius, ignore=[pid_particle_pair[0]])
            self.tx.update_particle(
                (pid_particle_pair[0],
                 _gfrd.Particle(newpos,
                                pid_particle_pair[1].radius,
                                pid_particle_pair[1].D,
                                pid_particle_pair[1].sid)))
        except NoSpace:
            if __debug__:
                log.info('propagation move rejected.')

    def attempt_single_reactions(self, sid):
        reaction_types = self.nr.query_reaction_rule(sid)
        if not reaction_types:
            return None  # no reaction

        rnd = self.rng.uniform() / self.dt

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

    def fire_reaction1(self, pid_particle_pair, rt):
        oldpos = pid_particle_pair[1].position

        if len(rt.products) == 0:
            self.tx.remove_particle(pid_particle_pair[0])
            self.reactions.append(Reaction(rt, [pid_particle_pair], []))
        elif len(rt.products) == 1:
            product_species = rt.products[0]
            radius = product_species.radius

            if self.world.check_overlap((oldpos, radius),
                                        pid_particle_pair[0]):
                if __debug__:
                    log.info('no space for product particle.')
                raise NoSpace()

            # self.clear_volume(oldpos, radius, ignore=[pid_particle_pair[0]])
                
            self.tx.remove_particle(pid_particle_pair[0])
            newparticle = self.create_particle(product_species.id, oldpos)

            self.reactions.append(Reaction(rt, [pid_particle_pair], [newparticle]))
        elif len(rt.products) == 2:
            product_species1 = rt.products[0]
            product_species2 = rt.products[1]
            
            D1 = product_species1.D
            D2 = product_species2.D
            D12 = D1 + D2
            
            radius1 = product_species1.radius
            radius2 = product_species2.radius
            radius12 = radius1 + radius2

            for i in xrange(self.dissociation_retry_moves):
                rnd = self.rng.uniform()
                pair_distance = _gfrd.drawR_gbd(rnd, radius12, self.dt, D12)

                unit_vector = random_unit_vector()
                vector = unit_vector * pair_distance # * (1.0 + 1e-10) # safety
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                newpos1 = oldpos + vector * (D1 / D12)
                newpos2 = oldpos - vector * (D2 / D12)
            
                newpos1 = self.world.apply_boundary(newpos1)
                newpos2 = self.world.apply_boundary(newpos2)

                # accept the new positions if there is enough space.
                if (not self.world.check_overlap(
                        (newpos1, radius1), pid_particle_pair[0])) and \
                   (not self.world.check_overlap(
                        (newpos2, radius2), pid_particle_pair[0])):
                    break
            else:
                if __debug__:
                    log.info('no space for product particles.')
                raise NoSpace()

            #self.clear_volume(newpos1, radius1, ignore=[pid_particle_pair[0]])
            #self.clear_volume(newpos2, radius2, ignore=[pid_particle_pair[0]])

            # move accepted
            self.tx.remove_particle(pid_particle_pair[0])

            newparticle1 = self.create_particle(product_species1.id, newpos1)
            newparticle2 = self.create_particle(product_species2.id, newpos2)

            self.reactions.append(Reaction(rt, [pid_particle_pair], 
                                         [newparticle1, newparticle2]))

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reaction_events += 1

    def fire_reaction2(self, pid_particle_pair1, pid_particle_pair2, rt):
        if len(rt.products) == 1:
            product_species = rt.products[0]

            D1 = pid_particle_pair1[1].D
            D2 = pid_particle_pair2[1].D

            pos2t = self.world.cyclic_transpose(pid_particle_pair1[1].position,
                                                pid_particle_pair2[1].position)
            new_pos = (D2 * pid_particle_pair1[1].position + D1 * pos2t) / (D1 + D2)
            new_pos = self.world.apply_boundary(new_pos)

            if self.world.check_overlap((new_pos, product_species.radius),
                                        pid_particle_pair1[0],
                                        pid_particle_pair2[0]):
                raise NoSpace()
            #self.clear_volume(new_pos, product_species.radius,
            #                 ignore=[pid_particle_pair1[0],
            #                         pid_particle_pair2[0]])

            self.tx.remove_particle(pid_particle_pair1[0])
            self.tx.remove_particle(pid_particle_pair2[0])
            newparticle = self.tx.new_particle(product_species.id, new_pos)

            try:
                self.particles_to_step.remove(pid_particle_pair2[0])
            except ValueError:  
                pass     # particle2 already stepped, which is fine.

            self.reaction_events += 1

            self.reactions.append(Reaction(
                rt, [pid_particle_pair1, pid_particle_pair2], 
                [newparticle]))

            return
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'


class BDSimulatorCoreBase(object):
    '''
    BDSimulatorCore borrows the following from the main simulator:
    - species_list
    - reaction_types list (both 1 and 2)
    
    '''
    def __init__(self, main):
        self.main = weakref.proxy(main)

        self.t = 0.0
        self.dt = 0.0

        self.dt_factor = DEFAULT_DT_FACTOR

        self.step_counter = 0

        self.reaction_events = 0

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

    def get_next_time(self):
        return self.t + self.dt

    def stop(self, t):
        # dummy
        self.t = t

    def determine_dt(self):
        self.dt = self.calculate_bd_dt(self.main.world.species, self.dt_factor)

    def poststep(self):
        tmp = zip(*itertools.chain(tx.modified_particles, tx.added_particles))
        if len(tmp) == 2:
            # tmp[0] -> (pid, pid, pid, ...)
            # tmp[1] -> (particle, particle, particle, ...)
            for p in tmp[1]:
                self.clear_volume(p.position, p.radius, tmp[0])

    def step(self):
        self.step_counter += 1

        tx = self.main.world.create_transaction()
        ppg = BDPropagator(self.world, tx, self.main.network_rules,
                     myrandom.rng, self.dt, self.main.dissociation_retry_moves,
                     [pid for pid, _ in self.world])
        while ppg():
            pass

        self.poststep(tx)

        self.t += self.dt

    def clear_volume(self, pos, radius, ignore=[]):
        '''
        This method is a customization point for implementing
        BD in protective domains.
        '''
        pass

    def check(self):
        # particles don't overlap

        for pid, particle in self.world:
            assert not self.world.check_overlap(
                (particle.position, particle.radius), pid)


class BDSimulatorCore(BDSimulatorCoreBase):
    def __init__(self, main):
        BDSimulatorCoreBase.__init__(self, main)

    def initialize(self):
        BDSimulatorCoreBase.initialize(self)


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
