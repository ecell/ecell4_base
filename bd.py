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
    def __init__(self, tx, network_rules, rng, dt, dissociation_retry_moves, particle_ids):
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
        pid_particle_pair = self.tx.get_particle(pid)
        sid = pid_particle_pair[1].sid

        rt1 = self.attempt_single_reactions(sid)
        if rt1:
            try:
                self.fire_reaction1(pid_particle_pair, rt1)
            except NoSpace:
                if __debug__:
                    log.info('fire_reaction1 rejected.')
            return True

        D = pid_particle_pair[1].D
        if D == 0.0:
            return True

        species = self.tx.get_species(sid)
        surface = self.tx.get_surface(species.surface_id)
        displacement = _gfrd.draw_bd_displacement(surface, math.sqrt(2.0 * D * self.dt), self.rng)
        newpos = self.tx.apply_boundary(pid_particle_pair[1].position + displacement)
        neighbors = self.tx.check_overlap(
            (newpos, pid_particle_pair[1].radius), pid_particle_pair[0])

        if neighbors:
            if len(neighbors) >= 2:
                if __debug__:
                    log.info('collision two or more particles; move rejected')
                return True

            closest = neighbors[0][0]
            reactions = list(self.nr.query_reaction_rule(sid, closest[1].sid))
            if len(reactions) == 0:
                return True
            else:
                assert len(reactions) == 1

            rt = reactions[0]

            if rt.k != 0.0:
                radius12 = pid_particle_pair[1].radius + closest[1].radius
                D12 = D + closest[1].D
                p = self.getP_acct(rt, D12, radius12)
                rnd = self.rng.uniform(0, 1.)

                if p > rnd:
                    try:
                        self.fire_reaction2(pid_particle_pair, closest, rt)
                    except NoSpace:
                        if __debug__:
                            log.info('fire_reaction2 move rejected')
                    return True

            else:
                if __debug__:
                    log.info('collision move rejected')

            return True
    
        try:
            self.tx.update_particle(
                (pid_particle_pair[0],
                 _gfrd.Particle(newpos,
                                pid_particle_pair[1].radius,
                                pid_particle_pair[1].D,
                                pid_particle_pair[1].sid)))
        except NoSpace:
            if __debug__:
                log.info('propagation move rejected.')
        return True

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
        if __debug__:
            log.info('fire_reaction1')
        oldpos = pid_particle_pair[1].position

        if len(rt.products) == 0:
            self.tx.remove_particle(pid_particle_pair[0])
            self.reactions.append((rt, (pid_particle_pair, None), []))
        elif len(rt.products) == 1:
            product_species = rt.products[0]
            radius = product_species.radius

            if self.tx.check_overlap((oldpos, radius),
                                        pid_particle_pair[0]):
                if __debug__:
                    log.info('no space for product particle.')
                raise NoSpace()

            self.tx.remove_particle(pid_particle_pair[0])
            newparticle = self.create_particle(product_species.id, oldpos)

            self.reactions.append((rt, (pid_particle_pair, None), [newparticle]))
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
            
                newpos1 = self.tx.apply_boundary(newpos1)
                newpos2 = self.tx.apply_boundary(newpos2)

                # accept the new positions if there is enough space.
                if (not self.tx.check_overlap(
                        (newpos1, radius1), pid_particle_pair[0])) and \
                   (not self.tx.check_overlap(
                        (newpos2, radius2), pid_particle_pair[0])):
                    break
            else:
                if __debug__:
                    log.info('no space for product particles.')
                raise NoSpace()

            # move accepted
            self.tx.remove_particle(pid_particle_pair[0])

            newparticle1 = self.create_particle(product_species1.id, newpos1)
            newparticle2 = self.create_particle(product_species2.id, newpos2)

            self.reactions.append((rt, (pid_particle_pair, None), 
                                         [newparticle1, newparticle2]))

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

    def fire_reaction2(self, pid_particle_pair1, pid_particle_pair2, rt):
        if __debug__:
            log.info('fire_reaction2')
        if len(rt.products) == 1:
            product_species = self.tx.get_species(rt.products[0])

            D1 = pid_particle_pair1[1].D
            D2 = pid_particle_pair2[1].D

            pos2t = self.tx.cyclic_transpose(pid_particle_pair2[1].position,
                                                pid_particle_pair1[1].position)
            new_pos = self.tx.apply_boundary(
                (D2 * pid_particle_pair1[1].position + D1 * pos2t) / (D1 + D2))

            if self.tx.check_overlap((new_pos, product_species.radius),
                                        pid_particle_pair1[0],
                                        pid_particle_pair2[0]):
                raise NoSpace()

            self.tx.remove_particle(pid_particle_pair1[0])
            self.tx.remove_particle(pid_particle_pair2[0])
            newparticle = self.tx.new_particle(product_species.id, new_pos)

            try:
                self.particles_to_step.remove(pid_particle_pair2[0])
            except ValueError:  
                pass     # particle2 already stepped, which is fine.

            self.reactions.append(
                (rt, (pid_particle_pair1, pid_particle_pair2), 
                [newparticle]))
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'

    def propagate_all(self):
        while self():
            pass

def calculate_bd_dt(species_list):
    D_list = []
    radius_list = []
    D_max = 0.
    radius_min = numpy.inf
    for species in species_list:
        if D_max < species.D:
            D_max = species.D
        if radius_min > species.radius:
            radius_min = species.radius
    return (radius_min * 2) ** 2 / (D_max * 2)

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

    def get_next_time(self):
        return self.t + self.dt

    def stop(self, t):
        # dummy
        self.t = t

    def determine_dt(self):
        self.dt = self.dt_factor * \
               calculate_bd_dt(self.main.world.species)
        if __debug__:
            log.debug('bd dt = %g' % self.dt)

    def step(self):
        self.step_counter += 1

        tx = self.main.world.create_transaction()
        ppg = BDPropagator(tx, self.main.network_rules,
                     myrandom.rng, self.dt, self.main.dissociation_retry_moves,
                     [pid for pid, _ in self.main.world])
        while ppg():
            pass

        self.reaction_events += len(ppg.reactions)
        self.t += self.dt

    def check(self):
        for pp in self.tx:
            assert not self.tx.check_overlap(pp)

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
