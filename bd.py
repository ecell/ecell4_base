#!/usr/env python

import math

import numpy

from utils import *

from gfrdbase import *
import _gfrd

import logging

import itertools

__all__ = [
    'calculate_bd_dt',
    'BDSimulatorCore',
    'BDSimulator',
    ]

log = logging.getLogger('ecell')

DEFAULT_DT_FACTOR = 1e-5

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


class BDSimulatorCore(object):
    '''
    BDSimulatorCore borrows the following from the main simulator:
    - species_list
    - reaction_types list (both 1 and 2)
    
    '''
    def __init__(self, world, rng, network_rules, dissociation_retry_moves):
        self.world = world
        self.rng = rng
        self.network_rules = network_rules
        self.dissociation_retry_moves = dissociation_retry_moves

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
               calculate_bd_dt(self.world.species)
        if __debug__:
            log.debug('bd dt = %g' % self.dt)

    def step(self):
        self.step_counter += 1

        def increment_reaction_events(rr):
            self.reaction_events += 1

        ppg = _gfrd.BDPropagator(self.world, self.network_rules,
                     self.rng, self.dt, self.dissociation_retry_moves,
                     increment_reaction_events, self.world.particle_ids)
        ppg.propagate_all()

        self.t += self.dt

    def check(self):
        for pp in self.tx:
            assert not self.tx.check_overlap(pp)

class BDSimulator(ParticleSimulatorBase):
    def __init__(self, world, rng, network_rules):
        ParticleSimulatorBase.__init__(self, world, rng, network_rules)
        self.is_dirty = True
        self.core = BDSimulatorCore(self.world, self.rng, self.network_rules,
                                    self.dissociation_retry_moves)

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

    def set_dt_factor(self, dt_factor):
        self.core.dt_factor = dt_factor

    def get_dt_factor(self):
        return self.core.dt_factor

    dt_factor = property(get_dt_factor, set_dt_factor)


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
