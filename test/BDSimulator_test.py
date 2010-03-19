#!/usr/bin/env python

import unittest

import numpy

from _gfrd import World
from bd import *

log.setLevel(logging.WARNING)


class BDSimulatorTestCase(unittest.TestCase):

    def setUp(self):
        self.m = ParticleModel()
        self.S = self.m.new_species_type('S', 2e-11, 5e-8)
        self.A = self.m.new_species_type('A', 0, 1e-8)
        self.B = self.m.new_species_type('B', 2e-11, 5e-9)
        self.m.set_all_repulsive()
        world = World(1e-5, 10)
        self.s = BDSimulator(world)
        self.s.set_model(self.m)

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        self.failIf(self.s == None)

    
    def test_one_particle(self):
        self.s.place_particle(self.S, [0.0,0.0,0.0])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_two_particles(self):
        self.s.place_particle(self.S, [0.0,0.0,0.0])
        self.s.place_particle(self.S, [5e-6,5e-6,5e-6])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_three_particles(self):
        self.s.place_particle(self.S, [0.0,0.0,0.0])
        self.s.place_particle(self.S, [5e-6,5e-6,5e-6])
        self.s.place_particle(self.S, [1e-7,1e-7,1e-7])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_immobile_is_immobile(self):
        particleA = self.s.place_particle(self.A, [0.0,0.0,0.0])
        self.s.place_particle(self.B, [1.5000001e-8,0.0,0.0])

        initial_position = particleA[1].position

        for i in range(10):
            self.s.step()
            #print particleA[1].position
        
        new_position = particleA[1].position
        dist = self.s.distance(initial_position, new_position)

        self.failIf(dist != 0, 'initial pos: %s,\tnew pos: %s' %
                    (initial_position, new_position))



if __name__ == "__main__":
    unittest.main()
