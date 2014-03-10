#!/usr/bin/env python

import unittest

import numpy

import model
import _gfrd
import gfrdbase
import bd
import logging
import myrandom

bd.log.setLevel(logging.WARNING)


class BDSimulatorTestCase(unittest.TestCase):

    def setUp(self):
        self.m = model.ParticleModel(1e-5)
        self.S = model.Species('S', 2e-11, 5e-8)
        self.A = model.Species('A', 0, 1e-8)
        self.B = model.Species('B', 2e-11, 5e-9)
        self.m.add_species_type(self.S)
        self.m.add_species_type(self.A)
        self.m.add_species_type(self.B)
        self.m.set_all_repulsive()
        self.w = gfrdbase.create_world(self.m, 10)
        self.nrw = _gfrd.NetworkRulesWrapper(self.m.network_rules)
        self.s = bd.BDSimulator(self.w, myrandom.rng, self.nrw)

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        self.failIf(self.s == None)

    
    def test_one_particle(self):
        gfrdbase.place_particle(self.w, self.S, [0.0,0.0,0.0])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_two_particles(self):
        gfrdbase.place_particle(self.w, self.S, [0.0,0.0,0.0])
        gfrdbase.place_particle(self.w, self.S, [5e-6,5e-6,5e-6])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_three_particles(self):
        gfrdbase.place_particle(self.w, self.S, [0.0,0.0,0.0])
        gfrdbase.place_particle(self.w, self.S, [5e-6,5e-6,5e-6])
        gfrdbase.place_particle(self.w, self.S, [1e-7,1e-7,1e-7])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_immobile_is_immobile(self):
        particleA = gfrdbase.place_particle(self.w, self.A, [0.0,0.0,0.0])
        gfrdbase.place_particle(self.w, self.B, [1.5000001e-8,0.0,0.0])

        initial_position = particleA[1].position

        for i in range(10):
            self.s.step()
            #print particleA[1].position
        
        new_position = particleA[1].position
        dist = self.w.distance(initial_position, new_position)

        self.failIf(dist != 0, 'initial pos: %s,\tnew pos: %s' %
                    (initial_position, new_position))



if __name__ == "__main__":
    unittest.main()
