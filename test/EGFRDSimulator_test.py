#!/usr/bin/env python

import logging

import unittest

import numpy

from egfrd import *

log.setLevel(logging.WARNING)


class EGFRDSimulatorTestCase(unittest.TestCase):

    def setUp(self):
        self.m = ParticleModel()
        self.S = self.m.new_species_type('S', 2e-11, 5e-8)
        self.SS = self.m.new_species_type('SS', 1e-12, 5e-9)
        self.A = self.m.new_species_type('A', 0, 1e-8)
        self.B = self.m.new_species_type('B', 2e-11, 5e-9)
        self.m.set_all_repulsive()
        world = World(1e-5, 10) 
        self.s = EGFRDSimulator(world)
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


    def test_three_particles_in_contact(self):
        self.s.place_particle(self.S, [0.0,0.0,0.0])
        self.s.place_particle(self.S, [1e-7,0.0,0.0])

        # dummy
        self.s.place_particle(self.S, [2e-7,0.0,0.0])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_four_particles_close(self):
        self.s.place_particle(self.SS, [2e-8,0.0,0.0])
        self.s.place_particle(self.SS, [3.003e-8,0.0,0.0])

        self.s.place_particle(self.SS, [0.994e-8-5e-10, 0.0, 0.0])

        t = self.s.t
        for i in range(10):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_immobile_is_immobile(self):
        particleA = self.s.place_particle(self.A, [0.0,0.0,0.0])
        self.s.place_particle(self.B, [1.5000001e-8,0.0,0.0])

        initial_position = particleA[1].position

        for i in range(10):
            self.s.step()
        
        new_position = particleA[1].position
        dist = self.s.distance(initial_position, new_position)

        self.failIf(dist != 0, 'initial pos: %s,\tnew pos: %s' %
                    (initial_position, new_position))

    def test_pair_with_immobile(self):
        self.s.place_particle(self.A, [0.0, 0.0, 0.0])
        self.s.place_particle(self.B, [1.51e-8, 0.0, 0.0])

        for i in range(2):
            self.s.step()

    def test_pair_with_immobile_switched_order(self):
        self.s.place_particle(self.B, [1.51e-8, 0.0, 0.0])
        self.s.place_particle(self.A, [0.0, 0.0, 0.0])

        for i in range(2):
            self.s.step()


'''Alternative user interface.

'''
class EGFRDSimulatorTestCaseBase(unittest.TestCase):
    def setUpBase(self):
        self.L = 1e-6

        self.D = 1e-12
        self.radius = 5e-9

        self.A = Species('A', self.D, self.radius)
        self.B = Species('B', self.D, self.radius)
        self.C = Species('C', self.D, self.radius)

        self.kf_1 = 4000
        self.kf_2 = 5e-19
        self.kb_1 = 4000
        self.kb_2 = 4000

        w = World(self.L, 10) 
        self.s = EGFRDSimulator(w)

    def tearDown(self):
        pass


class EGFRDSimulatorCytosoleTestCase(EGFRDSimulatorTestCaseBase):
    def setUp(self):
        EGFRDSimulatorTestCaseBase.setUpBase(self)

        self.m = ParticleModel()

        A = self.A
        B = self.B
        C = self.C

        self.m.add_species(A)
        self.m.add_species(B)
        self.m.add_species(C)

        self.m.add_reaction([A],    [B],    self.kf_1)
        self.m.add_reaction([B],    [A],    self.kb_1)
        self.m.add_reaction([A, B], [C],    self.kf_2)
        self.m.add_reaction([C],    [A, B], self.kb_2)
        self.m.add_reaction([C],    [],     self.kf_1)

        self.s.set_model(self.m)

        self.s.throw_in_particles(A, 2)
        self.s.throw_in_particles(B, 2)

    def test_run(self):
        for i in range(10):
            self.s.step()


class EGFRDSimulatorMembraneTestCase(EGFRDSimulatorTestCaseBase):
    def setUp(self):
        EGFRDSimulatorTestCaseBase.setUpBase(self)

        self.m = ParticleModel()

        A = self.A
        B = self.B
        C = self.C

        L = self.L
        thickness = self.radius

        m1 = self.m.add_planar_surface(origin=[L/2, L/2, 2*L/10],
                                     vector_x=[1, 0, 0],
                                     vector_y=[0, 1, 0],
                                     Lx=L/2,
                                     Ly=L/2,
                                     Lz=thickness,
                                     name='m1')
        self.m.add_species(A, m1)
        self.m.add_species(B, m1)
        self.m.add_species(C, m1)

        self.m.add_reaction([(A, m1)],          [(B, m1)],          self.kf_1)
        self.m.add_reaction([(B, m1)],          [(A, m1)],          self.kb_1)
        self.m.add_reaction([(A, m1), (B, m1)], [(C, m1)],          self.kf_2)
        self.m.add_reaction([(C, m1)],          [(A, m1), (B, m1)], self.kb_2)
        self.m.add_reaction([(C, m1)],          [],                 self.kf_1)

        self.s.set_model(self.m)

        self.s.throw_in_particles(A, 2, m1)
        self.s.throw_in_particles(B, 2, m1)

    def test_run(self):
        for i in range(10):
            self.s.step()


class EGFRDSimulatorDnaTestCase(EGFRDSimulatorTestCaseBase):
    def setUp(self):
        EGFRDSimulatorTestCaseBase.setUpBase(self)

        self.m = ParticleModel()

        A = self.A
        B = self.B
        C = self.C

        L = self.L
        radius = self.radius

        d = self.m.add_cylindrical_surface(origin=[L/2, L/2, L/2],
                                         radius=radius,
                                         orientation=[0, 1, 0],
                                         size=L/2,
                                         name='d')

        self.m.add_species(A, d)
        self.m.add_species(B, d)
        self.m.add_species(C, d)

        self.m.add_reaction([(A, d)],         [(B, d)],         self.kf_1)
        self.m.add_reaction([(B, d)],         [(A, d)],         self.kb_1)
        self.m.add_reaction([(A, d), (B, d)], [(C, d)],         self.kf_2)
        self.m.add_reaction([(C, d)],         [(A, d), (B, d)], self.kb_2)
        self.m.add_reaction([(C, d)],         [],               self.kf_1)

        self.s.set_model(self.m)

        self.s.throw_in_particles(A, 2, d)
        self.s.throw_in_particles(B, 2, d)

    def test_run(self):
        for i in range(10):
            self.s.step()


if __name__ == "__main__":
    unittest.main()
