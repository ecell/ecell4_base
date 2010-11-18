#!/usr/bin/env python

import logging

import unittest

import numpy

import _gfrd

from egfrd import *
from visualization import vtklogger

import model
import gfrdbase
import myrandom

log.setLevel(logging.WARNING)


class EGFRDSimulatorTestCase(unittest.TestCase):

    def setUp(self):
        self.m = model.ParticleModel(1e-5)
        self.S = model.Species('S', 2e-11, 5e-8)
        self.SS = model.Species('SS', 1e-12, 5e-9)
        self.A = model.Species('A', 0, 1e-8)
        self.B = model.Species('B', 2e-11, 5e-9)
        self.m.add_species_type(self.S)
        self.m.add_species_type(self.SS)
        self.m.add_species_type(self.A)
        self.m.add_species_type(self.B)
        self.m.set_all_repulsive()
        self.w = gfrdbase.create_world(self.m)
        self.nrw = _gfrd.NetworkRulesWrapper(self.m.network_rules)
        self.s = EGFRDSimulator(self.w, myrandom.rng, self.nrw)

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        self.failIf(self.s == None)


    def test_one_particle(self):
        place_particle(self.s.world, self.S, [0.0,0.0,0.0])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_two_particles(self):
        place_particle(self.s.world, self.S, [0.0,0.0,0.0])
        place_particle(self.s.world, self.S, [5e-6,5e-6,5e-6])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_three_particles(self):
        place_particle(self.s.world, self.S, [0.0,0.0,0.0])
        place_particle(self.s.world, self.S, [5e-6,5e-6,5e-6])
        place_particle(self.s.world, self.S, [1e-7,1e-7,1e-7])

        t = self.s.t
        for i in range(5):
            self.s.step()
        self.failIf(t == self.s.t)


    def test_three_particles_in_contact(self):
        place_particle(self.s.world, self.S, [0.0,0.0,0.0])
        place_particle(self.s.world, self.S, [1e-7,0.0,0.0])

        # dummy
        place_particle(self.s.world, self.S, [2e-7,0.0,0.0])

        r = model.create_unimolecular_reaction_rule(self.S, self.A, 1e10)
        self.m.network_rules.add_reaction_rule(r)

        t = self.s.t
        for i in range(5):
            self.s.step()
            
            # Check if species ids are consistent after unimolecular 
            # multi reaction.
            for species in self.s.world.species:
                for pid in self.s.world.get_particle_ids(species.id):
                    particle = self.s.world.get_particle(pid)[1]
                    self.failIf(particle.sid != species.id)
        self.failIf(t == self.s.t)

    def test_four_particles_close(self):
        place_particle(self.s.world, self.SS, [2e-8,0.0,0.0])
        place_particle(self.s.world, self.SS, [3.003e-8,0.0,0.0])

        place_particle(self.s.world, self.SS, [0.994e-8-5e-10, 0.0, 0.0])

        t = self.s.t
        for i in range(10):
            self.s.step()
        self.failIf(t == self.s.t)

    def test_immobile_is_immobile(self):
        particleA = place_particle(self.s.world, self.A, [0.0,0.0,0.0])
        place_particle(self.s.world, self.B, [1.5000001e-8,0.0,0.0])

        initial_position = particleA[1].position

        for i in range(10):
            self.s.step()
        
        new_position = particleA[1].position
        dist = self.w.distance(initial_position, new_position)

        self.failIf(dist != 0, 'initial pos: %s,\tnew pos: %s' %
                    (initial_position, new_position))

    def test_pair_with_immobile(self):
        place_particle(self.s.world, self.A, [0.0, 0.0, 0.0])
        place_particle(self.s.world, self.B, [1.51e-8, 0.0, 0.0])

        for i in range(2):
            self.s.step()

    def test_pair_with_immobile_switched_order(self):
        place_particle(self.s.world, self.B, [1.51e-8, 0.0, 0.0])
        place_particle(self.s.world, self.A, [0.0, 0.0, 0.0])

        for i in range(2):
            self.s.step()


class EGFRDSimulatorTestCaseBase(unittest.TestCase):
    """Base class for TestCases below.

    """
    def create_model(self):
        self.L = 1e-6

        self.D = 1e-12
        self.radius = 5e-9

        self.m = model.ParticleModel(self.L)

        self.A = model.Species('A', self.D, self.radius)
        self.B = model.Species('B', self.D, self.radius)
        self.C = model.Species('C', self.D, self.radius)

        self.kf_1 = 4000
        self.kf_2 = 5e-19
        self.kb_1 = 4000
        self.kb_2 = 4000

    def add_planar_surface(self):
        m1 = model.create_planar_surface('m1',
                                         [0, 0, 0],
                                         [1, 0, 0],
                                         [0, 1, 0],
                                         self.L,
                                         self.L)

        self.m.add_structure(m1)

    def add_cylindrical_surface(self):
        radius = self.radius
        d = model.create_cylindrical_surface('d',
                                             [0, 0, 0],
                                             radius,
                                             [0, 1, 0],
                                             self.L)

        self.m.add_structure(d)

    def add_species(self):
        self.m.add_species_type(self.A)
        self.m.add_species_type(self.B)
        self.m.add_species_type(self.C)

    def create_simulator(self):
        self.w = create_world(self.m)
        self.s = EGFRDSimulator(self.w)

    def add_reactions(self):
        A = self.A
        B = self.B
        C = self.C

        r = model.create_unimolecular_reaction_rule(A, B, self.kf_1)
        self.m.network_rules.add_reaction_rule(r)
        r = model.create_unimolecular_reaction_rule(B, A, self.kb_1)
        self.m.network_rules.add_reaction_rule(r)
        r = model.create_binding_reaction_rule(A, B, C, self.kf_2)
        self.m.network_rules.add_reaction_rule(r)
        r = model.create_unbinding_reaction_rule(C, A, B, self.kb_2)
        self.m.network_rules.add_reaction_rule(r)
        r = model.create_decay_reaction_rule(C, self.kb_1)
        self.m.network_rules.add_reaction_rule(r)

    def add_particles(self, n):
        throw_in_particles(self.w, self.A, n)
        throw_in_particles(self.w, self.B, n)

    def tearDown(self):
        pass


class CytosoleTestCase(EGFRDSimulatorTestCaseBase):
    """Events happening in the "world".

    """
    def setUp(self):
        self.create_model()
        self.add_species() 
        self.create_simulator() 
        self.add_reactions()
        self.add_particles(2)

    def test_run(self):
        for i in range(10):
            self.s.step()

    def test_vtklogger(self):
        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')
        for i in range(10):
            vtk_logger.log()
            self.s.step()
        vtk_logger.stop()
        vtk_logger.cleanup()

class PlanarSurfaceTestCase(EGFRDSimulatorTestCaseBase):
    """Events happening *on* a planar surface.

    """
    def setUp(self):
        self.create_model()
        self.add_planar_surface()

        # All species on planar surface.
        self.A["structure"] = "m1"
        self.B["structure"] = "m1"
        self.C["structure"] = "m1"

        self.add_species() 
        self.create_simulator() 
        self.add_reactions()
        self.add_particles(2)

    def test_run(self):
        for i in range(10):
            self.s.step()

    def test_vtklogger(self):
        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')
        for i in range(10):
            vtk_logger.log()
            self.s.step()
        vtk_logger.stop()
        vtk_logger.cleanup()


class CylindricalSurfaceTestCase(EGFRDSimulatorTestCaseBase):
    """Events happening *on* a cylindrical surface.

    """
    def setUp(self):
        self.create_model()
        self.add_cylindrical_surface()

        # All species on cylindrical surface.
        self.A["structure"] = "d"
        self.B["structure"] = "d"
        self.C["structure"] = "d"

        self.add_species() 
        self.create_simulator() 
        self.add_reactions()
        self.add_particles(2)

    def test_run(self):
        for i in range(10):
            self.s.step()

    def test_vtklogger(self):
        vtk_logger = vtklogger.VTKLogger(self.s, 'vtk_temp_data')
        for i in range(10):
            vtk_logger.log()
            self.s.step()
        vtk_logger.stop()
        vtk_logger.cleanup()


if __name__ == "__main__":
    unittest.main()
