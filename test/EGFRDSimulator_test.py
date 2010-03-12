#!/usr/bin/env python

import logging

import unittest

import numpy

from egfrd import *

log.setLevel( logging.WARNING )


class EGFRDSimulatorTestCase( unittest.TestCase ):

    def setUp(self):
        self.m = ParticleModel()
        self.S = self.m.new_species_type( 'S', 2e-11, 5e-8 )
        self.SS = self.m.new_species_type( 'SS', 1e-12, 5e-9 )
        self.A = self.m.new_species_type( 'A', 0, 1e-8 )
        self.B = self.m.new_species_type( 'B', 2e-11, 5e-9 )
        self.m.set_all_repulsive()
        world = World( 1e-5, 10 ) 
        self.s = EGFRDSimulator( world )
        self.s.setModel( self.m )

    def tearDown(self):
        pass
    
    def test_instantiation( self ):
        self.failIf( self.s == None )


    def test_OneParticle( self ):
        self.s.placeParticle( self.S, [0.0,0.0,0.0] )

        t = self.s.t
        for i in range( 5 ):
            self.s.step()
        self.failIf( t == self.s.t )

    def test_TwoParticles( self ):
        self.s.placeParticle( self.S, [0.0,0.0,0.0] )
        self.s.placeParticle( self.S, [5e-6,5e-6,5e-6] )

        t = self.s.t
        for i in range( 5 ):
            self.s.step()
        self.failIf( t == self.s.t )

    def test_ThreeParticles( self ):
        self.s.placeParticle( self.S, [0.0,0.0,0.0] )
        self.s.placeParticle( self.S, [5e-6,5e-6,5e-6] )
        self.s.placeParticle( self.S, [1e-7,1e-7,1e-7] )

        t = self.s.t
        for i in range( 5 ):
            self.s.step()
        self.failIf( t == self.s.t )


    def test_ThreeParticlesInContact( self ):
        self.s.placeParticle( self.S, [0.0,0.0,0.0] )
        self.s.placeParticle( self.S, [1e-7,0.0,0.0] )

        # dummy
        self.s.placeParticle( self.S, [2e-7,0.0,0.0] )

        t = self.s.t
        for i in range( 5 ):
            self.s.step()
        self.failIf( t == self.s.t )

    def test_FourParticlesClose( self ):
        self.s.placeParticle( self.SS, [ 2e-8,0.0,0.0 ] )
        self.s.placeParticle( self.SS, [ 3.003e-8,0.0,0.0 ] )

        self.s.placeParticle( self.SS, [ 0.994e-8-5e-10, 0.0, 0.0 ] )

        t = self.s.t
        for i in range( 10 ):
            self.s.step()
        self.failIf( t == self.s.t )

    def test_immobile_is_immobile( self ):
        particleA = self.s.placeParticle( self.A, [0.0,0.0,0.0] )
        self.s.placeParticle( self.B, [1.5000001e-8,0.0,0.0] )

        initialPosition = particleA[1].position

        for i in range( 10 ):
            self.s.step()
        
        newPosition = particleA[1].position
        dist = self.s.distance( initialPosition, newPosition )

        self.failIf( dist != 0, 'initial pos: %s,\tnew pos: %s' %
                     ( initialPosition, newPosition ) )


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

        self.m.addSpecies(A)
        self.m.addSpecies(B)
        self.m.addSpecies(C)

        self.m.addReaction([A],    [B],    self.kf_1)
        self.m.addReaction([B],    [A],    self.kb_1)
        self.m.addReaction([A, B], [C],    self.kf_2)
        self.m.addReaction([C],    [A, B], self.kb_2)
        self.m.addReaction([C],    [],     self.kf_1)

        self.s.setModel(self.m)

        self.s.throwInParticles(A, 2)
        self.s.throwInParticles(B, 2)

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

        m1 = self.m.addPlanarSurface(origin=[L/2, L/2, 2*L/10],
                                     vectorX=[1, 0, 0],
                                     vectorY=[0, 1, 0],
                                     Lx=L/2,
                                     Ly=L/2,
                                     Lz=thickness,
                                     name='m1')
        self.m.addSpecies(A, m1)
        self.m.addSpecies(B, m1)
        self.m.addSpecies(C, m1)

        self.m.addReaction([(A, m1)],          [(B, m1)],          self.kf_1)
        self.m.addReaction([(B, m1)],          [(A, m1)],          self.kb_1)
        self.m.addReaction([(A, m1), (B, m1)], [(C, m1)],          self.kf_2)
        self.m.addReaction([(C, m1)],          [(A, m1), (B, m1)], self.kb_2)
        self.m.addReaction([(C, m1)],          [],                 self.kf_1)

        self.s.setModel(self.m)

        self.s.throwInParticles(A, 2, m1)
        self.s.throwInParticles(B, 2, m1)

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

        d = self.m.addCylindricalSurface(origin=[L/2, L/2, L/2],
                                         radius=radius,
                                         orientation=[0, 1, 0],
                                         size=L/2,
                                         name='d')

        self.m.addSpecies(A, d)
        self.m.addSpecies(B, d)
        self.m.addSpecies(C, d)

        self.m.addReaction([(A, d)],         [(B, d)],         self.kf_1)
        self.m.addReaction([(B, d)],         [(A, d)],         self.kb_1)
        self.m.addReaction([(A, d), (B, d)], [(C, d)],         self.kf_2)
        self.m.addReaction([(C, d)],         [(A, d), (B, d)], self.kb_2)
        self.m.addReaction([(C, d)],         [],               self.kf_1)

        self.s.setModel(self.m)

        self.s.throwInParticles(A, 2, d)
        self.s.throwInParticles(B, 2, d)

    def test_run(self):
        for i in range(10):
            self.s.step()


if __name__ == "__main__":
    unittest.main()
