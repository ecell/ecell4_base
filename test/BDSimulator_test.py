#!/usr/bin/env python

import unittest

import numpy

from _gfrd import World
from bd import *

class BDSimulatorTestCase( unittest.TestCase ):

    def setUp(self):
        self.m = ParticleModel()
        self.S = self.m.new_species_type( 'S', 2e-11, 5e-8 )
        self.A = self.m.new_species_type( 'A', 0, 1e-8 )
        self.B = self.m.new_species_type( 'B', 2e-11, 5e-9 )
        self.m.set_all_repulsive()
        world = World(1e-5, 10)
        self.s = BDSimulator(world)
        self.s.setModel(self.m)

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

    def test_immobile_is_immobile( self ):
        particleA = self.s.placeParticle( self.A, [0.0,0.0,0.0] )
        self.s.placeParticle( self.B, [1.5000001e-8,0.0,0.0] )

        initialPosition = particleA[1].position

        for i in range( 10 ):
            self.s.step()
            print particleA[1].position
        
        newPosition = particleA[1].position
        dist = self.s.distance( initialPosition, newPosition )

        self.failIf( dist != 0, 'initial pos: %s,\tnew pos: %s' %
                     ( initialPosition, newPosition ) )



if __name__ == "__main__":
    unittest.main()
