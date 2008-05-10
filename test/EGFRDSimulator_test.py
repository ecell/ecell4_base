#!/usr/bin/env python

import logging

import unittest

import numpy

from egfrd import *


#log.setLevel( logging.WARNING )

class EGFRDSimulatorTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation( self ):
        s = EGFRDSimulator()
        self.failIf( s == None )

    
    def test_OneParticle( self ):
        s = EGFRDSimulator()
        s.setWorldSize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( S, [0.0,0.0,0.0] )

        t = s.t
        for i in range( 5 ):
            s.step()
        self.failIf( t == s.t )

    def test_TwoParticles( self ):
        s = EGFRDSimulator()
        s.setWorldSize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( S, [0.0,0.0,0.0] )
        s.placeParticle( S, [5e-6,5e-6,5e-6] )

        t = s.t
        for i in range( 5 ):
            s.step()
        self.failIf( t == s.t )

    def test_ThreeParticles( self ):
        s = EGFRDSimulator()
        s.setWorldSize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( S, [0.0,0.0,0.0] )
        s.placeParticle( S, [5e-6,5e-6,5e-6] )
        s.placeParticle( S, [1e-7,1e-7,1e-7] )

        t = s.t
        for i in range( 5 ):
            s.step()
        self.failIf( t == s.t )


    def test_ThreeParticlesInContact( self ):
        s = EGFRDSimulator()
        s.setWorldSize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( S, [0.0,0.0,0.0] )
        s.placeParticle( S, [1e-7,0.0,0.0] )

        # dummy
        s.placeParticle( S, [2e-7,0.0,0.0] )

        t = s.t
        for i in range( 5 ):
            s.step()
        self.failIf( t == s.t )

    def test_FourParticlesClose( self ):

        #log.setLevel( logging.DEBUG )

        s = EGFRDSimulator()
        s.setWorldSize( 1e-5 )
        S = Species( 'S', 1e-12, 5e-9 )
        s.addSpecies( S )

        s.placeParticle( S, [ 2e-8,0.0,0.0 ] )
        s.placeParticle( S, [ 3.003e-8,0.0,0.0 ] )

        s.placeParticle( S, [ 0.994e-8-5e-10, 0.0, 0.0 ] )
        #s.placeParticle( S, [ 0, 2e-8, 0.0 ] )


        t = s.t
        for i in range( 10 ):
            s.step()
        self.failIf( t == s.t )

        #log.setLevel( logging.WARNING )


    def test_immobile_is_immobile( self ):
        s = EGFRDSimulator()
        s.setWorldSize( 1e-5 )
        A = Species( 'A', 0, 1e-8 )
        s.addSpecies( A )
        B = Species( 'B', 2e-11, 5e-9 )
        s.addSpecies( B )

        particleA = s.placeParticle( A, [0.0,0.0,0.0] )
        s.placeParticle( B, [1.5000001e-8,0.0,0.0] )

        initialPosition = particleA.getPos().copy()

        for i in range( 10 ):
            s.step()
        
        newPosition = particleA.getPos().copy()
        dist = s.distance( initialPosition, newPosition )

        self.failIf( dist != 0, 'initial pos: %s,\tnew pos: %s' %
                     ( initialPosition, newPosition ) )



if __name__ == "__main__":
    unittest.main()
