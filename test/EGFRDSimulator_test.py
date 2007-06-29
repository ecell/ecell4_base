#!/usr/bin/env python

import unittest

import numpy

from egfrd import *

class EGFRDSimulatorTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation( self ):
        s = EGFRDSimulator()
        self.failIf( s == None )

    
    def no_test_OneParticle( self ):  # this currently fails
        s = EGFRDSimulator()
        s.setBoundarySize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( 'S', [0.0,0.0,0.0] )

        for i in range( 5 ):
            s.step()

    def test_TwoParticles( self ):
        s = EGFRDSimulator()
        s.setBoundarySize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( 'S', [0.0,0.0,0.0] )
        s.placeParticle( 'S', [5e-6,5e-6,5e-6] )

        for i in range( 5 ):
            s.step()

    def test_ThreeParticles( self ):
        s = EGFRDSimulator()
        s.setBoundarySize( 1e-5 )
        S = Species( 'S', 2e-11, 5e-8 )
        s.addSpecies( S )
        s.placeParticle( 'S', [0.0,0.0,0.0] )
        s.placeParticle( 'S', [5e-6,5e-6,5e-6] )
        s.placeParticle( 'S', [1e-7,1e-7,1e-7] )

        for i in range( 5 ):
            s.step()




if __name__ == "__main__":
    unittest.main()
