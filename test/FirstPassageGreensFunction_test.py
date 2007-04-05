#!/usr/bin/env python

import unittest

import _gfrd as mod

class FirstPassageGreensFunctionTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def testInstantiation(self):
        D = 1e-12
        gf = mod.FirstPassageGreensFunction( D )
        self.failIf( gf == None )
        
if __name__ == "__main__":
    unittest.main()
