import unittest
from ecell4 import *

class ParticleIDTest(unittest.TestCase):
    def test_constructor(self):
        pid = ParticleID()
        pid = ParticleID((1, 2))

    def test_lot_serial(self):
        pid = ParticleID()
        self.assertEqual(pid.lot(), 0)
        self.assertEqual(pid.serial(), 0)

        pid = ParticleID((1, 2))
        self.assertEqual(pid.lot(), 1)
        self.assertEqual(pid.serial(), 2)

class ParticleTest(unittest.TestCase):
    def test_constructor(self):
        p = Particle(Species("A"), Real3(1,2,3), 1.0e-5, 1.0e-12)

    def test_variables(self):
        p = Particle(Species("A"), Real3(1,2,3), 1.0e-5, 1.0e-12)
        self.assertEqual(p.species().serial(), "A")
        self.assertEqual(p.position(), Real3(1,2,3))
        self.assertEqual(p.radius(), 1.0e-5)
        self.assertEqual(p.D(), 1.0e-12)
