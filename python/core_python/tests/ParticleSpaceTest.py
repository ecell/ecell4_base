from ecell4.core import *

import unittest


class ParticleSpaceTest(unittest.TestCase):

    def setUp(self):
        self.ParticleSpace = ParticleSpaceVectorImpl

    def test_constructor(self):
        space = self.ParticleSpace(Position3(1e-6, 1e-6, 1e-6))

    def test_edge_lengths(self):
        space = self.ParticleSpace(Position3(1e-6, 1e-6, 1e-6))
        lengths = space.edge_lengths()
        self.assertEqual(lengths[0], 1e-6)
        self.assertEqual(lengths[1], 1e-6)
        self.assertEqual(lengths[2], 1e-6)

    def test_update_particle(self):
        space = self.ParticleSpace(Position3(1e-6, 1e-6, 1e-6))
        self.assertEqual(len(space.list_particles()), 0)
        space.update_particle(
            ParticleID(1, 2),
            Particle(Species("A"), Position3(0, 0, 0), 2.5e-9, 0))
        self.assertEqual(len(space.list_particles()), 1)
        pid, p = space.list_particles()[0]
        self.assertEqual(pid.lot(), 1)
        self.assertEqual(pid.serial(), 2)


if __name__ == '__main__':
    unittest.main()
