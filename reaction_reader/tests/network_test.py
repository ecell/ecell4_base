import unittest

import network
from decorator2 import create_species, create_reaction_rule


class NetworkTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_check_stoichiometry(self):
        sp1 = create_species(
            "A(l,r^1).A(l^1,r^2).B(l^2,r^3).A(l^3,r^4)"
            + ".B(l^4,r^5).A(l^5,r^6).B(l^6,r)")
        self.assertTrue(
            network.check_stoichiometry(sp1, dict(A=4, B=3)))
        self.assertTrue(
            network.check_stoichiometry(sp1, dict(C=0)))
        self.assertFalse(
            network.check_stoichiometry(sp1, dict(A=3, B=3)))
        self.assertFalse(
            network.check_stoichiometry(sp1, dict(A=4, B=2)))


if __name__ == '__main__':
    unittest.main()
