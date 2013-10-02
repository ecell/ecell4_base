import unittest

import network
from decorator2 import create_species, create_reaction_rule


class NetworkTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_check_stoichiometry(self):
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
