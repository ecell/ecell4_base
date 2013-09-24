import unittest

import species
from decorator2 import create_reaction_rule


class ReactionRuleTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test_creation(self):
        rr1 = create_reaction_rule("A + B > C")
        self.assertTrue(rr1 is not None)


if __name__ == '__main__':
    unittest.main()
