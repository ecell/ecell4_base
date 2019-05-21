import unittest
from ecell4_base.core import *

class ReactionRuleDescriptorPyfuncTest(unittest.TestCase):

    def test_clone(self):
        m = NetworkModel()
        rr = create_binding_reaction_rule(Species("A"), Species("B"), Species("C"), 0.0)
        desc = ReactionRuleDescriptorPyfunc(lambda r, p, v, t, rc, pc: 0.1 * r[0] * r[1], "test")
        rr.set_descriptor(desc)
        m.add_reaction_rule(rr)

        self.assertTrue(rr.has_descriptor())
        self.assertTrue(m.reaction_rules()[0].has_descriptor())
