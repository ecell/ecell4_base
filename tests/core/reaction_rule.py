import unittest
import copy
from ecell4_base.core import *
from . import setUpEqualities

class ReactionRuleDescriptor(unittest.TestCase):
    def setUp(self):
        setUpEqualities(self)

    def test_mass_action(self):
        d = ReactionRuleDescriptorMassAction(1.0)
        self.assertEqual(d.k(), 1.0)
        self.assertEqual(d.get_k(), Quantity_Real(1.0, ''))

        d.set_k(2.0)
        self.assertEqual(d.k(), 2.0)
        self.assertEqual(d.get_k(), Quantity_Real(2.0, ''))

        d.set_k(Quantity_Real(3.0, 'm'))
        self.assertEqual(d.k(), 3.0)
        self.assertEqual(d.get_k(), Quantity_Real(3.0, 'm'))

        self.assertEqual(d.reactant_coefficients(), [])
        self.assertEqual(d.product_coefficients(), [])

    def test_pyfunc(self):
        d = ReactionRuleDescriptorPyfunc(
                lambda reactants, products, volume, t, reactant_coefs, product_coefs: 3.14,
                "pyfunc_descriptor")
        self.assertEqual(d.as_string(), "pyfunc_descriptor")

        self.assertEqual(d.propensity([], [], 1.0, 1.0e-5), 3.14)

        self.assertEqual(d.reactant_coefficients(), [])
        self.assertEqual(d.product_coefficients(), [])


class ReactionRuleTest(unittest.TestCase):
    def setUp(self):
        setUpEqualities(self)

    def test_constructor(self):
        reactant = Species("A")
        product = Species("B")

        rr = ReactionRule()
        rr = ReactionRule([reactant], [product])
        rr = ReactionRule([reactant], [product], 1.0)
        rr = ReactionRule([reactant], [product], Quantity_Real(1.0))
        rr = ReactionRule([reactant], [product], Quantity_Real(1.0, "/s"))

    def test_getters(self):
        reactant = Species("A")
        product = Species("B")
        quantity = Quantity_Real(1.5, "/s")
        rr = ReactionRule([reactant], [product], quantity)

        self.assertEqual(rr.k(), 1.5)
        self.assertEqual(rr.get_k(), quantity)
        self.assertEqual(rr.reactants(), [reactant])
        self.assertEqual(rr.products(), [product])
        self.assertEqual(rr.as_string(), "A>B|1.5")

        self.assertEqual(rr.count([]), 0)
        self.assertEqual(rr.count([reactant]), 1)
        self.assertEqual(rr.count([product]), 0)
        self.assertEqual(rr.count([reactant, product]), 0)

    def test_policy(self):
        rr = ReactionRule()
        self.assertEqual(rr.policy(), ReactionRulePolicy.STRICT)

        rr.set_policy(ReactionRulePolicy.IMPLICIT)
        self.assertEqual(rr.policy(), ReactionRulePolicy.IMPLICIT)

        rr.set_policy(ReactionRulePolicy.DESTROY)
        self.assertEqual(rr.policy(), ReactionRulePolicy.DESTROY)


    def test_descriptor(self):
        rr = ReactionRule()
        self.assertFalse(rr.has_descriptor())
        self.assertEqual(rr.get_descriptor(), None)

        # rr.set_descriptor()

    def test_pickling(self):
        k = Quantity_Real(2.0, 'micrometer ** 2 / second')
        policy = DESTROY

        rr = ReactionRule()
        rr.add_reactant(Species('A'))
        rr.add_reactant(Species('B'))
        rr.add_product(Species('C'))
        rr.set_k(k)
        rr.set_policy(policy)
        rr.set_attribute("key1", "value1")
        rr.set_attribute("key2", Quantity_Real(2.0, "units"))
        rr.set_attribute("key3", True)

        self.assertEqual(len(rr.reactants()), 2)
        self.assertEqual(len(rr.products()), 1)
        self.assertEqual(rr.get_k(), k)
        self.assertEqual(rr.k(), k.magnitude)
        self.assertFalse(rr.has_descriptor())
        self.assertEqual(rr.policy(), policy)
        self.assertTrue(rr.has_attribute("key1"))
        self.assertEqual(rr.get_attribute("key1"), "value1")
        self.assertTrue(rr.has_attribute("key2"))
        self.assertEqual(rr.get_attribute("key2"), Quantity_Real(2.0, "units"))
        self.assertTrue(rr.has_attribute("key3"))
        self.assertEqual(rr.get_attribute("key3"), True)

        another = copy.copy(rr)
        self.assertTrue(rr is not another)
        self.assertEqual(rr, another)

        self.assertEqual(len(another.reactants()), 2)
        self.assertEqual(len(another.products()), 1)
        self.assertEqual(another.get_k(), k)
        self.assertEqual(another.k(), k.magnitude)
        self.assertFalse(another.has_descriptor())
        self.assertEqual(another.policy(), policy)
        self.assertTrue(another.has_attribute("key1"))
        self.assertEqual(another.get_attribute("key1"), "value1")
        self.assertTrue(another.has_attribute("key2"))
        self.assertEqual(another.get_attribute("key2"), Quantity_Real(2.0, "units"))
        self.assertTrue(another.has_attribute("key3"))
        self.assertEqual(another.get_attribute("key3"), True)
