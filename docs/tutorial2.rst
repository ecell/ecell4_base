
Tutorial 2 (Model Basics)
=========================

This is a tutorial for E-Cell4. Here, we represent the way of modeling
with E-Cell 4.

In E-Cell4, NetworkModel consists of a set of ReactionRules (and Species
attributes). First, instantiate NetworkModel.

.. code:: python

    from ecell4.core import *
    m = NetworkModel()
Here, you can add five kinds of ReactionRules.

1. create\_binding\_reaction\_rule
2. create\_degradation\_reaction\_rule
3. create\_synthesis\_reaction\_rule
4. create\_unbinding\_reaction\_rule
5. create\_unimolecular\_reaction\_rule

Now, for the simple binding/unbinding equilibrium model, let's use 1 and
4.

.. code:: python

    rr1 = create_binding_reaction_rule(Species("A"), Species("B"), Species("A.B"), 1.0)
    rr2 = create_unbinding_reaction_rule(Species("A.B"), Species("A"), Species("B"), 1.0)
A ReactionRule consists of reactants, products, and a kinetic rate. Of
course, you can change the kinetic rate.

.. code:: python

    print len(rr1.reactants()), len(rr1.products()), rr1.k()
    print len(rr2.reactants()), len(rr2.products()), rr2.k()
    rr2.set_k(2.0)
    print rr2.k()

.. parsed-literal::

    2 1 1.0
    1 2 1.0
    2.0


Next, add these ReactionRules into a NetworkModel.

.. code:: python

    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)
    print m.num_reaction_rules()
    print m.has_reaction_rule(rr1), m.has_reaction_rule(rr2)

.. parsed-literal::

    2
    True True


NetworkModel returns appropriate ReactionRules from the given Species
through APIs named query\_reaction\_rules.

.. code:: python

    print m.query_reaction_rules(Species("A.B"))
    print m.query_reaction_rules(Species("A"))
    print m.query_reaction_rules(Species("A"), Species("B"))
    print m.query_reaction_rules(Species("B"), Species("A"))

.. parsed-literal::

    [<ecell4.core.ReactionRule object at 0x7f38e7745978>]
    []
    [<ecell4.core.ReactionRule object at 0x7f38e7745990>]
    [<ecell4.core.ReactionRule object at 0x7f38e77459a8>]


NetworkModel also contains Species attributes. These attributes are
indispensable for particle and lattice simulations, but not necessarily
needed for gillespie and ode.

.. code:: python

    sp1 = Species("A")
    sp1.set_attribute("radius", "0.0025")
    sp1.set_attribute("D", "1")
    m.add_species_attribute(sp1)
    # m.add_species_attribute(Species("A", "0.0025", "1"))
NetworkModel attributes a Species based on the registered Species.

.. code:: python

    print m.has_species_attribute(Species("A"))
    sp2 = m.apply_species_attributes(Species("A"))
    print sp2.has_attribute("radius"), sp2.has_attribute("D")
    print sp2.get_attribute("radius"), sp2.get_attribute("D")

.. parsed-literal::

    True
    True True
    0.0025 1


