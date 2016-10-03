
4. How to Run a Simulation
==========================

In sections 2 and 3, we explain the way to build a model and to setup
the intial state. Now, it is the time to run a simulation. Corresponding
to ``World`` classes, six ``Simulator`` classes are there:
``spatiocyte.SpatiocyteSimulator``, ``egfrd.EGFRDSimulator``,
``bd.BDSimulator``, ``meso.MesoscopicSimulator``,
``gillespie.GillespieSimulator``, and ``ode.ODESimulator``. Each
``Simulator`` class only accepts the corresponding type of ``World``,
but all of them allow the same ``Model``.

.. code:: python

    import ecell4

4.1. How to Setup a Simulator
-----------------------------

Except for the initialization (so-called constructor function) with
arguments specific to the algorithm, all ``Simulator``\ s have the same
APIs.

.. code:: python

    from ecell4.core import *
    from ecell4.gillespie import GillespieWorld, GillespieSimulator
    from ecell4.ode import ODEWorld, ODESimulator
    from ecell4.spatiocyte import SpatiocyteWorld, SpatiocyteSimulator
    from ecell4.bd import BDWorld, BDSimulator
    from ecell4.meso import MesoscopicWorld, MesoscopicSimulator
    from ecell4.egfrd import EGFRDWorld, EGFRDSimulator

Before constructing a ``Simulator``, parepare a ``Model`` and a
``World`` corresponding to the type of ``Simulator``.

.. code:: python

    from ecell4 import species_attributes, reaction_rules, get_model
    
    with species_attributes():
        A | B | C | {'D': '1', 'radius': '0.005'}
    
    with reaction_rules():
        A + B == C | (0.01, 0.3)
    
    m = get_model()

.. code:: python

    w1 = GillespieWorld()
    w2 = ODEWorld()
    w3 = SpatiocyteWorld()
    w4 = BDWorld()
    w5 = MesoscopicWorld()
    w6 = EGFRDWorld()

``Simulator`` requires both ``Model`` and ``World`` in this order at the
construction.

.. code:: python

    sim1 = GillespieSimulator(m, w1)
    sim2 = ODESimulator(m, w2)
    sim3 = SpatiocyteSimulator(m, w3)
    sim4 = BDSimulator(m, w4)
    sim5 = MesoscopicSimulator(m, w5)
    sim6 = EGFRDSimulator(m, w6)

Once you bind a ``Model`` to a ``World``, only the ``World`` is needed
to create a ``Simulator``.

.. code:: python

    w1.bind_to(m)
    w2.bind_to(m)
    w3.bind_to(m)
    w4.bind_to(m)
    w5.bind_to(m)
    w6.bind_to(m)

.. code:: python

    sim1 = GillespieSimulator(w1)
    sim2 = ODESimulator(w2)
    sim3 = SpatiocyteSimulator(w3)
    sim4 = BDSimulator(w4)
    sim5 = MesoscopicSimulator(w5)
    sim6 = EGFRDSimulator(w6)

Of course, the ``Model`` and ``World`` bound to a ``Simulator`` can be
drawn from ``Simulator`` in the way below:

.. code:: python

    print(sim1.model(), sim1.world())
    print(sim2.model(), sim2.world())
    print(sim3.model(), sim3.world())
    print(sim4.model(), sim4.world())
    print(sim5.model(), sim5.world())
    print(sim6.model(), sim6.world())


.. parsed-literal::

    (<ecell4.core.Model object at 0x7fb7197c9bd0>, <ecell4.gillespie.GillespieWorld object at 0x7fb7197c9b58>)
    (<ecell4.ode.ODENetworkModel object at 0x7fb7197c9bd0>, <ecell4.ode.ODEWorld object at 0x7fb7197c9bb8>)
    (<ecell4.core.Model object at 0x7fb7197c9bd0>, <ecell4.spatiocyte.SpatiocyteWorld object at 0x7fb7197c9b58>)
    (<ecell4.core.Model object at 0x7fb7197c9bd0>, <ecell4.bd.BDWorld object at 0x7fb7197c9bb8>)
    (<ecell4.core.Model object at 0x7fb7197c9bd0>, <ecell4.meso.MesoscopicWorld object at 0x7fb7197c9ae0>)
    (<ecell4.core.Model object at 0x7fb7197c9bd0>, <ecell4.egfrd.EGFRDWorld object at 0x7fb7197c9bb8>)


After updating the ``World`` by yourself, you must initialize the
internal state of a ``Simulator`` before running simulation.

.. code:: python

    w1.add_molecules(Species('C'), 60)
    w2.add_molecules(Species('C'), 60)
    w3.add_molecules(Species('C'), 60)
    w4.add_molecules(Species('C'), 60)
    w5.add_molecules(Species('C'), 60)
    w6.add_molecules(Species('C'), 60)

.. code:: python

    sim1.initialize()
    sim2.initialize()
    sim3.initialize()
    sim4.initialize()
    sim5.initialize()
    sim6.initialize()

Algorithms with a fixed step interval also require ``dt``.

.. code:: python

    sim2.set_dt(1e-6)  # ODESimulator. This is optional
    sim4.set_dt(1e-6)  # BDSimulator

4.2. Running Simulation
-----------------------

For running simulation, ``Simulator`` provides two APIs, ``step`` and
``run``.

``step()`` advances a simulation for the time that the ``Simulator``
expects, ``next_time()``.

.. code:: python

    print(sim1.t(), sim1.next_time(), sim1.dt())
    print(sim2.t(), sim2.next_time(), sim2.dt())  # => (0.0, 1e-6, 1e-6)
    print(sim3.t(), sim3.next_time(), sim3.dt())
    print(sim4.t(), sim4.next_time(), sim4.dt())  # => (0.0, 1e-6, 1e-6)
    print(sim5.t(), sim5.next_time(), sim5.dt())
    print(sim6.t(), sim6.next_time(), sim6.dt())  # => (0.0, 0.0, 0.0)


.. parsed-literal::

    (0.0, 0.013012847696983692, 0.013012847696983692)
    (0.0, 1e-06, 1e-06)
    (0.0, 1.6666666666666667e-05, 1.6666666666666667e-05)
    (0.0, 1e-06, 1e-06)
    (0.0, 0.0045567279714810015, 0.0045567279714810015)
    (0.0, 0.0, 0.0)


.. code:: python

    sim1.step()
    sim2.step()
    sim3.step()
    sim4.step()
    sim5.step()
    sim6.step()

.. code:: python

    print(sim1.t())
    print(sim2.t())  # => 1e-6
    print(sim3.t())
    print(sim4.t())  # => 1e-6
    print(sim5.t())
    print(sim6.t())  # => 0.0


.. parsed-literal::

    0.013012847697
    1e-06
    1.66666666667e-05
    1e-06
    0.00455672797148
    0.0


``last_reactions()`` returns a list of pairs of ``ReactionRule`` and
``ReactionInfo`` which occurs at the last step. Each algorithm have its
own implementation of ``ReactionInfo``. See
``help(module.ReactionInfo)`` for details.

.. code:: python

    print(sim1.last_reactions())
    # print(sim2.last_reactions())
    print(sim3.last_reactions())
    print(sim4.last_reactions())
    print(sim5.last_reactions())
    print(sim6.last_reactions())


.. parsed-literal::

    [(<ecell4.core.ReactionRule object at 0x7fb7197c9bd0>, <ecell4.gillespie.ReactionInfo object at 0x7fb7197c9ae0>)]
    []
    []
    []
    []


``step(upto)`` advances a simulation for ``next_time`` if ``next_time``
is less than ``upto``, or for ``upto`` otherwise. ``step(upto)`` returns
whether the time does **NOT** reach the limit, ``upto``.

.. code:: python

    print(sim1.step(1.0), sim1.t())
    print(sim2.step(1.0), sim2.t())
    print(sim3.step(1.0), sim3.t())
    print(sim4.step(1.0), sim4.t())
    print(sim5.step(1.0), sim5.t())
    print(sim6.step(1.0), sim6.t())


.. parsed-literal::

    (True, 0.05959787878069947)
    (True, 2e-06)
    (True, 3.3333333333333335e-05)
    (True, 2e-06)
    (True, 0.00646251309677119)
    (True, 0.0)


Thus, for running a simulation just until the time, ``upto``, call
``step(upto)`` while it returns ``True``.

.. code:: python

    while sim1.step(1.0): pass
    while sim2.step(0.001): pass
    while sim3.step(0.001): pass
    while sim4.step(0.001): pass
    while sim5.step(1.0): pass
    while sim6.step(0.001): pass

.. code:: python

    print(sim1.t())  # => 1.0
    print(sim2.t())  # => 0.001
    print(sim3.t())  # => 0.001
    print(sim4.t())  # => 0.001
    print(sim5.t())  # => 1.0
    print(sim6.t())  # => 0.01


.. parsed-literal::

    1.0
    0.001
    0.001
    0.001
    1.0
    0.001


This is just what ``run`` does. ``run(tau)`` advances a simulation upto
``t()+tau``.

.. code:: python

    sim1.run(1.0)
    sim2.run(0.001)
    sim3.run(0.001)
    sim4.run(0.001)
    sim5.run(1.0)
    sim6.run(0.001)

.. code:: python

    print(sim1.t())  # => 2.0
    print(sim2.t())  # => 0.002
    print(sim3.t())  # => 0.002
    print(sim4.t())  # => 0.002
    print(sim5.t())  # => 2.0
    print(sim6.t())  # => 0.02


.. parsed-literal::

    2.0
    0.002
    0.002
    0.002
    2.0
    0.002


``num_steps`` returns the number of steps during the simulation.

.. code:: python

    print(sim1.num_steps())
    print(sim2.num_steps())
    print(sim3.num_steps())
    print(sim4.num_steps())
    print(sim5.num_steps())
    print(sim6.num_steps())


.. parsed-literal::

    39
    2001
    120
    2001
    998
    532


4.3. Capsulizing Algorithm into a Factory Class
-----------------------------------------------

Owing to the portability of a ``Model`` and consistent APIs of
``World``\ s and ``Simulator``\ s, it is very easy to write a script
common in algorithms. However, when switching the algorithm, still we
have to rewrite the name of classes in the code, one by one.

To avoid the trouble, E-Cell4 also provides a ``Factory`` class for each
algorithm. ``Factory`` encapsulates ``World`` and ``Simulator`` with
their arguments needed for the construction. By using ``Factory`` class,
your script could be portable and robust agaist changes in the
algorithm.

.. code:: python

    from ecell4.gillespie import GillespieFactory
    from ecell4.ode import ODEFactory
    from ecell4.spatiocyte import SpatiocyteFactory
    from ecell4.bd import BDFactory
    from ecell4.meso import MesoscopicFactory
    from ecell4.egfrd import EGFRDFactory

``Factory`` just provides two functions, ``create_world`` and
``create_simulator``.

.. code:: python

    def singlerun(f, m):
        w = f.create_world(Real3(1, 1, 1))
        w.bind_to(m)
        w.add_molecules(Species('C'), 60)
        
        sim = f.create_simulator(w)
        sim.run(0.01)
        print(sim.t(), w.num_molecules(Species('C')))

``singlerun`` above is free from the algorithm. Thus, by just switching
``Factory``, you can easily compare the results.

.. code:: python

    singlerun(GillespieFactory(), m)
    singlerun(ODEFactory(), m)
    singlerun(SpatiocyteFactory(), m)
    singlerun(BDFactory(bd_dt_factor=1), m)
    singlerun(MesoscopicFactory(), m)
    singlerun(EGFRDFactory(), m)


.. parsed-literal::

    (0.01, 59)
    (0.01, 59)
    (0.01, 60)
    (0.01, 60)
    (0.01, 60)
    (0.01, 59)


When you need to provide several parameters to initialize ``World`` or
``Simulator``, ``run_simulation`` also accepts ``Factory`` instead of
``solver``.

.. code:: python

    from ecell4.util import run_simulation
    print(run_simulation(0.01, model=m, y0={'C': 60}, return_type='array', factory=GillespieFactory())[-1])
    print(run_simulation(0.01, model=m, y0={'C': 60}, return_type='array', factory=ODEFactory())[-1])
    print(run_simulation(0.01, model=m, y0={'C': 60}, return_type='array', factory=SpatiocyteFactory())[-1])
    print(run_simulation(0.01, model=m, y0={'C': 60}, return_type='array', factory=BDFactory(bd_dt_factor=1))[-1])
    print(run_simulation(0.01, model=m, y0={'C': 60}, return_type='array', factory=MesoscopicFactory())[-1])
    print(run_simulation(0.01, model=m, y0={'C': 60}, return_type='array', factory=EGFRDFactory())[-1])


.. parsed-literal::

    [0.01, 0.0, 0.0, 60.0]
    [0.01, 0.17972919304001073, 0.17972919304001067, 59.82027080696036]
    [0.01, 0.0, 0.0, 60.0]
    [0.01, 0.0, 0.0, 60.0]
    [0.01, 0.0, 0.0, 60.0]
    [0.01, 0.0, 0.0, 60.0]

