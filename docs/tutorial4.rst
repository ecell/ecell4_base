
Tutorial 4 (Simulator)
======================

This is a tutorial for E-Cell4. Here, we explain how to handle
Simulators.

Each World has its corresponding Simulator.

.. code:: python

    from ecell4.core import *
    # from ecell4.gillespie import GillespieWorld as world_type, GillespieSimulator as simulator_type
    # from ecell4.ode import ODEWorld as world_type, ODESimulator as simulator_type
    from ecell4.lattice import LatticeWorld as world_type, LatticeSimulator as simulator_type
    # from ecell4.bd import BDWorld as world_type, BDSimulator as simulator_type
Simulator needs a NetworkModel and World at the instantiation.

.. code:: python

    m = NetworkModel()
    m.add_species_attribute(Species("A", "0.0025", "1"))
    m.add_reaction_rule(create_degradation_reaction_rule(Species("A"), 0.693 / 1))
    
    w = world_type(Position3(1, 1, 1))
    w.bind_to(m)
    w.add_molecules(Species("A"), 60)
    
    sim = simulator_type(m, w)
    sim.set_dt(0.01) #XXX: Optional
A Simulator has getters for a simulation time, a step interval, and the
next-event time. In principle, a Simulator returns the World's time as
its simulation time, and does the sum of the current time and a step
interval as the next-event time.

.. code:: python

    print sim.num_steps()
    print sim.t(), w.t()
    print sim.next_time(), sim.t() + sim.dt()

.. parsed-literal::

    0
    0.0 0.0
    0.01 0.01


A Simulator can return the connected model and world. They are not
copies, but the shared objects.

.. code:: python

    print sim.model(), sim.world()

.. parsed-literal::

    <ecell4.core.NetworkModel object at 0x7f603c71d948> <ode.ODEWorld object at 0x7f603c71d978>


If you change a World after connecting it to a Simulator, you have to
call ``initialize()`` manually before ``step()``. The call will update
the internal state of the Simulator. (NOTE: This requirement will be
removed for the beta version.)

.. code:: python

    sim.world().add_molecules(Species("A"), 60) # w.add_molecules(Species("A"), 60)
    sim.initialize() #XXX: this should be called automatically
.. code:: python

    w.save('test.h5')
Simulator has two types of ``step`` functions. First, with no argument,
``step()`` increments the time until ``next_time()``.

.. code:: python

    print "%.3e %.3e" % (sim.t(), sim.next_time())
    sim.step()
    print "%.3e %.3e" % (sim.t(), sim.next_time())

.. parsed-literal::

    0.000e+00 1.000e-02
    1.000e-02 2.000e-02


With an argument, ``step(upto)`` increments the time if ``upto`` is less
than ``next_time()`` and returns ``True``. If not, it increments the
time for ``upto`` and returns ``False``. (If the current time ``t()`` is
less than ``upto``, it does nothing and returns ``False``.

.. code:: python

    print "%.3e %.3e" % (sim.t(), sim.next_time())
    print sim.step(0.1)
    print "%.3e %.3e" % (sim.t(), sim.next_time())

.. parsed-literal::

    1.000e-02 2.000e-02
    False
    1.000e-01 1.100e-01


For a discrete-step simulation, the loop can be written like:

.. code:: python

    w.load('test.h5')
    sim.initialize()
.. code:: python

    next_time, dt = 0.0, 1e-2
    for _ in range(5):
        while sim.step(next_time): pass
        next_time += dt
    
        print "%.3e %.3e %d %g" % (sim.t(), sim.dt(), sim.num_steps(), w.num_molecules(Species("A")))

.. parsed-literal::

    0.000e+00 1.000e-02 2 120
    1.000e-02 1.000e-02 3 119.171
    2.000e-02 1.000e-02 4 118.348
    3.000e-02 1.000e-02 5 117.531
    4.000e-02 1.000e-02 6 116.719


