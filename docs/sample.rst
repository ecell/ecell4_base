
E-Cell4 Test1
=============

This is a simple test for E-Cell System version 4.

First, you need to import the Core module:

.. code:: python

    from ecell4.core import *
A simple equilibrium model
(:math:`\mathrm{A}+\mathrm{B}\leftrightarrows\mathrm{AB}`) with an
association rate :math:`k_a` and a dissociation rate :math:`k_d` is
tested here.

For the parameters, an equilibrium state
:math:`k_a\frac{N}{V}U^2=k_d\left(1-U\right)` gives
:math:`k_{a}=k_d\cdot\frac{1-U}{U^2}/\frac{N}{V}`.

.. code:: python

    L, N, kd, U = 1e-6, 60, 0.1, 0.5
    D, radius = "1e-12", "1e-8"
    volume = L * L * L
    ka = kd * volume * (1 - U) / (U * U * N)
    sp1, sp2, sp3 = Species("A"), Species("B"), Species("AB")
    sp1.set_attribute("radius", radius)
    sp1.set_attribute("D", D)
    sp2.set_attribute("radius", radius)
    sp2.set_attribute("D", D)
    sp3.set_attribute("radius", radius)
    sp3.set_attribute("D", D)
    rr1, rr2 = create_binding_reaction_rule(sp1, sp2, sp3, ka), create_unbinding_reaction_rule(sp3, sp1, sp2, kd)
Create a model:

.. code:: python

    m = NetworkModel()
    m.add_species_attribute(sp1)
    m.add_species_attribute(sp2)
    m.add_species_attribute(sp3)
    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)
Here, we use the Gillespie or ODE module.

.. code:: python

    from ecell4.gillespie import GillespieWorld as world_type, GillespieSimulator as simulator_type
    # from ecell4.ode import ODEWorld as world_type, ODESimulator as simulator_type
.. code:: python

    w = world_type(volume)
    w.add_molecules(sp1, N)
    w.add_molecules(sp2, N)
    sim = simulator_type(m, w)
Run a simulation:

.. code:: python

    next_time, dt = 0.0, 0.1
    data = [(w.t(), w.num_molecules(sp1), w.num_molecules(sp2), w.num_molecules(sp3))]
    for i in range(110):
        next_time += dt
        while (sim.step(next_time)): pass
        data.append((w.t(), w.num_molecules(sp1), w.num_molecules(sp2), w.num_molecules(sp3)))
Plot with Matplotlib:

.. code:: python

    import matplotlib.pylab as plt
    import numpy
    data = numpy.array(data)
    plt.plot(data.T[0], data.T[1], "g-", label=sp1.name())
    plt.plot(data.T[0], data.T[3], "r-", label=sp3.name())
    plt.xlabel("Time")
    plt.ylabel("Number Of Molecules")
    plt.xlim(data.T[0][0], data.T[0][-1])
    plt.legend(loc="best", shadow=True)
    plt.show()


.. image:: https://raw.githubusercontent.com/ecell/ecell4/master/docs/images/timeseries.png


A simple VTK viewer is also available:

.. code:: python

    from ecell4.bd import BDWorld
    w = BDWorld(Position3(L, L, L))
    w.add_molecules(sp1, int(round(data[-1][1])))
    w.add_molecules(sp2, int(round(data[-1][2])))
    w.add_molecules(sp3, int(round(data[-1][3])))
.. code:: python

    import vtk_test
    print "Number Of Molecules: A=%g, B=%g, AB=%g" % (w.num_molecules(sp1), w.num_molecules(sp2), w.num_molecules(sp3))
    vtk_test.show(w, 500, 500)

.. parsed-literal::

    Number Of Molecules: A=31, B=31, AB=29


.. image:: https://raw.githubusercontent.com/ecell/ecell4/master/docs/images/particle.png
