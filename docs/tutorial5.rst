
Tutorial 5 (Example 1)
======================

This is a tutorial for E-Cell4. Here, we show a simple equilibrium model
as an example.

.. code:: python

    %matplotlib inline
.. code:: python

    from ecell4.core import *
    from ecell4.gillespie import GillespieWorld as world_type, GillespieSimulator as simulator_type
    # from ecell4.ode import ODEWorld as world_type, ODESimulator as simulator_type
    # from ecell4.lattice import LatticeWorld as world_type, LatticeSimulator as simulator_type
    # from ecell4.bd import BDWorld as world_type, BDSimulator as simulator_type
Set up parameters:

.. code:: python

    L, N, kd, U, D, radius = 1.0, 60, 0.1, 0.5, "1", "0.01"
    volume = L * L * L
    ka = kd * volume * (1 - U) / (U * U * N)
    
    sp1, sp2, sp3 = Species("A", radius, D), Species("B", radius, D), Species("A.B", radius, D)
    rr1, rr2 = create_binding_reaction_rule(sp1, sp2, sp3, ka), create_unbinding_reaction_rule(sp3, sp1, sp2, kd)
Create a model:

.. code:: python

    m = NetworkModel()
    m.add_species_attribute(sp1)
    m.add_species_attribute(sp2)
    m.add_species_attribute(sp3)
    m.add_reaction_rule(rr1)
    m.add_reaction_rule(rr2)
Create a world and simulator:

.. code:: python

    w = world_type(Real3(L, L, L))
    w.bind_to(m)
    w.add_molecules(Species("A"), N)
    w.add_molecules(Species("B"), N)
    
    sim = simulator_type(m, w)
    sim.set_dt(1e-3) #XXX: This is too large to get the accurate result with BDSimulator.
Run a simulation:

.. code:: python

    next_time, dt = 0.0, 0.05
    data = [(w.t(), w.num_molecules(sp1), w.num_molecules(sp2), w.num_molecules(sp3))]
    for i in range(100):
        next_time += dt
        while (sim.step(next_time)): pass
        data.append((w.t(), w.num_molecules(sp1), w.num_molecules(sp2), w.num_molecules(sp3)))
Plot with Matplotlib:

.. code:: python

    import matplotlib.pylab as plt
    from numpy import array
    
    data = array(data)
    plt.plot(data.T[0], data.T[1], "r-", label=sp1.serial())
    plt.plot(data.T[0], data.T[2], "g--", label=sp2.serial())
    plt.plot(data.T[0], data.T[3], "b:", label=sp3.serial())
    plt.xlabel("Time")
    plt.ylabel("Number Of Molecules")
    plt.xlim(data.T[0][0], data.T[0][-1])
    plt.legend(loc="best", shadow=True)
    plt.show()


.. image:: tutorial5_files/tutorial5_12_0.png


