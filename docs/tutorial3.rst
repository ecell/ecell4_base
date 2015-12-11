
3. How to Setup the Initial Condition
=====================================

Here, we explain the basics of ``World`` classes. In E-Cell4, six types
of World classes are offically supported now:
``spatiocyte.SpatiocyteWorld``, ``egfrd.EGFRDWorld``, ``bd.BDWorld``,
``meso.MesoscopicWorld``, ``gillespie.GillespieWorld``, and
``ode.ODEWorld``.

In the most of softwares, the initial condition is supposed to be a part
of ``Model``. But, in E-Cell4, the initial condition must be set up as
``World`` separately from ``Model``. ``World`` stores an information
about the state at a time point, such as a current time, the number of
molecules, coordinate of molecules, structures, and random number
generator. Meanwhile, ``Model`` contains the type of interactions
between molecules and the common properties of molecules.

.. code:: python

    import ecell4

3.1. Common APIs of World
-------------------------

Even though ``World`` describes the spatial representation specific to
the corresponding algorithm, it has compatible APIs. In this section,
the common interfaces of six ``World`` classes are introduced.

.. code:: python

    from ecell4.core import *
    from ecell4.gillespie import GillespieWorld
    from ecell4.ode import ODEWorld
    from ecell4.spatiocyte import SpatiocyteWorld
    from ecell4.bd import BDWorld
    from ecell4.meso import MesoscopicWorld
    from ecell4.egfrd import EGFRDWorld

``World`` classes accept different sets of arguments in the constructor,
which determine the parameters specific to the algorithm. However, at
least, all ``World`` classes can be instantiated only with their size,
named ``edge_lengths``. The type of ``edge_lengths`` is ``Real3``, which
represents a triplet of ``Real``\ s. In E-Cell4, all 3-dimensional
positions are treated as a ``Real3``.

.. code:: python

    edge_lengths = Real3(1, 2, 3)
    w1 = GillespieWorld(edge_lengths)
    w2 = ODEWorld(edge_lengths)
    w3 = SpatiocyteWorld(edge_lengths)
    w4 = BDWorld(edge_lengths)
    w5 = MesoscopicWorld(edge_lengths)
    w6 = EGFRDWorld(edge_lengths)

``World`` has getter methods for the size and volume.

.. code:: python

    print(tuple(w1.edge_lengths()), w1.volume())
    print(tuple(w2.edge_lengths()), w2.volume())
    print(tuple(w3.edge_lengths()), w3.volume())
    print(tuple(w4.edge_lengths()), w4.volume())
    print(tuple(w5.edge_lengths()), w5.volume())
    print(tuple(w6.edge_lengths()), w6.volume())


.. parsed-literal::

    ((1.0, 2.0, 3.0), 6.0)
    ((1.0, 2.0, 3.0), 6.0)
    ((1.0, 2.0, 3.0), 6.0)
    ((1.0, 2.0, 3.0), 6.0)
    ((1.0, 2.0, 3.0), 6.0)
    ((1.0, 2.0, 3.0), 6.0)


Next, let's add molecules into the World. Here, you must give
``Species`` attributed with "radius" and "D" for ``EGFRDWorld``,
``BDWorld`` or ``SpatiocyteWorld`` to tell the shape of molecules.
Positions of the molecules are randomly determined in the ``World`` if
needed.

.. code:: python

    sp1 = Species("A", "0.0025", "1")
    w1.add_molecules(sp1, 10)
    w2.add_molecules(sp1, 10)
    w3.add_molecules(sp1, 10)
    w4.add_molecules(sp1, 10)
    w5.add_molecules(sp1, 10)
    w6.add_molecules(sp1, 10)

Once binding a ``NetworkModel`` to the ``World``, you don't need to give
attributes explicitly. The ``World`` will ask attributes to the bound
``NetworkModel``.

.. code:: python

    m = NetworkModel()
    m.add_species_attribute(Species("A", "0.0025", "1"))
    m.add_species_attribute(Species("B", "0.0025", "1"))
    
    w1.bind_to(m)
    w2.bind_to(m)
    w3.bind_to(m)
    w4.bind_to(m)
    w5.bind_to(m)
    w6.bind_to(m)
    w1.add_molecules(Species("B"), 20)
    w2.add_molecules(Species("B"), 20)
    w3.add_molecules(Species("B"), 20)
    w4.add_molecules(Species("B"), 20)
    w5.add_molecules(Species("B"), 20)
    w6.add_molecules(Species("B"), 20)

Similarly, ``remove_molecules`` and ``num_molecules_exact`` are also
available.

.. code:: python

    w1.remove_molecules(Species("B"), 5)
    w2.remove_molecules(Species("B"), 5)
    w3.remove_molecules(Species("B"), 5)
    w4.remove_molecules(Species("B"), 5)
    w5.remove_molecules(Species("B"), 5)
    w6.remove_molecules(Species("B"), 5)

.. code:: python

    print(w1.num_molecules_exact(Species("A")), w2.num_molecules_exact(Species("A")), w3.num_molecules_exact(Species("A")), w4.num_molecules_exact(Species("A")), w5.num_molecules_exact(Species("A")), w6.num_molecules_exact(Species("A")))
    print(w1.num_molecules_exact(Species("B")), w2.num_molecules_exact(Species("B")), w3.num_molecules_exact(Species("B")), w4.num_molecules_exact(Species("B")), w5.num_molecules_exact(Species("B")), w6.num_molecules_exact(Species("B")))


.. parsed-literal::

    (10, 10, 10, 10, 10, 10)
    (15, 15, 15, 15, 15, 15)


``num_molecules`` also count the number of molecules, but returns all
the number matched with the given ``Species`` in the rule-based way.
When all ``Species`` in the ``World`` has no site and bond,
``num_molecules`` is almost same with ``num_molecules_exact``.

.. code:: python

    print(w1.num_molecules(Species("A")), w2.num_molecules(Species("A")), w3.num_molecules(Species("A")), w4.num_molecules(Species("A")), w5.num_molecules(Species("A")), w6.num_molecules(Species("A")))
    print(w1.num_molecules(Species("B")), w2.num_molecules(Species("B")), w3.num_molecules(Species("B")), w4.num_molecules(Species("B")), w5.num_molecules(Species("B")), w6.num_molecules(Species("B")))


.. parsed-literal::

    (10, 10, 10, 10, 10, 10)
    (15, 15, 15, 15, 15, 15)


``World`` also owns a simulation time.

.. code:: python

    print(w1.t(), w2.t(), w3.t(), w4.t(), w5.t(), w6.t())
    w1.set_t(1.0)
    w2.set_t(1.0)
    w3.set_t(1.0)
    w4.set_t(1.0)
    w5.set_t(1.0)
    w6.set_t(1.0)
    print(w1.t(), w2.t(), w3.t(), w4.t(), w5.t(), w6.t())


.. parsed-literal::

    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)


Finally, you can ``save`` and ``load`` the state of a ``World``
into/from a HDF5 file.

.. code:: python

    w1.save("data/gillespie.h5")
    w2.save("data/ode.h5")
    w3.save("data/spatiocyte.h5")
    w4.save("data/bd.h5")
    w5.save("data/meso.h5")
    w6.save("data/egfrd.h5")
    del w1, w2, w3, w4, w5, w6

.. code:: python

    w1 = GillespieWorld()
    w2 = ODEWorld()
    w3 = SpatiocyteWorld()
    w4 = BDWorld()
    w5 = MesoscopicWorld()
    w6 = EGFRDWorld()
    print(w1.t(), tuple(w1.edge_lengths()), w1.volume(), w1.num_molecules(Species("A")), w1.num_molecules(Species("B")))
    print(w2.t(), tuple(w2.edge_lengths()), w2.volume(), w2.num_molecules(Species("A")), w2.num_molecules(Species("B")))
    print(w3.t(), tuple(w3.edge_lengths()), w3.volume(), w3.num_molecules(Species("A")), w3.num_molecules(Species("B")))
    print(w4.t(), tuple(w4.edge_lengths()), w4.volume(), w4.num_molecules(Species("A")), w4.num_molecules(Species("B")))
    print(w5.t(), tuple(w5.edge_lengths()), w5.volume(), w5.num_molecules(Species("A")), w5.num_molecules(Species("B")))
    print(w6.t(), tuple(w6.edge_lengths()), w6.volume(), w6.num_molecules(Species("A")), w6.num_molecules(Species("B")))


.. parsed-literal::

    (0.0, (1.0, 1.0, 1.0), 1.0, 0, 0)
    (0.0, (1.0, 1.0, 1.0), 1.0, 0, 0)
    (0.0, (1.0, 1.0, 1.0), 1.0, 0, 0)
    (0.0, (1.0, 1.0, 1.0), 1.0, 0, 0)
    (0.0, (1.0, 1.0, 1.0), 1.0, 0, 0)
    (0.0, (1.0, 1.0, 1.0), 1.0, 0, 0)


.. code:: python

    w1.load("data/gillespie.h5")
    w2.load("data/ode.h5")
    w3.load("data/spatiocyte.h5")
    w4.load("data/bd.h5")
    w5.load("data/meso.h5")
    w6.load("data/egfrd.h5")
    print(w1.t(), tuple(w1.edge_lengths()), w1.volume(), w1.num_molecules(Species("A")), w1.num_molecules(Species("B")))
    print(w2.t(), tuple(w2.edge_lengths()), w2.volume(), w2.num_molecules(Species("A")), w2.num_molecules(Species("B")))
    print(w3.t(), tuple(w3.edge_lengths()), w3.volume(), w3.num_molecules(Species("A")), w3.num_molecules(Species("B")))
    print(w4.t(), tuple(w4.edge_lengths()), w4.volume(), w4.num_molecules(Species("A")), w4.num_molecules(Species("B")))
    print(w5.t(), tuple(w5.edge_lengths()), w5.volume(), w5.num_molecules(Species("A")), w5.num_molecules(Species("B")))
    print(w6.t(), tuple(w6.edge_lengths()), w6.volume(), w6.num_molecules(Species("A")), w6.num_molecules(Species("B")))
    del w1, w2, w3, w4, w5, w6


.. parsed-literal::

    (1.0, (1.0, 2.0, 3.0), 6.0, 10, 15)
    (1.0, (1.0, 2.0, 3.0), 6.0, 10, 15)
    (1.0, (1.0, 2.0, 3.0), 6.0, 10, 15)
    (1.0, (1.0, 2.0, 3.0), 6.0, 10, 15)
    (1.0, (1.0, 2.0, 3.0), 6.0, 10, 15)
    (1.0, (1.0, 2.0, 3.0), 6.0, 10, 15)


All the ``World`` classes also accept a HDF5 file name as an unique
argument of the constructor.

.. code:: python

    print(GillespieWorld("data/gillespie.h5").t())
    print(ODEWorld("data/ode.h5").t())
    print(SpatiocyteWorld("data/spatiocyte.h5").t())
    print(BDWorld("data/bd.h5").t())
    print(MesoscopicWorld("data/meso.h5").t())
    print(EGFRDWorld("data/egfrd.h5").t())


.. parsed-literal::

    1.0
    1.0
    1.0
    1.0
    1.0
    1.0


3.2. How to Get Positions of Molecules
--------------------------------------

``World`` has the common functions to access coordinates of molecules
too. However, as the resolution in representation of molecules differs,
``World`` can automatically make up for the gap with random numbers. For
example, ``GillespieWorld`` has almost no information about the
coordinate of molecules. Thus, it returns a list of positions randomly
distributed in the ``World`` size.

.. code:: python

    w1 = GillespieWorld()
    w2 = ODEWorld()
    w3 = SpatiocyteWorld()
    w4 = BDWorld()
    w5 = MesoscopicWorld()
    w6 = EGFRDWorld()
    
    sp1 = Species("A", "0.0025", "1")
    w1.add_molecules(sp1, 3)
    w2.add_molecules(sp1, 3)
    w3.add_molecules(sp1, 3)
    w4.add_molecules(sp1, 3)
    w5.add_molecules(sp1, 3)
    w6.add_molecules(sp1, 3)

``list_particles`` and ``list_particles_exact`` return a list of pairs
of ``ParticleID`` and ``Particle`` in the ``World``.

.. code:: python

    print(w1.list_particles_exact(Species("A")))
    # print(w2.list_particles_exact(Species("A")))  # ODEWorld has no member named list_particles
    print(w3.list_particles_exact(Species("A")))
    print(w4.list_particles_exact(Species("A")))
    print(w5.list_particles_exact(Species("A")))
    print(w6.list_particles_exact(Species("A")))


.. parsed-literal::

    [(<ecell4.core.ParticleID object at 0x7f63db093930>, <ecell4.core.Particle object at 0x7f63db093a38>), (<ecell4.core.ParticleID object at 0x7f63db0939f0>, <ecell4.core.Particle object at 0x7f63db093ac8>), (<ecell4.core.ParticleID object at 0x7f63db093a68>, <ecell4.core.Particle object at 0x7f63db093af8>)]
    [(<ecell4.core.ParticleID object at 0x7f63db093a38>, <ecell4.core.Particle object at 0x7f63db093a68>), (<ecell4.core.ParticleID object at 0x7f63db0939f0>, <ecell4.core.Particle object at 0x7f63db0939d8>), (<ecell4.core.ParticleID object at 0x7f63db093ac8>, <ecell4.core.Particle object at 0x7f63db093ab0>)]
    [(<ecell4.core.ParticleID object at 0x7f63db093a68>, <ecell4.core.Particle object at 0x7f63db093ac8>), (<ecell4.core.ParticleID object at 0x7f63db0939f0>, <ecell4.core.Particle object at 0x7f63db093930>), (<ecell4.core.ParticleID object at 0x7f63db0939d8>, <ecell4.core.Particle object at 0x7f63db093ae0>)]
    [(<ecell4.core.ParticleID object at 0x7f63db093ac8>, <ecell4.core.Particle object at 0x7f63db0939d8>), (<ecell4.core.ParticleID object at 0x7f63db0939f0>, <ecell4.core.Particle object at 0x7f63db093a38>), (<ecell4.core.ParticleID object at 0x7f63db093930>, <ecell4.core.Particle object at 0x7f63db093af8>)]
    [(<ecell4.core.ParticleID object at 0x7f63db0939d8>, <ecell4.core.Particle object at 0x7f63db093930>), (<ecell4.core.ParticleID object at 0x7f63db0939f0>, <ecell4.core.Particle object at 0x7f63db093a68>), (<ecell4.core.ParticleID object at 0x7f63db093a38>, <ecell4.core.Particle object at 0x7f63db093ab0>)]


``Particle`` consists of ``species``, ``position``, ``radius`` and
``D``.

.. code:: python

    print([tuple(p.position()) for pid, p in w6.list_particles_exact(Species("A"))])


.. parsed-literal::

    [(0.729419355513528, 0.9559871961828321, 0.7664110218174756), (0.4843062642030418, 0.15541509888134897, 0.2066904455423355), (0.6487654014490545, 0.5129862017929554, 0.2914334700908512)]


.. code:: python

    pos = Real3(0.5, 0.5, 0.5)
    (pid1, p1), suc1 = w1.new_particle(sp1, pos)
    (pid2, p2), suc2 = w2.new_particle(sp1, pos)
    (pid3, p3), suc3 = w3.new_particle(sp1, pos)
    (pid4, p4), suc4 = w4.new_particle(sp1, pos)
    (pid5, p5), suc5 = w5.new_particle(sp1, pos)
    (pid6, p6), suc6 = w6.new_particle(sp1, pos)

.. code:: python

    # print(w1.get_particle(pid1))
    # print(w2.get_particle(pid2))
    print(w3.get_particle(pid3))
    print(w4.get_particle(pid4))
    # print(w5.get_particle(pid5))
    print(w6.get_particle(pid6))


.. parsed-literal::

    (<ecell4.core.ParticleID object at 0x7f63db093bd0>, <ecell4.core.Particle object at 0x7f63db093be8>)
    (<ecell4.core.ParticleID object at 0x7f63db093bd0>, <ecell4.core.Particle object at 0x7f63db093b70>)
    (<ecell4.core.ParticleID object at 0x7f63db093bd0>, <ecell4.core.Particle object at 0x7f63db0939d8>)


.. code:: python

    print(pid1, p1, suc1)
    print(pid2, p2, suc2)
    print(pid3, p3, suc3)
    print(pid4, p4, suc4)
    print(pid5, p5, suc5)
    print(pid6, p6, suc6)


.. parsed-literal::

    (<ecell4.core.ParticleID object at 0x7f63db093ab0>, <ecell4.core.Particle object at 0x7f63db093af8>, True)
    (<ecell4.core.ParticleID object at 0x7f63db093930>, <ecell4.core.Particle object at 0x7f63db093ac8>, True)
    (<ecell4.core.ParticleID object at 0x7f63db093a68>, <ecell4.core.Particle object at 0x7f63db093b28>, True)
    (<ecell4.core.ParticleID object at 0x7f63db0939f0>, <ecell4.core.Particle object at 0x7f63db093b58>, True)
    (<ecell4.core.ParticleID object at 0x7f63db093ae0>, <ecell4.core.Particle object at 0x7f63db093b88>, True)
    (<ecell4.core.ParticleID object at 0x7f63db093b40>, <ecell4.core.Particle object at 0x7f63db093bb8>, True)


