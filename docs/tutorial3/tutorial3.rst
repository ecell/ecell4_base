
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
positions are treated as a ``Real3``. See also `8. More about 1. Brief
Tour of E-Cell4
Simulations <8.%20More%20about%201.%20Brief%20Tour%20of%20E-Cell4%20Simulations.ipynb>`__.

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

    w1.save("gillespie.h5")
    w2.save("ode.h5")
    w3.save("spatiocyte.h5")
    w4.save("bd.h5")
    w5.save("meso.h5")
    w6.save("egfrd.h5")
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

    w1.load("gillespie.h5")
    w2.load("ode.h5")
    w3.load("spatiocyte.h5")
    w4.load("bd.h5")
    w5.load("meso.h5")
    w6.load("egfrd.h5")
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

    print(GillespieWorld("gillespie.h5").t())
    print(ODEWorld("ode.h5").t())
    print(SpatiocyteWorld("spatiocyte.h5").t())
    print(BDWorld("bd.h5").t())
    print(MesoscopicWorld("meso.h5").t())
    print(EGFRDWorld("egfrd.h5").t())


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
too.

.. code:: python

    w1 = GillespieWorld()
    w2 = ODEWorld()
    w3 = SpatiocyteWorld()
    w4 = BDWorld()
    w5 = MesoscopicWorld()
    w6 = EGFRDWorld()

First, you can place a molecule at the certain position with
``new_particle``.

.. code:: python

    sp1 = Species("A", "0.0025", "1")
    pos = Real3(0.5, 0.5, 0.5)
    (pid1, p1), suc1 = w1.new_particle(sp1, pos)
    (pid2, p2), suc2 = w2.new_particle(sp1, pos)
    (pid3, p3), suc3 = w3.new_particle(sp1, pos)
    (pid4, p4), suc4 = w4.new_particle(sp1, pos)
    (pid5, p5), suc5 = w5.new_particle(sp1, pos)
    (pid6, p6), suc6 = w6.new_particle(sp1, pos)

``new_particle`` returns a particle created and whether it's succeeded
or not. However the resolution in representation of molecules differs.
For example, ``GillespieWorld`` has almost no information about the
coordinate of molecules. Thus, it simply ignores the given position, and
just counts up the number of molecules here.

``ParticleID`` is a pair of ``Integer``\ s named ``lot`` and ``serial``.

.. code:: python

    print(pid6.lot(), pid6.serial())
    print(pid6 == ParticleID((0, 1)))


.. parsed-literal::

    (0, 1L)
    True


Particle simulators, i.e. ``spatiocyte``, ``bd`` and ``egfrd``, provide
an interface to access a particle by its id. ``has_particle`` returns if
a particles exists or not for the given ``ParticleID``.

.. code:: python

    # print(w1.has_particle(pid1))
    # print(w2.has_particle(pid2))
    print(w3.has_particle(pid3))  # => True
    print(w4.has_particle(pid4))  # => True
    # print(w5.has_particle(pid5))
    print(w6.has_particle(pid6))  # => True


.. parsed-literal::

    True
    True
    True


After checking the existency, you can get the partcle by
``get_particle`` as follows.

.. code:: python

    # pid1, p1 = w1.get_particle(pid1)
    # pid2, p2 = w2.get_particle(pid2)
    pid3, p3 = w3.get_particle(pid3)
    pid4, p4 = w4.get_particle(pid4)
    # pid5, p5 = w5.get_particle(pid5)
    pid6, p6 = w6.get_particle(pid6)

``Particle`` consists of ``species``, ``position``, ``radius`` and
``D``.

.. code:: python

    # print(p1.species().serial(), tuple(p1.position()), p1.radius(), p1.D())
    # print(p2.species().serial(), tuple(p2.position()), p2.radius(), p2.D())
    print(p3.species().serial(), tuple(p3.position()), p3.radius(), p3.D())
    print(p4.species().serial(), tuple(p4.position()), p4.radius(), p4.D())
    # print(p5.species().serial(), tuple(p5.position()), p5.radius(), p5.D())
    print(p6.species().serial(), tuple(p6.position()), p6.radius(), p6.D())


.. parsed-literal::

    (u'A', (0.5062278801751902, 0.5080682368868706, 0.5), 0.0025, 1.0)
    (u'A', (0.5, 0.5, 0.5), 0.0025, 1.0)
    (u'A', (0.5, 0.5, 0.5), 0.0025, 1.0)


In the case of ``spatiocyte``, a particle position is automatically
round to the center of the voxel nearest to the given position.

You can even move the position of the particle. ``update_particle``
replace the particle specified with the given ``ParticleID`` with the
given ``Particle`` and return ``False``. If no corresponding particle is
found, create new particle and return ``True``. If you give a
``Particle`` with the different type of ``Species``, the ``Species`` of
the ``Particle`` will be also changed.

.. code:: python

    newp = Particle(sp1, Real3(0.3, 0.3, 0.3), 0.0025, 1)
    # print(w1.update_particle(pid1, newp))
    # print(w2.update_particle(pid2, newp))
    print(w3.update_particle(pid3, newp))
    print(w4.update_particle(pid4, newp))
    # print(w5.update_particle(pid5, newp))
    print(w6.update_particle(pid6, newp))


.. parsed-literal::

    False
    False
    False


``list_particles`` and ``list_particles_exact`` return a list of pairs
of ``ParticleID`` and ``Particle`` in the ``World``. ``World``
automatically makes up for the gap with random numbers. For example,
``GillespieWorld`` returns a list of positions randomly distributed in
the ``World`` size.

.. code:: python

    print(w1.list_particles_exact(sp1))
    # print(w2.list_particles_exact(sp1))  # ODEWorld has no member named list_particles
    print(w3.list_particles_exact(sp1))
    print(w4.list_particles_exact(sp1))
    print(w5.list_particles_exact(sp1))
    print(w6.list_particles_exact(sp1))


.. parsed-literal::

    [(<ecell4.core.ParticleID object at 0x7f29138bfc48>, <ecell4.core.Particle object at 0x7f29138bfa80>)]
    [(<ecell4.core.ParticleID object at 0x7f29138bfc48>, <ecell4.core.Particle object at 0x7f29138bfc90>)]
    [(<ecell4.core.ParticleID object at 0x7f29138bfc48>, <ecell4.core.Particle object at 0x7f29138bfa38>)]
    [(<ecell4.core.ParticleID object at 0x7f29138bfc48>, <ecell4.core.Particle object at 0x7f29138bfa80>)]
    [(<ecell4.core.ParticleID object at 0x7f29138bfc48>, <ecell4.core.Particle object at 0x7f29138bfc90>)]


You can remove a specific particle with ``remove_particle``.

.. code:: python

    # w1.remove_particle(pid1)
    # w2.remove_particle(pid2)
    w3.remove_particle(pid3)
    w4.remove_particle(pid4)
    # w5.remove_particle(pid5)
    w6.remove_particle(pid6)
    # print(w1.has_particle(pid1))
    # print(w2.has_particle(pid2))
    print(w3.has_particle(pid3))  # => False
    print(w4.has_particle(pid4))  # => False
    # print(w5.has_particle(pid5))
    print(w6.has_particle(pid6))  # => False


.. parsed-literal::

    False
    False
    False


3.3. Lattice-based Coordinate
-----------------------------

In addition to the common interface, each ``World`` can have their own
interfaces. As an example, we explain methods to handle lattice-based
coordinate here. ``SpatiocyteWorld`` is based on a space discretized to
hexiagonal close packing lattices, ``LatticeSpace``.

.. code:: python

    w = SpatiocyteWorld(Real3(1, 2, 3), voxel_radius=0.01)
    w.bind_to(m)

The size of a single lattice, called ``Voxel``, can be obtained by
``voxel_radius()``. ``SpatiocyteWorld`` has methods to get the numbers
of rows, columns, and layers. These sizes are automatically calculated
based on the given ``edge_lengths`` at the construction.

.. code:: python

    print(w.voxel_radius())  # => 0.01
    print(tuple(w.shape()))  # => (62, 152, 116)
    print(w.col_size(), w.row_size(), w.layer_size())  # => (62, 152, 116)
    print(w.size())  # => 1093184 = 62 * 152 * 116


.. parsed-literal::

    0.01
    (62, 152, 116)
    (62, 152, 116)
    1093184


A position in the lattice-based space is treated as an ``Integer3``,
column, row and layer, called a global coordinate. Thus,
``SpatiocyteWorld`` provides the function to convert the ``Real3`` into
a lattice-based coordinate.

.. code:: python

    p1 = Real3(0.5, 0.5, 0.5)
    g1 = w.position2global(p1)
    p2 = w.global2position(g1)
    print(tuple(g1))  # => (31, 25, 29)
    print(tuple(p2))  # => (0.5062278801751902, 0.5080682368868706, 0.5)


.. parsed-literal::

    (31, 25, 29)
    (0.5062278801751902, 0.5080682368868706, 0.5)


In ``SpatiocyteWorld``, the global coordinate is translated to a single
integer. It is just called a coordinate. You can also treat the
coordinate as in the same way with a global coordinate.

.. code:: python

    p1 = Real3(0.5, 0.5, 0.5)
    c1 = w.position2coordinate(p1)
    p2 = w.coordinate2position(c1)
    g1 = w.coord2global(c1)
    print(c1)  # => 278033
    print(tuple(p2))  # => (0.5062278801751902, 0.5080682368868706, 0.5)
    print(tuple(g1))  # => (31, 25, 29)


.. parsed-literal::

    278033
    (0.5062278801751902, 0.5080682368868706, 0.5)
    (31, 25, 29)


With these coordinates, you can handle a ``Voxel``, which represents a
``Particle`` object. Instead of ``new_particle``, ``new_voxel`` provides
the way to create a new ``Voxel`` with a coordinate.

.. code:: python

    c1 = w.position2coordinate(Real3(0.5, 0.5, 0.5))
    ((pid, v), is_succeeded) = w.new_voxel(Species("A"), c1)
    print(pid, v, is_succeeded)


.. parsed-literal::

    (<ecell4.core.ParticleID object at 0x7f29138bfa80>, <ecell4.core.Voxel object at 0x7f29138bfcd8>, True)


A ``Voxel`` consists of ``species``, ``coordinate``, ``radius`` and
``D``.

.. code:: python

    print(v.species().serial(), v.coordinate(), v.radius(), v.D())  # => (u'A', 278033, 0.0025, 1.0)


.. parsed-literal::

    (u'A', 278033, 0.0025, 1.0)


Of course, you can get a voxel and list voxels with ``get_voxel`` and
``list_voxels_exact`` similar to ``get_particle`` and
``list_particles_exact``.

.. code:: python

    print(w.num_voxels_exact(Species("A")))
    print(w.list_voxels_exact(Species("A")))
    print(w.get_voxel(pid))


.. parsed-literal::

    1
    [(<ecell4.core.ParticleID object at 0x7f29138bfcc0>, <ecell4.core.Voxel object at 0x7f29138bfcf0>)]
    (<ecell4.core.ParticleID object at 0x7f29138bfcc0>, <ecell4.core.Voxel object at 0x7f29138bfbb8>)


You can move and update the voxel with ``update_voxel`` corresponding to
``update_particle``.

.. code:: python

    c2 = w.position2coordinate(Real3(0.5, 0.5, 1.0))
    w.update_voxel(pid, Voxel(v.species(), c2, v.radius(), v.D()))
    pid, newv = w.get_voxel(pid)
    print(c2)  # => 278058
    print(newv.species().serial(), newv.coordinate(), newv.radius(), newv.D())  # => (u'A', 278058, 0.0025, 1.0)
    print(w.num_voxels_exact(Species("A")))  # => 1


.. parsed-literal::

    278058
    (u'A', 278058, 0.0025, 1.0)
    1


Finally, ``remove_voxel`` remove a voxel as ``remove_particle`` does.

.. code:: python

    print(w.has_voxel(pid))  # => True
    w.remove_voxel(pid)
    print(w.has_voxel(pid))  # => False


.. parsed-literal::

    True
    False


3.4 Structure
-------------

.. code:: python

    w1 = GillespieWorld()
    w2 = ODEWorld()
    w3 = SpatiocyteWorld()
    w4 = BDWorld()
    w5 = MesoscopicWorld()
    w6 = EGFRDWorld()

By using a ``Shape`` object, you can confine initial positions of
molecules to a part of ``World``. In the case below, 60 molecules are
positioned inside the given ``Sphere``. Diffusion of the molecules
placed here is **NOT** restricted in the ``Shape``. This ``Shape`` is
only for the initialization.

.. code:: python

    sp1 = Species("A", "0.0025", "1")
    sphere = Sphere(Real3(0.5, 0.5, 0.5), 0.3)
    w1.add_molecules(sp1, 60, sphere)
    w2.add_molecules(sp1, 60, sphere)
    w3.add_molecules(sp1, 60, sphere)
    w4.add_molecules(sp1, 60, sphere)
    w5.add_molecules(sp1, 60, sphere)
    w6.add_molecules(sp1, 60, sphere)

A property of ``Species``, ``'location'``, is available to restrict
diffusion of molecules. ``'location'`` is not fully supported yet, but
only supported in ``spatiocyte`` and ``meso``. ``add_structure`` defines
a new structure given as a pair of ``Species`` and ``Shape``.

.. code:: python

    membrane = SphericalSurface(Real3(0.5, 0.5, 0.5), 0.4)  # This is equivalent to call `Sphere(Real3(0.5, 0.5, 0.5), 0.4).surface()`
    w3.add_structure(Species("M"), membrane)
    w5.add_structure(Species("M"), membrane)

After defining a structure, you can bind molecules to the structure as
follows:

.. code:: python

    sp2 = Species("B", "0.0025", "0.1", "M")  # `'location'` is the fourth argument
    w3.add_molecules(sp2, 60)
    w5.add_molecules(sp2, 60)

The molecules bound to a ``Species`` named ``B`` diffuse on a structure
named ``M``, which has a shape of ``SphericalSurface`` (a hollow
sphere). In ``spatiocyte``, a structure is represented as a set of
particles with ``Species`` ``M`` occupying a voxel. It means that
molecules not belonging to the structure is not able to overlap the
voxel and it causes a collision. On the other hand, in ``meso``, a
structure means a list of subvolumes. Thus, a structure doesn't avoid an
incursion of other particles.

3.5. Random Number Generator
----------------------------

A random number generator is also a part of ``World``. All ``World``
except ``ODEWorld`` store a random number generator, and updates it when
the simulation needs a random value. On E-Cell4, only one class
``GSLRandomNumberGenerator`` is implemented as a random number
generator.

.. code:: python

    rng1 = GSLRandomNumberGenerator()
    print([rng1.uniform_int(1, 6) for _ in range(20)])


.. parsed-literal::

    [6, 1, 2, 6, 2, 3, 6, 5, 4, 5, 5, 4, 2, 5, 4, 2, 3, 3, 2, 2]


With no argument, the random number generator is always initialized with
a seed, ``0``.

.. code:: python

    rng2 = GSLRandomNumberGenerator()
    print([rng2.uniform_int(1, 6) for _ in range(20)])  # => same as above


.. parsed-literal::

    [6, 1, 2, 6, 2, 3, 6, 5, 4, 5, 5, 4, 2, 5, 4, 2, 3, 3, 2, 2]


You can initialize the seed with an integer as follows:

.. code:: python

    rng2 = GSLRandomNumberGenerator()
    rng2.seed(15)
    print([rng2.uniform_int(1, 6) for _ in range(20)])


.. parsed-literal::

    [6, 5, 2, 4, 1, 1, 3, 5, 2, 6, 4, 1, 2, 5, 2, 5, 1, 2, 2, 6]


When you call the ``seed`` function with no input, the seed is drawn
from the current time.

.. code:: python

    rng2 = GSLRandomNumberGenerator()
    rng2.seed()
    print([rng2.uniform_int(1, 6) for _ in range(20)])


.. parsed-literal::

    [6, 6, 2, 2, 4, 2, 6, 4, 6, 5, 6, 2, 4, 3, 1, 4, 6, 4, 2, 3]


``GSLRandomNumberGenerator`` provides several ways to get a random
number.

.. code:: python

    print(rng1.uniform(0.0, 1.0))
    print(rng1.uniform_int(0, 100))
    print(rng1.gaussian(1.0))


.. parsed-literal::

    0.0303352042101
    33
    0.893555545521


``World`` accepts a random number generator at the construction. As a
default, ``GSLRandomNumberGenerator()`` is used. Thus, when you don't
give a generator, behavior of the simulation is always same
(determinisitc).

.. code:: python

    rng = GSLRandomNumberGenerator()
    rng.seed()
    w1 = GillespieWorld(Real3(1, 1, 1), rng=rng)

You can access the ``GSLRandomNumberGenerator`` in a ``World`` through
``rng`` function.

.. code:: python

    print(w1.rng().uniform(0.0, 1.0))


.. parsed-literal::

    0.913418973563


``rng()`` returns a shared pointer to the ``GSLRandomNumberGenerator``.
Thus, in the example above, ``rng`` and ``w1.rng()`` point exactly the
same thing.
