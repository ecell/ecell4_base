
Tutorial 6 (World Advanced)
===========================

This is a tutorial for E-Cell4.

.. code:: python

    from ecell4.core import *
    from ecell4.gillespie import GillespieWorld
    from ecell4.ode import ODEWorld
    from ecell4.lattice import LatticeWorld
    from ecell4.bd import BDWorld
    
    edge_lengths = Real3(1, 2, 3)
.. code:: python

    m = NetworkModel()
    m.add_species_attribute(Species("A", "0.0025", "1"))
``ParticleSpace`` interfaces:

.. code:: python

    w = BDWorld(edge_lengths)
    # w = LatticeWorld(edge_lengths)
    w.bind_to(m)
Create a new particle. ``new_particle`` returns a pair of ``ParticleID``
and ``Particle`` with if the creation is succeeded or not.

.. code:: python

    ((pid, p), is_succeeded) = w.new_particle(Species("A"), Real3(0.5, 0.5, 0.5))
.. code:: python

    print pid, p, is_succeeded
    print w.num_molecules(Species("A"))

.. parsed-literal::

    <ecell4.core.ParticleID object at 0x7f114e57ea68> <ecell4.core.Particle object at 0x7f114e57e930> True
    1


``ParticleID`` is a pair of ``Integer``\ s named ``serial`` and ``lot``.
``Particle`` consists of ``species``, ``position``, ``radius``, ``D``.

.. code:: python

    print pid.serial(), pid.lot()
    print p.species().serial(), tuple(p.position()), p.radius(), p.D()

.. parsed-literal::

    1 0
    A (0.5, 0.5, 0.5) 0.0025 1.0


``ParticleSpace`` has interfaces to check if the ``ParticleID`` exists
(``has_particle``), to list all particles belonging to the ``Species``
(``list_particles``), and to get ``Particle`` assigned to the given
``ParticleID`` (``get_particle``).

.. code:: python

    print w.has_particle(pid)
    print w.list_particles(Species("A"))
    print w.get_particle(pid)

.. parsed-literal::

    True
    [(<ecell4.core.ParticleID object at 0x7f114e57e9d8>, <ecell4.core.Particle object at 0x7f114e57e978>)]
    (<ecell4.core.ParticleID object at 0x7f114e57e9d8>, <ecell4.core.Particle object at 0x7f114e57ea38>)


To move a ``Particle`` to the new position, use ``update_particle``.
``update_particle`` tries to replace the ``Particle`` with given one.

.. code:: python

    w.update_particle(pid, Particle(p.species(), Real3(0.5, 1.5, 2.5), p.radius(), p.D()))
    _, newp = w.get_particle(pid)
    print tuple(newp.position()), w.num_molecules(Species("A"))

.. parsed-literal::

    (0.5, 1.5, 2.5) 1


You can remove a ``Particle`` by ``remove_particle``.

.. code:: python

    w.remove_particle(pid)
    print w.has_particle(pid)

.. parsed-literal::

    False


.. code:: python

    del w
``LatticeSpace`` interfaces:

.. code:: python

    w = LatticeWorld(edge_lengths)
    w.bind_to(m)
``LatticeSpace`` has interfaces to give essential information about it.

.. code:: python

    print w.voxel_radius(), w.row_size(), w.col_size(), w.layer_size()

.. parsed-literal::

    0.01 150 61 115


Positions in ``LatticeSpace`` is represented as a single ``Integer``
named ``coordinate``, which is corresponding to ``Real3``
(``position``) in ``ParticleSpace``. ``coordinate2position`` and
``position2coordinate`` give the way to convert between them.

.. code:: python

    coord = w.position2coordinate(Real3(0.5, 0.5, 0.5))
    pos = w.coordinate2position(coord)
    new_coord = w.position2coordinate(pos)
    
    print coord, tuple(pos)
    print new_coord, tuple(w.coordinate2position(new_coord))
    print w.position2coordinate(Real3(0.5, 1.5, 2.5))

.. parsed-literal::

    260725 (0.48989794855663565, 0.48497422611928565, 0.5)
    260725 (0.48989794855663565, 0.48497422611928565, 0.5)
    791525


Interfaces similar to ``ParticleSpace`` are available in
``LatticeSpace``.

.. code:: python

    ((pid, v), is_succeeded) = w.new_voxel(Species("A"), coord)
.. code:: python

    print pid, v, is_succeeded
    print w.num_molecules(Species("A"))

.. parsed-literal::

    <ecell4.core.ParticleID object at 0x7f114e57e978> <ecell4.core.Voxel object at 0x7f114e57ea20> True
    1


.. code:: python

    print v.species().serial(), v.coordinate(), v.radius(), v.D()

.. parsed-literal::

    A 260725 0.0025 1.0


.. code:: python

    print w.has_particle(pid)
    print w.list_voxels(Species("A"))
    print w.get_voxel(pid)

.. parsed-literal::

    True
    [(<ecell4.core.ParticleID object at 0x7f114e57ea68>, <ecell4.core.Voxel object at 0x7f114e57ea38>)]
    (<ecell4.core.ParticleID object at 0x7f114e57ea68>, <ecell4.core.Voxel object at 0x7f114e57ea80>)


.. code:: python

    w.update_voxel(pid, Voxel(v.species(), 791525, v.radius(), v.D()))
    _, newv = w.get_voxel(pid)
    print newv.coordinate(), w.num_molecules(Species("A"))

.. parsed-literal::

    791525 1


.. code:: python

    w.remove_voxel(pid)
    print w.has_voxel(pid)

.. parsed-literal::

    False


``RandomNumberGenerator`` interfaces:

.. code:: python

    rng1, rng2 = GSLRandomNumberGenerator(), GSLRandomNumberGenerator()
    rng1.seed(0)
    rng2.seed(0)
    print rng1.uniform(0, 1), rng2.uniform(0, 1)
    print rng1.uniform_int(0, 100), rng2.uniform_int(0, 100)
    print rng1.gaussian(0.0, 1.0), rng2.gaussian(0.0, 1.0)

.. parsed-literal::

    0.999741748907 0.999741748907
    16 16
    0.133918608119 0.133918608119


``GillespieWorld``, ``LatticeWorld`` and ``BDWorld`` can be constructed
with a ``GSLRandomNumberGenerator``. ``rng()`` returns a shared object
of ``RandomNumberGenerator``. In below, ``w.rng()`` and ``rng1`` point
the same ``RandomNumberGenerator``.

.. code:: python

    w = GillespieWorld(Real3(1, 1, 1), rng1)
    # w = BDWorld(Real3(1, 1, 1), rng1)
    # w = LatticeWorld(Real3(1, 1, 1), 0.05, rng1) # The second argument is voxel_radius.
.. code:: python

    print w.rng().uniform(0, 1), rng2.uniform(0, 1), rng1.uniform(0, 1)

.. parsed-literal::

    0.231656542746 0.231656542746 0.484973614337


.. code:: python

    del rng1, rng2
    print w.rng().uniform(0, 1)
    del w

.. parsed-literal::

    0.957476956537


