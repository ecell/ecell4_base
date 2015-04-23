
Tutorial 7 (Model Advanced)
===========================

This is a tutorial for E-Cell4.

.. code:: python

    %matplotlib inline
.. code:: python

    from ecell4.core import *
    from ecell4.reaction_reader.decorator import species_attributes, reaction_rules
.. code:: python

    @species_attributes
    def attrgen(radius, D):
        K | {"radius": radius, "D": D}
        Kp | {"radius": radius, "D": D}
        Kpp | {"radius": radius, "D": D}
        KK | {"radius": radius, "D": D}
        PP | {"radius": radius, "D": D}
        K.KK | {"radius": radius, "D": D}
        Kp.KK | {"radius": radius, "D": D}
        Kpp.PP | {"radius": radius, "D": D}
        Kp.PP | {"radius": radius, "D": D}
    
    @reaction_rules
    def rulegen(kon1, koff1, kcat1, kon2, koff2, kcat2):
        (K + KK == K.KK | (kon1, koff1)
            > Kp + KK | kcat1
            == Kp.KK | (kon2, koff2)
            > Kpp + KK | kcat2)
    
        (Kpp + PP == Kpp.PP | (kon1, koff1)
            > Kp + PP | kcat1
            == Kp.PP | (kon2, koff2)
            > K + PP | kcat2)
.. code:: python

    m = NetworkModel()
.. code:: python

    for i, sp in enumerate(attrgen("0.0025", "1")):
        print i, sp.serial(), sp.get_attribute("radius"), sp.get_attribute("D")
        m.add_species_attribute(sp)

.. parsed-literal::

    0 K 0.0025 1
    1 Kp 0.0025 1
    2 Kpp 0.0025 1
    3 KK 0.0025 1
    4 PP 0.0025 1
    5 K.KK 0.0025 1
    6 KK.Kp 0.0025 1
    7 Kpp.PP 0.0025 1
    8 Kp.PP 0.0025 1


.. code:: python

    ka1, kd1, kcat1 = 0.04483455086786913, 1.35, 1.5
    ka2, kd2, kcat2 = 0.09299017957780264, 1.73, 15.0
    
    for i, rr in enumerate(rulegen(ka1, kd2, kcat1, ka2, kd2, kcat2)):
        reactants, products, k = rr.reactants(), rr.products(), rr.k()
        print i, map(lambda sp: sp.serial(), reactants), map(lambda sp: sp.serial(), products), k
        m.add_reaction_rule(rr)

.. parsed-literal::

    0 ['K', 'KK'] ['K.KK'] 0.0448345508679
    1 ['K.KK'] ['K', 'KK'] 1.73
    2 ['K.KK'] ['KK', 'Kp'] 1.5
    3 ['KK', 'Kp'] ['KK.Kp'] 0.0929901795778
    4 ['KK.Kp'] ['KK', 'Kp'] 1.73
    5 ['KK.Kp'] ['KK', 'Kpp'] 15.0
    6 ['Kpp', 'PP'] ['Kpp.PP'] 0.0448345508679
    7 ['Kpp.PP'] ['Kpp', 'PP'] 1.73
    8 ['Kpp.PP'] ['Kp', 'PP'] 1.5
    9 ['Kp', 'PP'] ['Kp.PP'] 0.0929901795778
    10 ['Kp.PP'] ['Kp', 'PP'] 1.73
    11 ['Kp.PP'] ['K', 'PP'] 15.0


.. code:: python

    from ecell4.gillespie import GillespieWorld as world_type, GillespieSimulator as simulator_type
    # from ecell4.ode import ODEWorld as world_type, ODESimulator as simulator_type
    
    w = world_type(Real3(1, 1, 1))
    w.bind_to(m)
    w.add_molecules(Species("K"), 120)
    w.add_molecules(Species("KK"), 30)
    w.add_molecules(Species("PP"), 30)
    sim = simulator_type(m, w)
.. code:: python

    next_time, dt = 0.0, 1.0
    data = [(w.t(),
        w.num_molecules(Species("K")) + w.num_molecules(Species("K.KK")),
        w.num_molecules(Species("Kp")) + w.num_molecules(Species("Kp.KK")) + w.num_molecules(Species("Kp.PP")),
        w.num_molecules(Species("Kpp")) + w.num_molecules(Species("Kpp.PP")))]
    for i in range(60):
        next_time += dt
        while (sim.step(next_time)): pass
        data.append((w.t(),
            w.num_molecules(Species("K")) + w.num_molecules(Species("K.KK")),
            w.num_molecules(Species("Kp")) + w.num_molecules(Species("Kp.KK")) + w.num_molecules(Species("Kp.PP")),
            w.num_molecules(Species("Kpp")) + w.num_molecules(Species("Kpp.PP"))))
.. code:: python

    import matplotlib.pylab as plt
    from numpy import array
    
    data = array(data)
    plt.plot(data.T[0], data.T[1], "r-", label="K")
    plt.plot(data.T[0], data.T[2], "g--", label="Kp")
    plt.plot(data.T[0], data.T[3], "b:", label="Kpp")
    plt.xlabel("Time")
    plt.ylabel("Number Of Molecules")
    plt.xlim(data.T[0][0], data.T[0][-1])
    plt.legend(loc="best", shadow=True)
    plt.show()


.. image:: tutorial7_files/tutorial7_9_0.png

