
Tutorial 1 (Species Basics)
===========================

This is a tutorial for E-Cell4. Here, we introduce how to use Species.

Import the core library, first:

.. code:: python

    from ecell4.core import *
This enables Species class.

.. code:: python

    print Species("A") == Species("B"), Species("A") == Species("A")

.. parsed-literal::

    False True


A Species can be translated into an unique string, named "serial":

.. code:: python

    print Species("A").serial(), Species("B").serial()

.. parsed-literal::

    A B


A Species can have attributes as a pair of strings.

.. code:: python

    sp = Species("A")
    sp.set_attribute("radius", "0.0025")
    sp.set_attribute("D", "1.0")
    print sp.has_attribute("radius"), sp.has_attribute("spam")
    print sp.get_attribute("radius"), sp.get_attribute("D")
    sp.remove_attribute("radius")
    print sp.has_attribute("radius")
    del sp

.. parsed-literal::

    True False
    0.0025 1.0
    False


Tips: Especially for "radius" and "D", A Species accepts these
attributes at the instantiation:

.. code:: python

    sp = Species("A", "0.0025", "1")
    print sp.get_attribute("radius"), sp.get_attribute("D")
    del sp

.. parsed-literal::

    0.0025 1


A Species consists of one or more UnitSpecies. UnitSpecies are
automatically sorted in the Species.

.. code:: python

    sp = Species()
    usp = UnitSpecies("C")
    print usp.serial()
    sp.add_unit(usp)
    sp.add_unit(UnitSpecies("A"))
    sp.add_unit(UnitSpecies("B"))
    print sp.serial(), sp.num_units()
    del usp, sp

.. parsed-literal::

    C
    A.B.C 3


A Species can be reproduced from serial. In the serial, all serials of
UnitSpecies are joined with the separator, ".".

.. code:: python

    sp = Species("C.A.B")
    print sp.serial()
    print Species("A.B.C") == Species("C.A.B")
    del sp

.. parsed-literal::

    A.B.C
    True


An UnitSpecies can have sites. Sites are also sorted automatically in
the UnitSpecies.

.. code:: python

    usp = UnitSpecies("A")
    usp.add_site("us", "u", "")
    usp.add_site("ps", "p", "_")
    usp.add_site("bs", "", "_")
    print usp.name(), usp.serial()
    del usp

.. parsed-literal::

    A A(bs^_,ps=p^_,us=u)


An UnitSpecies can be reproduced from its serial.

.. code:: python

    usp = UnitSpecies()
    usp.deserialize("A(bs^_, us = u, ps = p^_)")
    print usp.serial()
    del usp

.. parsed-literal::

    A(bs^_,ps=p^_,us=u)


Of course, a site of UnitSpecies is available even in a Species.

.. code:: python

    sp = Species("A(bs^1, ps=u).A(bs, ps=p^1)")
    print sp.serial(), sp.num_units()
    del sp

.. parsed-literal::

    A(bs,ps=p^1).A(bs^1,ps=u) 2


