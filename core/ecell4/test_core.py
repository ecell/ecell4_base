from ecell4.core import *


# a = Species("A")
# space = CompartmentSpaceVectorImpl(0.5)
# space.add_species(a)


def BuildEnvironmentTest():
    volume = 0.5
    sp1 = Species("A")
    sp2 = Species("B")

    # Build NetworkModel.
    rr = ReactionRule()
    rr.add_reactant(sp1)
    rr.add_product(sp2)
    rr.set_k(0.5)

    rr_rev = ReactionRule()
    rr_rev.add_reactant(sp2)
    rr_rev.add_product(sp1)
    rr_rev.set_k(0.2)

    model = NetworkModel()
    model.add_species(sp1)
    model.add_species(sp2)
    model.add_reaction_rule(rr)
    model.add_reaction_rule(rr_rev)

    # Build Space
    volume = 0.5
    space = CompartmentSpaceVectorImpl(volume)
    space.add_species(sp1)
    space.add_species(sp2)
    space.add_molecules(sp1, 8)
    space.add_molecules(sp2, 5)
    print space.num_species()

def CompartmentSpaceTest():
    print 'start ...'
    volume = 0.5
    sp1 = Species("A")
    sp2 = Species("B")
    space = CompartmentSpaceVectorImpl(0.5)
    if space.volume() == 0.5:
        pass
    else:
        print "volume() fail"
    space.set_volume(0.8)

    if space.volume() == 0.8:
        pass
    else:
        print "set_volume() fail"

    space.add_species(sp1)
    space.add_species(sp2)
    if space.num_species() == 2:
        pass
    else:
        print "num_species() failed"

    space.add_molecules(sp1, 100)
    space.add_molecules(sp2, 200)
    if space.num_molecules(sp1) == 100 and space.num_molecules(sp2) == 200:
        pass
    else:
        print "num_molecules fail"

    space.remove_molecules(sp1, 50)
    if space.num_molecules(sp1) == 50:
        pass
    else:
        print "remove_molecules fail"
    space.remove_species(sp1)
    if space.num_species() == 1:
        pass
    else:
        print "remove_speces failed"
    print '... done'

CompartmentSpaceTest()
# BuildEnvironmentTest()
