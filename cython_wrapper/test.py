import PySpecies
import PyCompartmentSpace
import PyODEWorld
import PyReactionRule
import PyNetworkModel

a = PySpecies.PySpecies("A")
world = PyODEWorld.PyOdeWorld(0.5)
world.add_species(a)

#print world.num_species()

def BuildEnvironmentTest():
    volume = 0.5
    sp1 = PySpecies.PySpecies("A")
    sp2 = PySpecies.PySpecies("B")

    # Build NetworkModel.
    rr = PyReactionRule.PyReactionRule()
    rr.add_reactant(sp1)
    rr.add_product(sp2)
    rr.set_k(0.5)

    rr_rev = PyReactionRule.PyReactionRule()
    rr_rev.add_reactant(sp2)
    rr_rev.add_product(sp1)
    rr_rev.set_k(0.2)

    model = PyNetworkModel.PyNetworkModel()
    model.add_species(sp1)
    model.add_species(sp2)
    model.add_reaction_rule(rr)
    model.add_reaction_rule(rr_rev)

    # Build Space
    volume = 0.5
    ode_world = PyODEWorld.PyOdeWorld(volume)
    ode_world.add_species(sp1)
    ode_world.add_species(sp2)
    ode_world.set_num_molecules(sp1, 0.8)
    ode_world.set_num_molecules(sp2, 0.2)
    print ode_world.num_species()


def CompartmentSpaceTest():
    volume = 0.5
    sp1 = PySpecies.PySpecies("A")
    sp2 = PySpecies.PySpecies("B")
    space= PyCompartmentSpace.PyCompartmentSpace(0.5)
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

#CompartmentSpaceTest()
BuildEnvironmentTest()
