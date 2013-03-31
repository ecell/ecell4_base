import PySpecies
import PyCompartmentSpace
import PyODEWorld

a = PySpecies.PySpecies("A")
world = PyODEWorld.PyOdeWorld(0.5)
world.add_species(a)

print world.num_species()

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

CompartmentSpaceTest()
