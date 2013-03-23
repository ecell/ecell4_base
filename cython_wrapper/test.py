import PySpecies
import PyODEWorld

a = PySpecies.PySpecies("A")
world = PyODEWorld.PyOdeWorld(0.5)
world.add_species(a)

print world.num_species()
