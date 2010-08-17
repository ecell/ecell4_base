import model
from _gfrd import create_planar_surface, create_cylindrical_surface, \
                  create_cuboidal_region, ReactionRule
from gfrdbase import throw_in_particles, place_particle

__all__ = [
    'throw_in',
    'place',
    'ParticleType',
    'EGFRDModel',
    ]

def throw_in(simulator, type, n):
    """Add n particles of a certain type to the simulator.

    Return nothing.

    simulator -- a EGFRDSimulator.
    type      -- a (ParticleType, Region)-tuple, or
                 a (ParticleType, Surface)-tuple, or
                 a ParticleType, for a ParticleType in the "world".

    Make sure to first add this ParticleType to the model with a call to 
    EGFRDModel.add_particle_type(), together with the mentioned Region 
    or Surface.

    Use this function together with an EGFRDModel, instead of 
    throw_in_particles from gfrdbase.py. 

    """
    species_type = simulator.world.model._get_species_type(type)
    throw_in_particles(simulator.world, species_type.id, n)

def place(simulator, type, position):
    """Place a particle of a certain type at a specific position in the 
    simulator.

    Return nothing.

    simulator -- a EGFRDSimulator.
    type      -- a (ParticleType, Region)-tuple, or
              -- a (ParticleType, Surface)-tuple, or
                 a ParticleType, for particles in the "world".

    Make sure to first add this ParticleType to the model with a call to 
    EGFRDModel.add_particle_type(), together with the mentioned Region 
    or Surface.

    Use this function together with an EGFRDModel, instead of 
    place_particle from gfrdbase.py. 

    """
    species_type = simulator.world.model._get_species_type(type)
    place_particle(s.world, species_type.id, position)


class ParticleType(object):
    def __init__(self, name, D=None, radius=None):
        """Define a new type of particles.

        name   -- the name of this ParticleType. Should not be specific 
                  to any Region or Surface.
        D      -- the default diffusion constant for this ParticleType.
                  Optional, but then you have to define a diffusion 
                  constant for this ParticleType in a specific Region 
                  or Surface in the call to EGFRDModel.add_particle_type().
        radius -- default radius for this ParticleType.
                  Optional, but then you have to define a radius for 
                  this ParticleType in a specific Region or Surface in 
                  the call to EGFRDModel.add_particle_type().

        """
        self.name = name
        self.D = D
        self.radius = radius


class EGFRDModel(model.ParticleModel):
    """Define the type of particles and reaction rules for an eGFRD 
    simulation.
    
    """
    def add_planar_surface(self, id, corner, unit_x, unit_y):
        """Create a new planar surface and add it to the model.

        Return the created surface.

        id     -- a descriptive name.
        corner -- the point [x, y, z] on the plane closest to [0, 0, 0].
        unit_x -- a unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along 
                  the plane.
        unit_y -- a unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along 
                  the plane and perpendicular to unit_x.
        
        Surfaces are not allowed to touch or overlap.

        """
        surface = create_planar_surface(id, corner, unit_x, unit_y,
                                        self.world_size, self.world_size)
        return self.add_structure(surface)

    def add_cylindrical_surface(self, id, corner, radius, orientation):
        """Create a new cylindrical surface and add it to the model.

        Return the created surface.

        id          -- a descriptive name.
        corner      -- the point [x, y, z] on the axis of the cylinder
                       closest to [0, 0, 0].
        radius      -- the radius of the cylinder.
        orientation -- the unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1]
                       along the axis of the cylinder.

        Surfaces are not allowed to touch or overlap.

        """
        surface = create_cylindrical_surface(id, corner, radius, 
                                             orientation, self.world_size)
        return self.add_structure(surface)

    def add_cuboidal_region(self, id, corner, diagonal):
        """Create a new cuboidal region and add it to the model.

        Return the created region.

        id       -- a descriptive name.
        corner   -- the point [x, y, z] of the cuboidal region closest 
                    to [0, 0, 0].
        diagonal -- the vector [x, y, z] from the corner closest to
                    [0, 0, 0], to the corner furthest away from
                    [0, 0, 0].

        """
        region = create_cuboidal_region(id, corner, diagonal)
        return self.add_structure(region)

    def add_particle_type(self, particle_type, structure=None, v=None, 
                          D=None, radius=None):
        """Add a ParticleType to the model (in a specific Region or 
        Surface).

        Return nothing.
        
        particle_type -- a ParticleType.
        structure     -- the Region or Surface in which this 
                         ParticleType can exist.
                         Optional. If you do not specify a Region or 
                         Surface the ParticleType is added to the "world".
        v             -- the drift term for this ParticleType on a 
                         CylindricalSurface (1D drift).
                         Optional.
        D             -- the diffusion constant for this ParticleType in 
                         this Region or Surface, if different from the 
                         default diffusion constant for this 
                         ParticleType.
                         Optional, but then you have to define a default 
                         diffusion constant for this ParticleType in 
                         it's constructor call.
        radius        -- the radius for this ParticleType in this 
                         Region or Surface, if different from the 
                         default radius for this ParticleType.
                         Optional, but then you have to define a default 
                         radius for this ParticleType in it's 
                         constructor call.

        Note: if a ParticleType should be able to exist in the "world" 
        as well as on one of the previously added Surfaces or Regions, 
        then it should be added twice. Ones with and ones without an 
        explicit structure argument.

        Internally this method creates a new _gfrd.SpeciesType and adds 
        it to the model.

        """
        if structure == None:
            structure = self.get_structure("world")
        else:
            assert self.get_structure(structure.id) == structure, \
                   ('Region or Surface %s does not exist. Use '
                    'add_planar_surface or add_cylindrical_surface.' % 
                    structure)

        if D == None:
            assert particle_type.D != None, \
                   ('Diffusion constant of ParticleType %s not specified.' % 
                    particle_type)
            D = particle_type.D

        if radius == None:
            assert particle_type.radius != None, \
                   'Radius of ParticleType %s not specified.' % particle_type
            radius = particle_type.radius

        name = self._get_species_type_name(particle_type, structure)
        species_type = model.Species(name, D, radius, structure.id)

        self.add_species_type(species_type)

    def add_reaction_rule(self, reactants, products, k):
        """Create a new reaction rule and add it to the model.

        Return the created reaction rule.

        ### Unimolecular and bimolecular reaction rules.
        # Syntax.
        add_reaction_rule([reactants], [products], rate)

        where each reactant and product is of type:
            (ParticleType, Region), or
            (ParticleType, Surface), or
            ParticleType, for a ParticleType in the "world".

        and all reactants and products should be in the same Region, 
        Surface or "world".

        A unimolecular reaction rule (a reaction rule with 1 reactant) 
        can have 0 (empty list), 1 or 2 products.
        A bimolecular reaction rule (a reaction rule with 2 reactants)
        can have 0 (empty list) or 1 products.

        If there is only 1 reactant or 1 product, the brackets around it 
        can be omitted. 

        # Notes.
        Units for the rate of a bimolecular reaction rule:
            [rate] = meters^3 / (molecules * second)

        (Order of magnitude rate: 1e-18)

        For each possible combination of reactants (as described in the 
        syntax) for which no bimolecular reaction rule is specified, a
        repulsive bimolecular reaction rule is added by default (k=0).


        ### Surface binding reaction rules.
        # Syntax.
        add_reaction_rule([reactant], [product], rate)
        or
        add_reaction_rule(reactant, product, rate)

        where the reactant is of type:
            (ParticleType, Region), or
            ParticleType, for a ParticleType in the "world".

        and the product is of type:
            (ParticleType, Surface), or
            (0, Surface), for a Surface that absorps this reactant.

        # Notes.
        For each possible combination of reactant (as described in the 
        syntax) and Surface, for which no surface binding reaction rule 
        is specified (with the Surface as the second argument of the 
        product), a reflective surface binding reaction rule is added 
        by default (k=0).


        ### Surface unbinding reaction rules.
        # Syntax.
        add_reaction_rule([reactant], [product], rate)
        or
        add_reaction_rule(reactant, product, rate)

        where the reactant is of type:
            (ParticleType, Surface)

        and the product is of type:
            (ParticleType, Region), or
            ParticleType, for a ParticleType in the "world".

        # Notes.
        Unbinding from a Surface to a Region or to the "world" is a 
        Poissonian process.

        """
        if not isinstance(reactants, list):
            reactants = [reactants]
        if not isinstance(products, list):
            products = [products]

        reactants = map(self._get_species_type, reactants)
        products  = map(self._get_species_type, products)

        rr = ReactionRule(reactants, products)
        rr['k'] = '%.16g' % k
        self.network_rules.add_reaction_rule(rr)
        return rr

    def _get_species_type_name(self, particle_type, structure):
        """Helper.

        Construct new name for SpeciesType from ParticleType name and 
        Region or Surface id.

        """
        if structure.id == "world":
            return particle_type.name
        else:
            return '(' + particle_type.name + ', ' + str(structure.id) + ')'

    def _get_species_type(self, key):
        """Helper.

        (ParticleType, Region)        -> SpeciesType
        or
        (ParticleType, Surface)       -> SpeciesType
        or
        ParticleType (in the "world") -> SpeciesType

        """
        if isinstance(key, ParticleType):
            particle_type = key
            structure = self.get_structure("world")
        else:
            # Unpack (ParticleType, Region)-key.
            particle_type = key[0]
            structure = key[1]

            if particle_type == 0:
                # This is the virtual product of a decay or surface 
                # absorption reaction.
                particle_type = DummyParticleType(structure)

        name = self._get_species_type_name(particle_type, structure)

        return self._get_species_type_by_name(name)

    def _get_species_type_by_name(self, name):
        """Helper.

        """
        for species_type in self.species_types:
            if species_type['name'] == name:
                return species_type

        raise RuntimeError('SpeciesType %s does not exist.' % (name))

