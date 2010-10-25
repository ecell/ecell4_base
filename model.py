import _gfrd
from _gfrd import create_cuboidal_region, create_cylindrical_surface, \
        create_planar_surface
import numpy

__all__ = [
    'Species',
    'ParticleModel',
    'create_unimolecular_reaction_rule',
    'create_decay_reaction_rule',
    'create_annihilation_reaction_rule',
    'create_binding_reaction_rule',
    'create_unbinding_reaction_rule',

    # From _gfrd. Should be part of the model class.
    'create_cuboidal_region',
    'create_cylindrical_surface',
    'create_planar_surface',
    ]


# Define _gfrd docstrigns here, much easier to format than in C++.
_gfrd.create_cuboidal_region.__doc__ = \
"""create_cuboidal_region(id, corner, diagonal)

Create and return a new cuboidal Region.

Arguments:
    - id
        a descriptive name.
    - corner
        the point [x, y, z] of the cuboidal Region closest to
        [0, 0, 0]. Units: [meters, meters, meters]
    - diagonal
        the vector [x, y, z] from the corner closest to [0, 0, 0], to 
        the corner furthest away from [0, 0, 0]. Units:
        [meters, meters, meters]

"""

_gfrd.create_cylindrical_surface.__doc__ = \
"""create_cylindrical_surface(id, corner, radius, orientation, length)

Create and return a new cylindrical Surface.

Arguments:
    - id
        a descriptive name.
    - corner
        the point [x, y, z] on the axis of the cylinder closest to 
        [0, 0, 0]. Units: [meters, meters, meters]
    - radius
        the radius of the cylinder. Units: meters.
    - orientation
        the unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along the 
        axis of the cylinder.
    - length
        the length of the cylinder. Should be equal to the world_size. 
        Units: meters.

Surfaces are not allowed to touch or overlap.

"""

_gfrd.create_planar_surface.__doc__ = \
"""create_planar_surface(id, corner, unit_x, unit_y, length_x, length_y)

Create and return a new planar Surface.

Arguments:
    - id
        a descriptive name.
    - corner
        the point [x, y, z] on the plane closest to [0, 0, 0]. Units: 
        [meters, meters, meters]
    - unit_x
        a unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along the 
        plane.
    - unit_y
        a unit vector [1, 0, 0], [0, 1, 0] or [0, 0, 1] along the plane 
        and perpendicular to unit_x.
    - length_x
        the length of the plane along the unit vector unit_x. Should be 
        equal to the world_size. Units: meters.
    - length_y
        the length of the plane along the unit vector unit_y. Should be 
        equal to the world_size. Units: meters.

Surfaces are not allowed to touch or overlap.

"""

_gfrd.Model.add_species_type.im_func.__doc__ = \
"""add_species_type(self, species)

Add a Species to the ParticleModel.

Arguments:
    - species
        a Species created with the function model.Species.

"""


def Species(name, D, radius=0, structure="world", drift=0):
    """Define a new Species (in/on a specific Region or Surface).

    Arguments:
        - name
            the name of this Species.
        - D
            the diffusion constant for this Species in/on this 
            Region or Surface. Units: meters^2/second.
        - radius
            the radius for this Species in/on this Region or Surface. 
            Units: meters.
        - structure
            the Region or Surface in/on which this Species can exist.  
            Optional. If you do not specify a Structure the Species is 
            added to the "world".
        - drift
            the drift term for this ParticleType on a 
            CylindricalSurface (1D drift). Units: meters/second. 
            Optional.

    If a certain Species should be able to exist in the "world" as 
    well as in/on one of the previously created Regions or Surfaces, 
    then two distinct Species should be created. One with and one 
    without an explicit Structure argument.

    """
    st = _gfrd.SpeciesType()
    st["name"] = str(name)
    st["D"] = str(D)
    st["v"] = str(drift)
    st["radius"] = str(radius)
    st["structure"] = structure
    return st


class ParticleModel(_gfrd.Model):
    """
    """
    def __init__(self, world_size):
        """Create a new ParticleModel.

        Arguments:
            - world_size
                the size of one side of the simulation "world". Units: 
                meters.

        The simulation "world" is always assumed to be a cube with 
        *periodic boundary conditions*, with 1 corner at [0, 0, 0] and 
        the corner furthest away from [0, 0, 0] being at
        [world_size, world_size, world_size].

        """
        _gfrd.Model.__init__(self)
        self.world_size = world_size
        self.structures = {}

        # Particles of a Species whose Surface is not specified will be 
        # added to the "world". Dimensions don't matter, except for 
        # visualization.
        x = numpy.repeat(world_size / 2, 3)
        region = _gfrd.CuboidalRegion('world', _gfrd.Box(x, x))
        self.add_structure(region)

    def add_structure(self, structure):
        """Add a Structure (Region or Surface) to the ParticleModel.

        Arguments:
            - structure
              a Region or Surface created with one of the functions
              model.create_<>_region or model.create_<>_surface.

        """
        assert isinstance(structure, _gfrd.Structure)
        self.structures[structure.id] = structure
        return structure

    def add_reaction_rule(self, reaction_rule):
        """Add a ReactionRule to the ParticleModel.

        Argument:
            - reaction rule
                a ReactionRule created by one of the functions
                model.create_<>_reaction_rule.

        """
        self.network_rules.add_reaction_rule(reaction_rule)

    def get_structure(self, id): 
        return self.structures[id]

    def set_all_repulsive(self):
        """Set all 'other' possible ReactionRules to be repulsive.

        By default an EGFRDSimulator will assume:
            - a repulsive bimolecular reaction rule (k=0) for each 
              possible combination of reactants for which no 
              bimolecular reaction rule is specified. 
          
        This method explicitly adds these ReactionRules to the 
        ParticleModel.

        """
        nr = self.network_rules
        # Maybe the user has defined a reaction rule for any 2 species since a 
        # previous call to this method, so remove *all* repulsive reaction 
        # rules first.
        for species1 in self.species_types:
            for species2 in self.species_types:
                gen = nr.query_reaction_rule(species1, species2)
                if gen is not None:
                    for reaction_rule in gen:
                        if float(reaction_rule['k']) == 0.0:
                            nr.remove_reaction_rule(reaction_rule)

        for species1 in self.species_types:
            for species2 in self.species_types:
                gen = nr.query_reaction_rule(species1, species2)
                if gen is None or len(set(gen)) == 0:
                    rr = _gfrd.ReactionRule([species1, species2], [])
                    rr['k'] = '0.0'
                    nr.add_reaction_rule(rr)

def create_unimolecular_reaction_rule(reactant, product, k):
    """Example: A -> B.

    Arguments:
        - reactant
            a Species.
        - product 
            a Species.
        - k
            reaction rate. Units: per second. (Rough order of magnitude: 
            1e-2 /s to 1e2 /s).

    The reactant and the product should be in/on the same 
    Region or Surface.

    There is no distinction between an intrinsic and an overall reaction 
    rate for a unimolecular ReactionRule.

    """
    rr = _gfrd.ReactionRule([reactant], [product])
    rr['k'] = '%.16g' % k
    return rr

def create_decay_reaction_rule(reactant, k):
    """Example: A -> 0.

    Arguments:
        - reactant
            a Species.
        - k
            reaction rate. Units: per second. (Rough order of magnitude: 
            1e-2 /s to 1e2 /s).

    There is no distinction between an intrinsic and an overall reaction 
    rate for a decay ReactionRule.

    """
    rr = _gfrd.ReactionRule([reactant], [])
    rr['k'] = '%.16g' % k
    return rr

def create_annihilation_reaction_rule(reactant1, reactant2, ka):
    """Example: A + B -> 0.

    Arguments:
        - reactant1
            a Species.
        - reactant2
            a Species.
        - ka
            intrinsic reaction rate. Units: meters^3 per second. (Rough 
            order of magnitude: 1e-16 m^3/s to 1e-20 m^3/s).

    The reactants should be in/on the same Region or Surface.

    ka should be an *intrinsic* reaction rate. You can convert an 
    overall reaction rate (kon) to an intrinsic reaction rate (ka) with 
    the function utils.k_a(kon, kD).

    By default an EGFRDSimulator will assume a repulsive 
    bimolecular reaction rule (ka=0) for each possible combination of 
    reactants for which no bimolecular reaction rule is specified. 
    You can explicitly add these reactions to the model with the method 
    model.ParticleModel.set_all_repulsive.

    """
    rr = _gfrd.ReactionRule([reactant1, reactant2], [])
    rr['k'] = '%.16g' % ka
    return rr

def create_binding_reaction_rule(reactant1, reactant2, product, ka):
    """Example: A + B -> C.

    Arguments:
        - reactant1
            a Species.
        - reactant2
            a Species.
        - product
            a Species.
        - ka
            intrinsic reaction rate. Units: meters^3 per second. (Rough 
            order of magnitude: 1e-16 m^3/s to 1e-20 m^3/s)

    The reactants and the product should be in/on the same 
    Region or Surface.

    A binding reaction rule can not have more than 1 product.

    ka should be an *intrinsic* reaction rate. You can convert an 
    overall reaction rate (kon) to an intrinsic reaction rate (ka) with 
    the function utils.k_a(kon, kD).

    By default an EGFRDSimulator will assume a repulsive 
    bimolecular reaction rule (ka=0) for each possible combination of 
    reactants for which no bimolecular reaction rule is specified. 
    You can explicitly add these reactions to the model with the method
    model.ParticleModel.set_all_repulsive.

    """
    rr = _gfrd.ReactionRule([reactant1, reactant2], [product])
    rr['k'] = '%.16g' % ka
    return rr

def create_unbinding_reaction_rule(reactant, product1, product2, kd):
    """Example: A -> B + C.

    Arguments:
        - reactant
            a Species.
        - product1
            a Species.
        - product2
            a Species.
        - intrinsic_rate

    The reactant and the products should be in/on the same 
    Region or Surface.

    An unbinding reaction rule can not have more than 2 products.

    kd should be an *intrinsic* reaction rate. You can convert an 
    overall reaction rate (koff) for this reaction rule to an intrinsic 
    reaction rate (kd) with the function utils.k_d(koff, kon, kD) or 
    utils.k_d_using_ka(koff, ka, kD).

    """
    rr = _gfrd.ReactionRule([reactant], [product1, product2])
    rr['k'] = '%.16g' % kd
    return rr

