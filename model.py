import _gfrd
import numpy

__all__ = [
    'ParticleModel',
    'create_unimolecular_reaction_rule',
    'create_decay_reaction_rule',
    'create_binding_reaction_rule',
    'create_unbinding_reaction_rule',
    ]

def Species(name, D, radius, surface="world"):
    st = _gfrd.SpeciesType()
    st["name"] = str(name)
    st["D"] = str(D)
    st["radius"] = str(radius)
    st["surface"] = surface
    return st


class ParticleModel(_gfrd.Model):
    def __init__(self, world_size):
        _gfrd.Model.__init__(self)
        self.regions = {}

        # Particles of a Species whose surface is not specified will be added 
        # to the world. Dimensions don't matter, except for visualization.
        x = numpy.repeat(world_size / 2, 3)
        region = _gfrd.CuboidalRegion('world', _gfrd.Box(x, x))
        self.add_structure(region)

    def add_planar_surface(self, id, origin, unit_x, unit_y, Lx, Ly):
        """Add a planar surface.

        id -- a descriptive name, should not be omitted.

        origin -- [x0, y0, z0] is the *center* of the planar surface.
        vector_x -- [x1, y1, z1] and
        vector_y -- [x2, y2, z2] are 2 perpendicular vectors that don't have 
        to be normalized that span the plane. For example [1,0,0] and [0,1,0]
        for a plane at z=0.

        Lx -- lx and 
        Ly -- ly are the distances from the origin of the plane along vector_x 
            or vector_y *to an edge* of the plane. PlanarSurfaces are finite.
        Lz -- dz, the thickness of the planar surface, can be omitted for Lz=0.

        """
        unit_x = _gfrd.normalize(unit_x, 1.)
        unit_y = _gfrd.normalize(unit_y, 1.)
        return self.add_structure(
            _gfrd.PlanarSurface(id,
                _gfrd.Plane(origin, unit_x, unit_y, Lx, Ly)))

    def add_cylindrical_surface(self, id, origin, radius, orientation, size):
        """Add a cylindrical surface.

        id -- a descriptive name, should not be omitted.

        origin -- [x0, y0, z0] is the *center* of the cylindrical surface.
        radius -- r is the radis of the cylinder.
        orientation -- [x1, y1, z1] is a vector that doesn't have to
            normalized that defines the orienation of the cylinder. For 
            example [0,0,1] for a for a cylinder along the z-axis.
        size -- lz is the distances from the origin of the cylinder along 
            the oriention vector to the end of the cylinder. So effectively
            the *half-length*. CylindricalSurfaces are finite.

        """
        return self.add_structure(
            _gfrd.CylindricalSurface(id,
                _gfrd.Cylinder(origin, radius, orientation, size)))

    def add_structure(self, surface):
        assert isinstance(surface, _gfrd.Structure)
        self.regions[surface.id] = surface
        return surface

    def get_structure(self, id): 
        return self.regions[id]

    def set_all_repulsive(self):
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

    def add_reaction(self, reactants, products, k): 
        rr = _gfrd.ReactionRule(reactants, products)
        rr['k'] = '%.16g' % k
        self.network_rules.add_reaction_rule(rr)
        return rr

    def dump_reaction_rule(self, reaction_rule):
        '''Pretty print reaction rule.

        ReactionRule.__str__ would be good, but we are actually getting a 
        ReactionRuleInfo or ReactionRuleCache object.

        '''
        buf = ('k=%.3g' % reaction_rule.k + ': ').ljust(15)
        for index, sid in enumerate(reaction_rule.rt.reactants):
            if index != 0:
                buf += ' + '
            reactant = self.get_species_type_by_id(sid)
            buf += reactant['id'].ljust(15)
        if len(reaction_rule.products) == 0:
            if reaction_rule.k != 0:
                buf += '..decays'
        else:
            buf += '-> '

        for index, sid in enumerate(reaction_rule.products):
            if index != 0:
                buf += ' + '
            product = self.get_species_type_by_id(sid)
            buf += product['id'].ljust(15)

        return buf + '\n'


def create_unimolecular_reaction_rule(s1, p1, k):
    rr = _gfrd.ReactionRule([s1, ], [p1, ])
    rr['k'] = '%.16g' % k
    return rr


def create_decay_reaction_rule(s1, k):
    rr = _gfrd.ReactionRule([s1, ], [])
    rr['k'] = '%.16g' % k
    return rr


def create_binding_reaction_rule(s1, s2, p1, k):
    rr = _gfrd.ReactionRule([s1, s2], [p1, ])
    rr['k'] = '%.16g' % k
    return rr


def create_unbinding_reaction_rule(s1, p1, p2, k):
    rr = _gfrd.ReactionRule([s1, ], [p1, p2])
    rr['k'] = '%.16g' % k
    return rr



