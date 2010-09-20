import _gfrd
from _gfrd import create_planar_surface, create_cylindrical_surface, \
        create_cuboidal_region
import numpy

__all__ = [
    'Species',
    'ParticleModel',
    'create_unimolecular_reaction_rule',
    'create_decay_reaction_rule',
    'create_annihilation_reaction_rule',
    'create_binding_reaction_rule',
    'create_unbinding_reaction_rule',
    'create_planar_surface',
    'create_cylindrical_surface',
    'create_cuboidal_region'
    ]



def Species(name, D, radius=0, structure="world", drift=0):
    st = _gfrd.SpeciesType()
    st["name"] = str(name)
    st["D"] = str(D)
    st["v"] = str(drift)
    st["radius"] = str(radius)
    st["surface"] = structure
    return st


class ParticleModel(_gfrd.Model):
    def __init__(self, world_size):
        _gfrd.Model.__init__(self)
        self.world_size = world_size
        self.structures = {}

        # Particles of a Species whose surface is not specified will be 
        # added to the "world". Dimensions don't matter, except for 
        # visualization.
        x = numpy.repeat(world_size / 2, 3)
        region = _gfrd.CuboidalRegion('world', _gfrd.Box(x, x))
        self.add_structure(region)

    def add_structure(self, structure):
        assert isinstance(structure, _gfrd.Structure)
        self.structures[structure.id] = structure
        return structure

    def get_structure(self, id): 
        return self.structures[id]

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

    def dump_reaction_rule(self, reaction_rule):
        '''Pretty print reaction rule.

        ReactionRule.__str__ would be good, but we are actually getting a 
        ReactionRuleInfo or ReactionRuleCache object.

        '''
        buf = ('k=%.3g' % reaction_rule.k + ': ').ljust(15)
        for index, sid in enumerate(reaction_rule.reactants):
            if index != 0:
                buf += ' + '
            reactant = self.get_species_type_by_id(sid)
            buf += reactant['name'].ljust(15)
        if len(reaction_rule.products) == 0:
            if reaction_rule.k != 0:
                buf += '..decays'
            else:
                buf += '..reflective'
        else:
            buf += '-> '

        for index, sid in enumerate(reaction_rule.products):
            if index != 0:
                buf += ' + '
            product = self.get_species_type_by_id(sid)
            buf += product['name'].ljust(15)

        return buf + '\n'


def create_unimolecular_reaction_rule(s1, p1, k):
    rr = _gfrd.ReactionRule([s1, ], [p1, ])
    rr['k'] = '%.16g' % k
    return rr

def create_decay_reaction_rule(s1, k):
    rr = _gfrd.ReactionRule([s1, ], [])
    rr['k'] = '%.16g' % k
    return rr

def create_annihilation_reaction_rule(s1, s2, k):
    rr = _gfrd.ReactionRule([s1, s2], [])
    rr['k'] = '%.16g' % k
    return rr

def create_binding_reaction_rule(s1, s2, p1, k):
    rr = _gfrd.ReactionRule([s1, s2], [p1])
    rr['k'] = '%.16g' % k
    return rr

def create_unbinding_reaction_rule(s1, p1, p2, k):
    rr = _gfrd.ReactionRule([s1, ], [p1, p2])
    rr['k'] = '%.16g' % k
    return rr



