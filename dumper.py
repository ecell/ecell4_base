from single import Single
from pair import Pair
from multi import Multi
from egfrd import EGFRDSimulator
from gillespie import GillespieSimulator
import _gfrd

# get methods return an iterator, dump methods return a string.
__all__ = [
    'get_species',
    'dump_species',
    'get_species_names',
    'dump_species_names',
    'get_particles',
    'dump_particles',
    'get_number_of_particles',
    'dump_number_of_particles',
    'get_domains',
    'dump_domains',
    'get_reaction_rules',
    'dump_reaction_rules',
    ]


def get_species(sim):
    """Return an iterator over the Species in the simulator.

    Arguments:
        - sim
            an EGFRDSimulator.

    """
    # sim.world.species returns a SpeciesRange over SpeciesInfo.
    return sim.world.model.species_types

def dump_species(sim):
    """Return a string containing the Species in the simulator.

    Arguments:
        - sim
            an EGFRDSimulator.

    """
    return '\n'.join((str(st) for st in get_species(sim)))

def get_species_names(sim):
    """Return an iterator over the names of the Species in the 
    simulator.

    Arguments:
        - sim
            an EGFRDSimulator.

    """
    return (st['name'] for st in get_species(sim))

def dump_species_names(sim):
    """Return a string containing the names of the Species in the 
    simulator.

    Arguments:
        - sim
            an EGFRDSimulator.

    """
    return ' '.join(get_species_names(sim))

def _get_species_type_by_name(sim, name):
    #Helper.
    for species_type in sim.world.model.species_types:
        if species_type['name'] == name:
            return species_type

    raise RuntimeError('SpeciesType %s does not exist.' % (name))

def _get_particles_by_sid(sim, sid):
    # Helper.
    for pid in sim.world.get_particle_ids(sid):
        particle = sim.world.get_particle(pid)[1]
        yield (pid, particle)

def get_particles(sim, identifier=None):
    """Return an iterator over the
    (particle identifier, particle)-pairs in the simulator.

    Arguments:
        - sim
            an EGFRDSimulator.
        - identifier
            a Species or the name of a Species. If none is specified, 
            all (particle identifier, particle)-pairs will be returned.

    """
    if isinstance(sim, GillespieSimulator):
        raise RuntimeError('GillespieSimulator does not keep track '
                           'of individual particles.')

    if identifier == None:
        return sim.world
    else:
        if isinstance(identifier, _gfrd.SpeciesType):
            sid = identifier
        elif isinstance(identifier, str):
            sid = _get_species_type_by_name(sim, identifier).id
        else:
            raise RuntimeError('Wrong identifier type.')
        return _get_particles_by_sid(sim, sid)

def dump_particles(sim, identifier=None):
    """Return a string containing the
    (particle identifier, particle)-pairs in the simulator.

    Arguments:
        - sim
            an EGFRDSimulator.
        - identifier
            a Species or the name of a Species. If none is specified, 
            all (particle identifier, particle)-pairs will be returned.

    """
    return '\n'.join((str(x) for x in get_particles(sim, identifier)))

def _get_number_of_particles_by_sid(sim, sid):
    # Helper.
    if isinstance(sim, EGFRDSimulator):
        return len(sim.world.get_particle_ids(sid))
    else:
        # Gillespie.
        species_index = sim.speciesDict[sid]
        return sim.stateArray[species_index]

def get_number_of_particles(sim, identifier=None):
    """Return the number of particles of a certain Species in the 
    simulator.

    Arguments:
        - sim
            either an EGFRDSimulator or a GillespieSimulator.
        - identifier
            a Species. Optional. If none is specified, a list of 
            (Species name, number of particles)-pairs will be returned.

    """
    if identifier == None:
        if isinstance(sim, EGFRDSimulator):
            return [(st["name"], _get_number_of_particles_by_sid(sim, st.id))
                    for st in sim.world.model.species_types]
        else:
            return sim.stateArray
    else:
        if isinstance(identifier, _gfrd.SpeciesType):
            sid = identifier.id
        #elif isinstance(identifier, str):
        #    sid = _get_species_type_by_name(sim, identifier).id
        else:
            raise RuntimeError('Wrong identifier type.')
        return _get_number_of_particles_by_sid(sim, sid)

def dump_number_of_particles(sim, identifier=None):
    """Return a string containing the number of particles of a certain 
    Species in the simulator.

    Arguments:
        - sim
            either an EGFRDSimulator or a GillespieSimulator.
        - identifier
            a Species. Optional. If none is specified, 
            a string of (Species name, number of particles)-pairs will 
            be returned.

    """
    return str(get_number_of_particles(sim, identifier))

def get_domains(egfrdsim):
    """Return an iterator over the protective domains in the simulator.

    """
    for did, domain in egfrdsim.domains.iteritems():
        shell_list = domain.shell_list
        pid_particle_pair_list = []

        if isinstance(domain, Single):
            pid_particle_pair_list = [domain.pid_particle_pair]
        elif isinstance(domain, Pair):
            pid_particle_pair_list = [
                domain.single1.pid_particle_pair,
                domain.single2.pid_particle_pair]
        elif isinstance(domain, Multi):
            pid_particle_pair_list = []
            for pid_particle_pair in domain.particle_container:
                pid_particle_pair_list.append(pid_particle_pair)

        yield ((did, domain), pid_particle_pair_list, shell_list) 

def dump_domains(egfrdsim):
    """Return an string containing the protective domains in the 
    simulator.

    """
    return '\n'.join((str(domain) for domain in get_domains(egfrdsim)))

def get_reaction_rules(model_or_simulator):
    """Return three lists with all the reaction rules defined in the 
    ParticleModel or EGFRDSimulator.

    The three lists are:
        - reaction rules of only one reactant.
        - reaction rules between two reactants with a reaction rate 
          larger than 0.
        - repulsive reaction rules between two reactants with a 
          reaction rate equal to 0.

    Arguments:
        - model_or_simulator
            a ParticleModel or EGFRDSimulator.

    """
    if isinstance(model_or_simulator, EGFRDSimulator):
        model = model_or_simulator.world.model
    else:
        model = model_or_simulator

    # Return 3 lists with different types of reaction rules. 
    reaction_rules_1 = []
    reaction_rules_2 = []
    repulsive_rules = []

    # Wrap the network_rules first, the iterator over the products 
    # of the unwrapped one fails when there are no products.
    network_rules = _gfrd.NetworkRulesWrapper(model.network_rules)

    for index_of_si1, si1 in enumerate(model.species_types):
        rri_vector = network_rules.query_reaction_rule(si1)
        for rr_info in rri_vector:
            reaction_rules_1.append(rr_info)
        for si2 in list(model.species_types)[index_of_si1:]:
            rri_vector = network_rules.query_reaction_rule(si1, si2)
            for rr_info in rri_vector:
                if rr_info.k > 0:
                    reaction_rules_2.append(rr_info)
                else:
                    repulsive_rules.append(rr_info)

    return reaction_rules_1, reaction_rules_2, repulsive_rules

def _dump_reaction_rule(model, reaction_rule):
    # Helper. Return ReactionRule as string.

    #ReactionRule.__str__ would be good, but we are actually getting a 
    #ReactionRuleInfo or ReactionRuleCache object.
    buf = ('k=%.3g' % reaction_rule.k + ': ').ljust(15)
    for index, sid in enumerate(reaction_rule.reactants):
        if index != 0:
            buf += ' + '
        reactant = model.get_species_type_by_id(sid)
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
        product = model.get_species_type_by_id(sid)
        buf += product['name'].ljust(15)

    return buf + '\n'

def dump_reaction_rules(model_or_simulator):
    """Return a formatted string containing all the reaction rules 
    defined in the ParticleModel or EGFRDSimulator.

    Arguments:
        - model_or_simulator
            a ParticleModel or EGFRDSimulator.

    """
    if isinstance(model_or_simulator, EGFRDSimulator):
        model = model_or_simulator.world.model
    else:
        model = model_or_simulator

    rr1, rr2, rrr = get_reaction_rules(model)

    reaction_rules_1 = [_dump_reaction_rule(model, rule) for rule in rr1]
    reaction_rules_2 = [_dump_reaction_rule(model, rule) for rule in rr2]
    repulsive_rules  = [_dump_reaction_rule(model, rule) for rule in rrr]

    if repulsive_rules == []:
        repulsive_rules_as_string = (
            'None.\n  '
            'An EGFRDSimulator assumes all other possible\n'
            'reaction rules to be repulsive. You can explicitly add\n'
            'these repulsive reaction rules to the model with the\n'
            'method ParticleModel.set_all_repulsive.')
    else:
        repulsive_rules_as_string = ''.join(repulsive_rules)

    return('\nMonomolecular reaction rules:\n' +
           ''.join(reaction_rules_1) +
           '\nBimolecular reaction rules:\n' +
           ''.join(reaction_rules_2) +
           '\nRepulsive bimolecular reaction rules:\n' +
           repulsive_rules_as_string
           )

