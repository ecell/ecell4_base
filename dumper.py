from egfrd import Single, Pair, Multi
__all__ = [
    'dump_species',
    'dump_particles',
    'dump_particles_by_sid',
    'dump_domains',
    ]

def dump_species(sim):
    return sim.world.species

def dump_particles_by_sid(sim, sid):
    for pid in sim.world.get_particle_ids(sid):
        particle = sim.world.get_particle(pid)[1]
        yield (pid, particle)

def dump_particles(sim):
    return sim.world

def dump_domains(egfrdsim):
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
