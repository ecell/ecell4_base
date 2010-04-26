import bd
from weakref import ref

from gfrdbase import *
from utils import *
import itertools

import _gfrd

import os

BDPropagator = eval(os.environ.get("BDPropagator", "_gfrd._BDPropagator"))

class _MultiParticleContainer(_gfrd._ParticleContainer):
    def __init__(self, world):
        _gfrd._ParticleContainer.__init__(self)
        self.world = world
        self.particles = {}

    def num_particles(self):
        return len(self.particles)
    num_particles = property(num_particles)

    def get_surface(self, id):
        return self.world.get_surface(id)

    def get_species(self, id):
        return self.world.get_species(id)

    def species(self):
        return self.world.species
    species = property(species)

    def new_particle(self, species_id, position):
        retval = self.world.new_particle(species_id, position)
        self.particles[retval[0]] = retval[1]
        return retval

    def update_particle(self, pid_particle_pair):
        self.particles[pid_particle_pair[0]] = pid_particle_pair[1]
        return self.world.update_particle(pid_particle_pair)

    def remove_particle(self, pid):
        del self.particles[pid]
        self.world.remove_particle(pid)

    def get_particle(self, pid):
        p = self.particles.get(pid, None)
        if p is None:
            raise NotFound
        return pid, p

    def check_overlap(self, sphere, *arg):
        if len(arg) == 0:
            ignores = ()
        elif len(arg) == 1:
            if isinstance(arg[0], _gfrd.ParticleID):
                ignores = (arg[0],)
            else:
                ignores = arg[0]
        elif len(arg) == 2:
            assert all(isinstance(a, _gfrd.ParticleID) for a in arg)
            ignores = arg

        retval = []
        for pp in self.particles.iteritems():
            if pp[0] in ignores:
                continue
            dist = _gfrd.distance(pp[1].position, sphere[0]) - pp[1].radius
            if dist < sphere[1]:
                retval.append((pp, dist))
        retval.sort(lambda a, b: cmp(a[1], b[1]))
        return retval

    def distance(self, x, y):
        return self.world.distance(x, y)

    def apply_boundary(self, x):
        return self.world.apply_boundary(x)

    def cyclic_transpose(self, x, y):
        return self.world.cyclic_transpose(x, y)

    def __iter__(self):
        return self.particles.iteritems()

    def create_transaction(self):
        return _gfrd.TransactionImpl(self)


MultiParticleContainer = _gfrd._MultiParticleContainer

class Multi(object):
    def __init__(self, domain_id, main, dt_factor):
        self.main = ref(main)
        self.domain_id = domain_id
        self.event_id = None
        self.last_event = None
        self.sphere_container = _gfrd.SphericalShellContainer(main.world.world_size, 3)
        self.particle_container = MultiParticleContainer(main.world)
        self.escaped = False
        self.dt_factor = dt_factor
        self.last_reaction = None

    def initialize(self, t):
        self.last_time = t
        self.start_time = t
        main = self.main()
        self.dt = self.dt_factor * bd.calculate_bd_dt(main.world.get_species(sid) for sid in main.world.species)

    def get_multiplicity(self):
        return self.particle_container.num_particles
    multiplicity = property(get_multiplicity)

    def within_shell(self, pp):
        return bool(self.sphere_container.get_neighbors_within_radius(pp[1].position, -pp[1].radius))

    def add_shell(self, shell_id_shell_pair):
        self.sphere_container.update(shell_id_shell_pair)

    def add_particle(self, pid_particle_pair):
        if __debug__:
            log.info("add_particle: %s\n", pid_particle_pair)
        self.particle_container.update_particle(pid_particle_pair)

    def step(self):
        self.escaped = False
        tx = self.particle_container.create_transaction()
        main = self.main()
        ppg = BDPropagator(tx, main.network_rules,
                     myrandom.rng, self.dt, main.dissociation_retry_moves,
                     [pid for pid, _ in self.particle_container])

        self.last_event = None
        while ppg():
            if ppg.reactions:
                self.last_event = _gfrd.EventType.MULTI_REACTION
                self.last_reaction = ppg.reactions[-1]
                break

        for pid_particle_pair in itertools.chain(
                tx.modified_particles, tx.added_particles):
            overlapped = main.world.check_overlap(pid_particle_pair[1].shape, pid_particle_pair[0])
            if overlapped:
                if __debug__:
                    log.info("collision occurred between particles of a multi and the outside: %s - %s.  moves will be rolled back." % (pid_particle_pair, list(overlapped)))
                tx.rollback()
                self.step()
                return

            if not self.within_shell(pid_particle_pair):
                self.last_event = _gfrd.EventType.MULTI_ESCAPE
                main.clear_volume(
                    pid_particle_pair[1].position,
                    pid_particle_pair[1].radius, ignore=[self.domain_id, ])

    def check(self):
        # shells are contiguous
        # FIXME: this code cannot detect a pair of shells that are isolated
        #        from others.
        for _, shell in self.shell_list:
            result = self.sphere_container.get_neighbors(shell.shape.position)
            # Check contiguity with nearest neighbor only (get_neighbors 
            # returns a sorted list).
            nearest = result[1]
            distance = nearest[1]
            assert distance - shell.shape.radius < 0.0,\
                'shells of %s are not contiguous.' % str(self.multiref())

        # all particles within the shell.
        for pid_particle_pair in self.particle_container:
            assert self.within_shell(pid_particle_pair),\
                'not all particles within the shell.'

        main = self.main()
        for shell_id, shell in self.shell_list:
            try:
                container = main.get_container(shell)
                container[shell_id]
            except:
                raise RuntimeError,\
                    'self.sim.main.sphere_container does not contain %s'\
                    % str(shell_id)

    def __repr__(self):
        return 'Multi[domain_id=%s, event_id=%s: %s: %s]' % (
            self.domain_id, self.event_id,
            ', '.join('(%s, %s)' % (p[0], repr(p[1])) for p in self.particle_container),
            ', '.join('(%s, %s)' % (s[0], repr(s[1])) for s in self.sphere_container)
            )

    def has_particle(self, pid):
        try:
            self.particle_container.get_particle(pid)
            return True
        except:
            return False

    def particles(self):
        return iter(self.particle_container)
    particles = property(particles)

    def num_shells(self):
        return len(self.sphere_container)
    num_shells = property(num_shells)    

    def shell_list(self):
        return iter(self.sphere_container)
    shell_list = property(shell_list)
