import bd
from weakref import ref

from gfrdbase import *
from utils import *
import itertools

import _gfrd
from constants import EventType

import os

class Multi(object):
    def __init__(self, domain_id, main, dt_factor):
        self.main = ref(main)
        self.domain_id = domain_id
        self.event_id = None
        self.last_event = None
        self.sphere_container = _gfrd.SphericalShellContainer(main.world.world_size, 3)
        self.particle_container = _gfrd.MultiParticleContainer(main.world)
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
        if __debug__:
            log.info("add shell to multi:\n  (%s,\n   %s)" %
                     (shell_id_shell_pair[0], shell_id_shell_pair[1]))
        self.sphere_container.update(shell_id_shell_pair)

    def add_particle(self, pid_particle_pair):
        if __debug__:
            log.info("add particle to multi:\n  (%s,\n   %s)" % 
                     (pid_particle_pair[0], pid_particle_pair[1]))
        self.particle_container.update_particle(pid_particle_pair)

    def step(self):
        self.escaped = False
        tx = self.particle_container.create_transaction()
        main = self.main()

        class check_reaction(object):
            def __init__(self):
                self.reactions = []

            def __call__(self, ri):
                self.reactions.append(ri)

        class clear_volume(object):
            def __init__(self, outer):
                self.outer_ = outer

            def __call__(self, shape, ignore0, ignore1=None):
                within_radius = bool(
                    self.outer_.sphere_container.get_neighbors_within_radius(
                        shape.position, -shape.radius))
                if not within_radius:
                    main = self.outer_.main()
                    if self.outer_.last_event == None:
                        self.outer_.last_event = EventType.MULTI_ESCAPE
                    main.clear_volume(shape.position, shape.radius, 
                                      ignore=[self.outer_.domain_id, ])
                    if ignore1 is None:
                        return not main.world.check_overlap(
                            (shape.position, shape.radius), ignore0)
                    else:
                        return not main.world.check_overlap(
                            (shape.position, shape.radius), ignore0, ignore1)
                return True

        cr = check_reaction()
        vc = clear_volume(self)

        ppg = _gfrd.BDPropagator(tx, main.network_rules,
                     myrandom.rng, self.dt, main.dissociation_retry_moves,
                     cr, vc, [pid for pid, _ in self.particle_container])

        self.last_event = None
        while ppg():
            if cr.reactions:
                self.last_reaction = cr.reactions[-1]
                if len(self.last_reaction.reactants) == 1:
                    self.last_event = EventType.MULTI_UNIMOLECULAR_REACTION
                else:
                    self.last_event = EventType.MULTI_BIMOLECULAR_REACTION
                break

        # for pid_particle_pair in itertools.chain(
        #         tx.modified_particles, tx.added_particles):
        #     overlapped = main.world.check_overlap(pid_particle_pair[1].shape, pid_particle_pair[0])
        #     if overlapped:
        #         if __debug__:
        #             log.info("collision occurred between particles of a multi and the outside: %s - %s.  moves will be rolled back." % (pid_particle_pair, list(overlapped)))
        #         tx.rollback()
        #         return

        #     if not self.within_shell(pid_particle_pair):
        #         if self.last_event == None:
        #             self.last_event = EventType.MULTI_ESCAPE
        #         main.clear_volume(
        #             pid_particle_pair[1].position,
        #             pid_particle_pair[1].radius, ignore=[self.domain_id, ])

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
                'shells of %s are not contiguous.' % str(self)

        # all particles within the shell.
        for pid_particle_pair in self.particle_container:
            assert self.within_shell(pid_particle_pair),\
                'not all particles within the shell.'

        main = self.main()
        for shell_id, shell in self.shell_list:
            container = main.get_container(shell)
            if not container.contains(shell_id):
                raise RuntimeError,\
                    'self.sim.main.sphere_container does not contain %s'\
                    % str(shell_id)
        for shell_id, shell in main.containers[0]:
            if shell.did == self.domain_id:
                if not self.sphere_container.contains(shell_id):
                    raise RuntimeError,\
                        'self.sphere_container does not contain %s'\
                        % str(shell_id)

    def __repr__(self):
        return 'Multi[domain_id=%s, event_id=%s,\n    %s,\n    %s]' % (
            self.domain_id, self.event_id,
            ',\n    '.join('(%s,\n     %s)' % (p[0], repr(p[1])) for p in self.particle_container),
            ',\n    '.join('(%s,\n     %s)' % (s[0], repr(s[1])) for s in self.sphere_container)
            )

    def has_particle(self, pid):
        return pid in [pid_particle_pair[0] for pid_particle_pair in self.particles]
        # try:
        #     self.particle_container.get_particle(pid)
        #     return True
        # except:
        #     return False

    def particles(self):
        return iter(self.particle_container)
    particles = property(particles)

    def num_shells(self):
        return len(self.sphere_container)
    num_shells = property(num_shells)    

    def shell_list(self):
        return iter(self.sphere_container)
    shell_list = property(shell_list)
