from bd import BDSimulatorCoreBase, DEFAULT_DT_FACTOR
from weakref import ref

from _gfrd import *
from gfrdbase import *

from utils import *


class MultiBDCore(BDSimulatorCoreBase):
    '''
    Used internally by Multi.
    '''
    def __init__(self, main, multi):

        BDSimulatorCoreBase.__init__(self, main)

        # this has to be ref, not proxy, since it is used for comparison.
        self.multiref = ref(multi)

        self.particle_matrix = ParticleContainer(self.main.world.world_size, 3)
        self.sphere_container = SphericalShellContainer(self.main.world.world_size, 3)
        self.escaped = False

    def update_particle(self, pid_particle_pair):
        self.particle_matrix.update(pid_particle_pair)
        self.main.move_particle(pid_particle_pair)

    def initialize(self):
        BDSimulatorCoreBase.initialize(self)

    def step(self):
        self.escaped = False
        BDSimulatorCoreBase.step(self)

    def clear_volume(self, pos, radius, ignore=[]):
        if not self.within_shell(pos, radius):
            self.escaped = True
            self.clear_outer_volume(pos, radius, ignore)

    def clear_outer_volume(self, pos, radius, ignore=[]):
        self.main.clear_volume(pos, radius, ignore=[self.multiref().domain_id, ])
        if self.main.get_particles_within_radius(pos, radius, ignore):
            raise NoSpace()

    def within_shell(self, pos, radius):
        result = self.sphere_container.get_neighbors_within_radius(pos, - radius)
        return bool(result)
        
    def check_overlap(self, pos, radius, ignore=[]):
        result = self.particle_matrix.get_neighbors_within_radius(pos, radius)
        for item in result:
            if item[0][0] not in ignore:
                return item
        return None

    def get_particles_within_radius(self, pos, radius, ignore=[]):
        result = self.particle_matrix.get_neighbors_within_radius(pos, radius)
        return [n for n in result if n[0][0] not in ignore]

    def check(self):
        BDSimulatorCoreBase.check(self)

        # shells are contiguous
        for (_, shell) in self.multiref().shell_list:
            result = self.sphere_container.get_neighbors(shell.shape.position)
            # Check contiguity with nearest neighbor only (get_neighbors 
            # returns a sorted list).
            nearest = result[1]
            distance = nearest[1]
            assert distance - shell.shape.radius < 0.0,\
                'shells of %s are not contiguous.' % str(self.multiref())

        # all particles within the shell.
        for pid in self.particle_list:
            p = self.main.world.get_particle(pid)[1]
            assert self.within_shell(p.position, p.radius),\
                'not all particles within the shell.'


class Multi(object):
    def __init__(self, domain_id, main):
        self.domain_id = domain_id
        self.event_id = None
        self.sim = MultiBDCore(main, self)
        self.pid_shell_id_map = {}

    def initialize(self, t):
        self.last_time = t
        self.start_time = t

        self.sim.initialize() # ??

    def get_dt(self):
        return self.sim.dt

    dt = property(get_dt)

    def get_multiplicity(self):
        return len(self.sim.particle_list)

    multiplicity = property(get_multiplicity)

    def __add_particle(self, pid_particle_pair):
        self.sim.add_particle(pid_particle_pair)
        return pid_particle_pair

    def __add_shell(self, position, size):
        shell_id_shell_pair = (
            self.sim.main.shell_id_generator(),
            SphericalShell(self.domain_id, Sphere(position, size)))
        self.sim.main.move_shell(shell_id_shell_pair)
        self.sim.sphere_container.update(shell_id_shell_pair)
        return shell_id_shell_pair

    def add_particle_and_shell(self, pid_particle_pair, shell_size):
        self.__add_particle(pid_particle_pair)
        shell_id_shell_pair = self.__add_shell(pid_particle_pair[1].position, shell_size)
        self.pid_shell_id_map[pid_particle_pair[0]] = shell_id_shell_pair[0]
        return pid_particle_pair, shell_id_shell_pair

    def check(self):
        self.sim.check()

        for (shell_id, shell) in self.shell_list:
            try:
                container = self.sim.main.get_container(shell)
                container[shell_id]
            except:
                raise RuntimeError,\
                    'self.sim.main.sphere_container does not contain %s'\
                    % str(shell_id)

    def __repr__(self):
        return 'Multi[%s: %s: event_id=%s]' % (
            self.domain_id,
            ', '.join(repr(p) for p in self.sim.particle_list),
            self.event_id)

    def get_shell_list(self):
        return self.sim.sphere_container
    shell_list = property(get_shell_list)


