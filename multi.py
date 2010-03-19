from bd import BDSimulatorCoreBase, DEFAULT_DT_FACTOR
from weakref import ref

from _gfrd import *
from gfrdbase import *

from utils import *


class MultiBDCore( BDSimulatorCoreBase ):
    '''
    Used internally by Multi.
    '''
    def __init__( self, main, multi ):

        BDSimulatorCoreBase.__init__( self, main )

        # this has to be ref, not proxy, since it is used for comparison.
        self.multiref = ref( multi )

        self.particleMatrix = ParticleContainer(self.main.world.world_size, self.main.world.matrix_size)
        self.sphere_container = SphericalShellContainer(self.main.world.world_size, self.main.world.matrix_size)
        self.escaped = False

    def updateParticle( self, pid_particle_pair ):
        self.particleMatrix.update( pid_particle_pair )
        self.main.moveParticle( pid_particle_pair )

    def initialize( self ):
        BDSimulatorCoreBase.initialize( self )

    def step( self ):
        self.escaped = False
        BDSimulatorCoreBase.step( self )

    def addParticle(self, pid_particle_pair):
        self.addToParticleList(pid_particle_pair[0])
        self.particleMatrix.update(pid_particle_pair)

    def removeParticle( self, pid_particle_pair):
        self.main.removeParticle(pid_particle_pair)
        self.removeFromParticleList(pid_particle_pair[0])
        del self.particleMatrix[pid_particle_pair[0]]

    def createParticle( self, sid, pos ):
        particle = self.main.createParticle( sid, pos )
        self.addParticle( particle )
        return particle

    def moveParticle( self, pid_particle_pair ):
        self.updateParticle( pid_particle_pair )

    def clearVolume( self, pos, radius, ignore=[] ):
        if not self.withinShell( pos, radius ):
            self.escaped = True
            self.clearOuterVolume( pos, radius, ignore )

    def clearOuterVolume( self, pos, radius, ignore=[] ):
        self.main.clearVolume( pos, radius, ignore=[self.multiref().domain_id,] )
        if self.main.getParticlesWithinRadius(pos, radius, ignore):
            raise NoSpace()

    def withinShell( self, pos, radius ):
        result = self.sphere_container.get_neighbors_within_radius( pos, - radius )
        return bool(result)
        
    def checkOverlap( self, pos, radius, ignore=[] ):
        result = self.particleMatrix.get_neighbors_within_radius( pos, radius )
        for item in result:
            if item[0][0] not in ignore:
                return item
        return None

    def getParticlesWithinRadius( self, pos, radius, ignore=[] ):
        result = self.particleMatrix.get_neighbors_within_radius( pos, radius )
        return [ n for n in result if n[0][0] not in ignore ]

    def check( self ):
        BDSimulatorCoreBase.check( self )

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
        for pid in self.particleList:
            p = self.main.world.get_particle(pid)[1]
            assert self.withinShell( p.position, p.radius ),\
                'not all particles within the shell.'


class Multi( object ):
    def __init__( self, domain_id, main ):
        self.domain_id = domain_id
        self.eventID = None
        self.sim = MultiBDCore( main, self )
        self.pid_shell_id_map = {}

    def initialize( self, t ):
        self.lastTime = t
        self.startTime = t

        self.sim.initialize() # ??

    def getDt( self ):
        return self.sim.dt

    dt = property( getDt )

    def getMultiplicity( self ):
        return len( self.sim.particleList )

    multiplicity = property( getMultiplicity )

    def __addParticle(self, pid_particle_pair):
        self.sim.addParticle(pid_particle_pair)
        return pid_particle_pair

    def __addShell(self, position, size):
        shell_id_shell_pair = (
            self.sim.main.shellIDGenerator(),
            SphericalShell(self.domain_id, Sphere(position, size)) )
        self.sim.main.moveShell(shell_id_shell_pair)
        self.sim.sphere_container.update(shell_id_shell_pair)
        return shell_id_shell_pair

    def addParticleAndShell(self, pid_particle_pair, shellSize):
        self.__addParticle(pid_particle_pair)
        shell_id_shell_pair = self.__addShell(pid_particle_pair[1].position, shellSize)
        self.pid_shell_id_map[pid_particle_pair[0]] = shell_id_shell_pair[0]
        return pid_particle_pair, shell_id_shell_pair

    def check( self ):
        self.sim.check()

        for (shell_id, shell) in self.shell_list:
            try:
                container = self.sim.main.get_container(shell)
                container[shell_id]
            except:
                raise RuntimeError,\
                    'self.sim.main.sphere_container does not contain %s'\
                    % str(shell_id)

    def __repr__( self ):
        return 'Multi[%s: %s: eventID=%s]' % (
            self.domain_id,
            ', '.join( repr( p ) for p in self.sim.particleList ),
            self.eventID )

    def getShellList(self):
        return self.sim.sphere_container
    shell_list = property(getShellList)


