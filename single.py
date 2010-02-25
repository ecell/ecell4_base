from _gfrd import *
from greens_function_wrapper import *

class Single( object ):
    """There are 2 main types of Singles:
        * NonInteractionSingle
        * InteractionSingle (when the particle is nearby a surface)

    Each type of Single defines a list of coordinates, see coordinate.py. For each 
    coordinate the Green's function is specified.

    """
    def __init__( self, domain_id, pid_particle_pair, shell_id_shell_pair, reactiontypes ):
        self.multiplicity = 1

        self.pid_particle_pair = pid_particle_pair
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.shell_list = [shell_id_shell_pair, ]

        self.eventID = None

        self.domain_id = domain_id

        self.updatek_tot()

    def getD( self ):
        return self.pid_particle_pair[1].D
    D = property( getD )

    def getMinRadius(self):
        return self.pid_particle_pair[1].radius
    minRadius = property(getMinRadius)

    def getShell(self):
        return self.shell_list[0]
    def setShell(self, value):
        self.shell_list[0] = value
    shell = property(getShell, setShell)

    def initialize( self, t ):
        '''
        Initialize this Single.

        The radius (shell size) is shrunken to the radius of the particle
        it represents.   
        self.lastTime is reset to the current time, and self.dt
        is set to zero.
        '''
        self.reset()
        self.lastTime = t

    def reset( self ):
        '''
        Reset the Single.

        Radius (shell size) is shrunken to the actual radius of the particle.
        self.dt is reset to 0.0.  Do not forget to reschedule this Single
        after calling this method.
        '''
        self.dt = 0.0
        self.eventType = EventType.SINGLE_ESCAPE

    def isReset( self ):
        return self.dt == 0.0 and self.eventType == EventType.SINGLE_ESCAPE

    def draw_reaction_time_tuple( self ):
        """Return a (reaction time, event type)-tuple.

        """
        if self.k_tot == 0:
            dt = numpy.inf
        elif self.k_tot == numpy.inf:
            dt = 0.0
        else:
            dt = (1.0 / self.k_tot) * math.log(1.0 / myrandom.uniform())
        return dt, EventType.SINGLE_REACTION

    def draw_interaction_time(self):
        """Todo.
        
        Note: we are not calling single.drawEventType() just yet, but 
        postpone it to the very last minute (when this event is executed in 
        fireSingle). So IV_EVENT can still be an iv escape or an iv 
        interaction.

        """
        pass

    def draw_escape_or_interaction_time_tuple(self):
        """Return an (escape or interaction time, event type)-tuple.

        Handles also all interaction events.
        
        """
        if self.getD() == 0:
            dt = numpy.inf
        else:
            dt = draw_time_wrapper(self.greens_function())

        event_type = EventType.SINGLE_ESCAPE
        return dt, event_type

    def determineNextEvent(self):
        """Return an (event time, event type)-tuple.

        """
        return min(self.draw_escape_or_interaction_time_tuple(),
                   self.draw_reaction_time_tuple())

    def updatek_tot( self ):
        self.k_tot = 0

        if not self.reactiontypes:
            return

        for rt in self.reactiontypes:
            self.k_tot += rt.k

    def drawReactionRule( self ):
        k_array = [ rt.k for rt in self.reactiontypes ]
        k_array = numpy.add.accumulate( k_array )
        k_max = k_array[-1]

        rnd = myrandom.uniform()
        i = numpy.searchsorted( k_array, rnd * k_max )

        return self.reactiontypes[i]

    def check( self ):
        pass

    def __repr__( self ):
        return 'Single[%s: %s: eventID=%s]' % ( self.domain_id, self.pid_particle_pair[0], self.eventID )


class NonInteractionSingle(Single):
    """1 Particle inside a shell, no other particles around. 

    There are 3 types of NonInteractionSingles:
        * SphericalSingle: spherical shell, 3D movement.
        * PlanarSurfaceSingle: cylindrical shell, 2D movement.
        * CylindricalSurfaceSingle: cylindrical shell, 1D movement.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id_shell_pair, 
                 reactiontypes):
        Single.__init__(self, domain_id, pid_particle_pair, shell_id_shell_pair,
                        reactiontypes)

    def getMobilityRadius(self):
        return self.shell[1].radius - self.pid_particle_pair[1].radius

    def drawNewPosition(self, dt, eventType):
        gf = self.greens_function()
        a = self.getMobilityRadius()
        r = draw_displacement_wrapper(gf, dt, eventType, a)
        displacement = self.displacement(r)
        assert abs(length(displacement) - abs(r)) <= 1e-15 * abs(r)
        return self.pid_particle_pair[1].position + displacement


class SphericalSingle(NonInteractionSingle):
    """1 Particle inside a (spherical) shell not on any surface.

        * Particle coordinate inside shell: r, theta, phi.
        * Coordinate: radial r.
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector: theta, phi.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes):
        # Create shell.
        position = pid_particle_pair[1].position
        radius = pid_particle_pair[1].radius
        shell = self.createNewShell(position, radius, domain_id)
        shell_id_shell_pair = (shell_id, shell)
        
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id_shell_pair, reactiontypes)

    def greens_function(self):
        gf = FirstPassageGreensFunction(self.getD())
        a = self.getMobilityRadius()
        gf.seta(a)
        return gf

    def createNewShell(self, position, radius, domain_id):
        return SphericalShell(position, radius, domain_id)

    def displacement(self, r):
        return randomVector(r)

    def __str__(self):
        return 'SphericalSingle' + Single.__str__(self)



