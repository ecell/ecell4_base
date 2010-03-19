from _gfrd import *
from greens_function_wrapper import *

class Single(object):
    """There are 2 main types of Singles:
        * NonInteractionSingle
        * InteractionSingle (when the particle is nearby a surface)

    Each type of Single defines a list of coordinates, see coordinate.py. For each 
    coordinate the Green's function is specified.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        self.multiplicity = 1

        self.pid_particle_pair = pid_particle_pair
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.surface = surface

        # Create shell.
        shell = self.createNewShell(pid_particle_pair[1].position, pid_particle_pair[1].radius, domain_id)

        self.shell_list = [(shell_id, shell), ]

        self.eventID = None

        self.domain_id = domain_id

        self.updatek_tot()

    def getD(self):
        return self.pid_particle_pair[1].D
    D = property(getD)

    def get_shell_id(self):
        return self.shell_list[0][0]
    shell_id = property(get_shell_id)

    def get_shell(self):
        return self.shell_list[0][1]
    shell = property(get_shell)

    def get_shell_id_shell_pair(self):
        return self.shell_list[0]
    def set_shell_id_shell_pair(self, value):
        self.shell_list[0] = value
    shell_id_shell_pair = property(get_shell_id_shell_pair, 
                                   set_shell_id_shell_pair)

    def initialize(self, t):
        '''
        Reset the Single.

        Radius (shell size) is shrunken to the actual radius of the particle.
        self.dt is reset to 0.0.  Do not forget to reschedule this Single
        after calling this method.
        '''
        self.dt = 0.0
        self.lastTime = t
        self.eventType = EventType.SINGLE_ESCAPE

    def isReset(self):
        return self.dt == 0.0 and self.eventType == EventType.SINGLE_ESCAPE

    def draw_reaction_time_tuple(self):
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

    def updatek_tot(self):
        self.k_tot = 0

        if not self.reactiontypes:
            return

        for rt in self.reactiontypes:
            self.k_tot += rt.k

    def drawReactionRule(self):
        k_array = [rt.k for rt in self.reactiontypes]
        k_array = numpy.add.accumulate(k_array)
        k_max = k_array[-1]

        rnd = myrandom.uniform()
        i = numpy.searchsorted(k_array, rnd * k_max)

        return self.reactiontypes[i]

    def check(self):
        pass

    def __str__(self):
        return 'Single[%s: %s: eventID=%s]' % (self.domain_id, self.pid_particle_pair[0], self.eventID)


class NonInteractionSingle(Single):
    """1 Particle inside a shell, no other particles around. 

    There are 3 types of NonInteractionSingles:
        * SphericalSingle: spherical shell, 3D movement.
        * PlanarSurfaceSingle: cylindrical shell, 2D movement.
        * CylindricalSurfaceSingle: cylindrical shell, 1D movement.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        Single.__init__(self, domain_id, pid_particle_pair, shell_id,
                        reactiontypes, surface)

    def get_mobility_radius(self):
        return self.shell_list[0][1].shape.radius - self.pid_particle_pair[1].radius

    def get_shell_size(self):
        return self.shell_list[0][1].shape.radius

    def drawNewPosition(self, dt, eventType):
        gf = self.greens_function()
        a = self.get_mobility_radius()
        r = draw_displacement_wrapper(gf, dt, eventType, a)
        displacement = self.displacement(r)
        if __debug__:
            scale = self.pid_particle_pair[1].radius
            if feq(length(displacement), abs(r), typical=scale) == False:
                raise AssertionError('displacement != abs(r): %g != %g.' % 
                                     (length(displacement), abs(r)))
        return self.pid_particle_pair[1].position + displacement


class SphericalSingle(NonInteractionSingle):
    """1 Particle inside a (spherical) shell not on any surface.

        * Particle coordinate inside shell: r, theta, phi.
        * Coordinate: radial r.
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector: theta, phi.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactiontypes, surface)

    def greens_function(self):
        gf = FirstPassageGreensFunction(self.getD())
        a = self.get_mobility_radius()
        gf.seta(a)
        return gf

    def createNewShell(self, position, radius, domain_id):
        return SphericalShell(domain_id, Sphere(position, radius))

    def displacement(self, r):
        return randomVector(r)

    def __str__(self):
        return 'Spherical' + Single.__str__(self)


class PlanarSurfaceSingle(NonInteractionSingle):
    """1 Particle inside a (cylindrical) shell on a PlanarSurface. (Hockey 
    pucks).

        * Particle coordinates on surface: x, y.
        * Domain: radial r. (determines x and y together with theta).
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector: theta.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactiontypes, surface)

    def greens_function(self):
        # Todo. 2D gf Abs Sym.
        #gf = FirstPassageGreensFunction2D(self.getD())

        gf = FirstPassageGreensFunction(self.getD())
        a = self.get_mobility_radius()
        gf.seta(a)
        return gf

    def createNewShell(self, position, radius, domain_id):
        # The size (thickness) of a hockey puck is not more than it has to be 
        # (namely the radius of the particle), so if the particle undergoes an 
        # unbinding reaction we still have to clear the target volume and the 
        # move may be rejected (NoSpace error).
        orientation = self.surface.shape.unit_z
        size = self.pid_particle_pair[1].radius
        return CylindricalShell(domain_id, Cylinder(position, radius, orientation, size))

    def displacement(self, r):
        x, y = randomVector2D(r)
        return x * self.surface.shape.unit_x + y * self.surface.shape.unit_y

    def __str__(self):
        return 'PlanarSurface' + Single.__str__(self)


class CylindricalSurfaceSingle(NonInteractionSingle):
    """1 Particle inside a (cylindrical) shell on a CylindricalSurface. 
    (Rods).

        * Particle coordinates on surface: z.
        * Domain: cartesian z.
        * Initial position: z = 0.
        * Selected randomly when drawing displacement vector: none.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactiontypes, surface)

    def greens_function(self):
        # Todo. 1D gf Abs Abs.
        #gf = FirstPassageGreensFunction1D(self.getD())
        gf = FirstPassageGreensFunction(self.getD())
        a = self.get_mobility_radius()
        gf.seta(a)
        return gf

    def createNewShell(self, position, size, domain_id):
        # Heads up. The cylinder's *size*, not radius, is changed when you 
        # make the cylinder bigger, because of the redefinition of setRadius.

        # The radius of a rod is not more than it has to be (namely the radius 
        # of the particle), so if the particle undergoes an unbinding reaction 
        # we still have to clear the target volume and the move may be 
        # rejected (NoSpace error).
        radius = self.pid_particle_pair[1].radius
        orientation = self.surface.shape.unit_z
        return CylindricalShell(domain_id, Cylinder(position, radius, orientation, size))

    def displacement(self, z):
        # z can be pos or min.
        return z * self.shell_list[0][1].shape.unit_z

    def get_mobility_radius(self):
        # Heads up.
        return self.shell_list[0][1].shape.size - self.pid_particle_pair[1].radius

    def get_shell_size(self):
        # Heads up.
        return self.shell_list[0][1].shape.size

    def __str__(self):
        return 'CylindricalSurface' + Single.__str__(self)


