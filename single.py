from _gfrd import *
from _greens_functions import *
from greens_function_wrapper import *
from constants import EventType
import utils

__all__ = [
    'CylindricalSurfaceSingle',
    'PlanarSurfaceSingle',
    'SphericalSingle',
    'Single',
    ]

class Single(object):
    """There are 2 main types of Singles:
        * NonInteractionSingle
        * InteractionSingle (when the particle is nearby a surface)

    Each type of Single defines a list of coordinates, see 
    coordinate.py. For each coordinate the Green's function is 
    specified.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        self.multiplicity = 1
        self.num_shells = 1

        self.pid_particle_pair = pid_particle_pair
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.last_time = 0.0
        self.dt = 0.0
        self.event_type = None

        self.surface = surface

        # Create shell.
        shell = self.create_new_shell(pid_particle_pair[1].position,
                                      pid_particle_pair[1].radius, domain_id)

        self.shell_list = [(shell_id, shell), ]

        self.event_id = None

        self.domain_id = domain_id

        self.updatek_tot()

    def getD(self):
        return self.pid_particle_pair[1].D
    D = property(getD)

    def getv(self):
        return 0 # TODO.
        return self.pid_particle_pair[1].v
    v = property(getv)

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

        Radius (shell size) is shrunken to the actual radius of the 
        particle.  self.dt is reset to 0.0.  Do not forget to reschedule 
        this Single after calling this method.
        '''
        self.dt = 0.0
        self.last_time = t
        self.event_type = EventType.SINGLE_ESCAPE

    def is_reset(self):
        return self.dt == 0.0 and self.event_type == EventType.SINGLE_ESCAPE

    def draw_reaction_time_tuple(self):
        """Return a (reaction time, event type)-tuple.

        """
        if self.k_tot <= 0:
            dt = numpy.inf
        elif self.k_tot == numpy.inf:
            dt = 0.0
        else:
            rnd = myrandom.uniform()
            if rnd == 0:
                dt = numpy.inf
            else:
                dt = (1.0 / self.k_tot) * (- math.log(rnd)) # log(1/x) == - log(x)
        return dt, EventType.SINGLE_REACTION

    def draw_interaction_time(self):
        """Todo.
        
        Note: we are not calling single.drawEventType() just yet, but 
        postpone it to the very last minute (when this event is executed 
        in fire_single). So IV_EVENT can still be an iv escape or an iv 
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

    def determine_next_event(self):
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

    def draw_reaction_rule(self):
        k_array = [rt.k for rt in self.reactiontypes]
        k_array = numpy.add.accumulate(k_array)
        k_max = k_array[-1]

        rnd = myrandom.uniform()
        i = numpy.searchsorted(k_array, rnd * k_max)

        return self.reactiontypes[i]

    def check(self):
        pass

    def __str__(self):
        pid = self.pid_particle_pair[0]
        sid = self.pid_particle_pair[1].sid
        name = self.world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        return 'Single[%s: %s, ST%s]' % (self.domain_id, pid, name)


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
        return self.get_shell_size() - self.pid_particle_pair[1].radius

    def get_shell_size(self):
        return self.shell_list[0][1].shape.radius

    def draw_new_position(self, dt, event_type):
        if event_type == EventType.SINGLE_ESCAPE:
            # Moving this checks to the Green's functions is not a good 
            # idea, because then you'd draw an unused random number.  
            # The same yields for the draw_new_com and draw_new_iv.  
            r = self.get_mobility_radius()
        else:
            gf = self.greens_function()
            r = draw_r_wrapper(gf, dt, self.get_mobility_radius())

        displacement = self.create_position_vector(r)

        # This should be checked in the unit test of random_vector.
        # if __debug__:
        #     scale = self.pid_particle_pair[1].radius
        #     if feq(length(displacement), abs(r), typical=scale) == False:
        #         raise AssertionError('displacement != abs(r): %g != %g.' % 
        #                              (length(displacement), abs(r)))

        # Add displacement to shape.position, not to particle.position.  
        # This distinction is important only in the case of an 
        # asymmetric 1D domain (r0 != 0, or drift), since draw_r always 
        # returns a position relative to the centre of the shell (r=0), 
        # not relative to r0.
        return self.shell.shape.position + displacement


class SphericalSingle(NonInteractionSingle):
    """1 Particle inside a (spherical) shell not on any surface.

        * Particle coordinate inside shell: r, theta, phi.
        * Coordinate: radial r.
        * Initial position: r = 0.
        * Selected randomly when drawing displacement vector:
          theta, phi.

    """
    def __init__(self, domain_id, pid_particle_pair, shell_id, reactiontypes, 
                 surface):
        NonInteractionSingle.__init__(self, domain_id, pid_particle_pair, 
                                      shell_id, reactiontypes, surface)

    def greens_function(self):
        return GreensFunction3DAbsSym(self.getD(),
                                          self.get_mobility_radius())

    def create_new_shell(self, position, radius, domain_id):
        return SphericalShell(domain_id, Sphere(position, radius))

    def create_position_vector(self, r):
        return random_vector(r)

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
        #gf = GreensFunction2DAbsSym(self.getD())
        return GreensFunction3DAbsSym(self.getD(),
                                          self.get_mobility_radius())

    def create_new_shell(self, position, radius, domain_id):
        # The half_length (thickness) of a hockey puck is not more than 
        # it has to be (namely the radius of the particle), so if the 
        # particle undergoes an unbinding reaction we still have to 
        # clear the target volume and the move may be rejected (NoSpace 
        # error).
        orientation = normalize(
            utils.crossproduct(self.surface.shape.unit_x,
                               self.surface.shape.unit_y))
        half_length = self.pid_particle_pair[1].radius
        return CylindricalShell(domain_id, Cylinder(position, radius, 
                                                    orientation, half_length))

    def create_position_vector(self, r):
        x, y = random_vector2D(r)
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
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        return GreensFunction1DAbsAbs(self.getD(), self.getv(), 0.0, -self.get_mobility_radius(), self.get_mobility_radius())

    def create_new_shell(self, position, half_length, domain_id):
        # The radius of a rod is not more than it has to be (namely the 
        # radius of the particle), so if the particle undergoes an 
        # unbinding reaction we still have to clear the target volume 
        # and the move may be rejected (NoSpace error).
        radius = self.pid_particle_pair[1].radius
        orientation = self.surface.shape.unit_z
        return CylindricalShell(domain_id, Cylinder(position, radius, 
                                                    orientation, half_length))

    def create_position_vector(self, z):
        if utils.feq(z, self.get_mobility_radius()):
            # Escape, can be either to the left or to the right.
            z = myrandom.choice(-1, 1) * z 
        return z * self.shell_list[0][1].shape.unit_z

    def get_shell_size(self):
        # Heads up. The cylinder's *half_length*, not radius, 
        # determines the size in case of a cylindrical surface.
        return self.shell_list[0][1].shape.half_length

    def __str__(self):
        return 'CylindricalSurface' + Single.__str__(self)


