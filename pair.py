from _gfrd import *
from constants import EventType
from _greens_functions import *
from greens_function_wrapper import *

__all__ = [
    'CylindricalSurfacePair',
    'PlanarSurfacePair',
    'SphericalPair',
    'Pair',
    ]

class Pair(object):
    """There are 3 types of pairs:
        * SphericalPair
        * PlanarSurfacePair
        * CylindricalSurfacePair

    """
    # CUTOFF_FACTOR is a threshold to choose between the real and 
    # approximate Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__(self, domain_id, com, single1, single2, shell_id, r0, 
                 shell_size, rt, surface):
        self.multiplicity = 2
        self.num_shells = 1

        self.single1 = single1
        self.single2 = single2 

        self.a_R, self.a_r = self.determine_radii(r0, shell_size)

        self.rt = rt

        self.event_id = None

        self.last_time = 0.0
        self.dt = 0.0
        self.event_type = None

        self.surface = surface

        # Create shell.
        shell = self.create_new_shell(com, shell_size, domain_id)

        self.shell_list = [(shell_id, shell), ]
        self.domain_id = domain_id

    def __del__(self):
        if __debug__:
            log.debug('del %s' % str(self))

    def get_com(self):
        return self.shell_list[0][1].shape.position
    com = property(get_com)

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

    def get_shell_size(self):
        return self.shell_list[0][1].shape.radius

    def get_D_tot(self):
        return self.single1.pid_particle_pair[1].D + \
               self.single2.pid_particle_pair[1].D
    D_tot = property(get_D_tot)

    def get_D_R(self):
        return (self.single1.pid_particle_pair[1].D *
                self.single2.pid_particle_pair[1].D) / self.D_tot
    D_R = property(get_D_R)

    def get_v_tot(self):
        return self.single2.pid_particle_pair[1].v - \
               self.single1.pid_particle_pair[1].v
    v_tot = property(get_v_tot)

    def get_v_R(self):
        return (self.single1.pid_particle_pair[1].v * 
                self.single2.pid_particle_pair[1].D +
                self.single2.pid_particle_pair[1].v *
                self.single1.pid_particle_pair[1].D) / self.D_tot
    v_R = property(get_v_R)

    def getSigma(self):
        return self.single1.pid_particle_pair[1].radius + \
               self.single2.pid_particle_pair[1].radius
    sigma = property(getSigma)

    def initialize(self, t):
        self.last_time = t
        self.dt = 0
        self.event_type = None

    def determine_radii(self, r0, shell_size):
        """Determine a_r and a_R.

        Todo. Make dimension (1D/2D/3D) specific someday. Optimization only.

        """
        single1 = self.single1
        single2 = self.single2
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius

        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D

        # Make sure that D1 != 0 to avoid division by zero in the followings.
        if D1 == 0:
            D1, D2 = D2, D1

        shell_size /= SAFETY

        D_tot = D1 + D2
        D_geom = math.sqrt(D1 * D2)

        assert r0 >= self.sigma, \
            '%s;  r0 %g < sigma %g' % (self, r0, self.sigma)

        # equalize expected mean t_r and t_R.
        if ((D_geom - D2) * r0) / D_tot + shell_size +\
                math.sqrt(D2 / D1) * (radius1 - shell_size) - radius2 >= 0:
            Da = D1
            Db = D2
            radiusa = radius1
            radiusb = radius2
        else:
            Da = D2
            Db = D1
            radiusa = radius2
            radiusb = radius1


        #aR
        a_R = (D_geom * (Db * (shell_size - radiusa) + \
                         Da * (shell_size - r0 - radiusa))) /\
              (Da * Da + Da * Db + D_geom * D_tot)

        #ar
        a_r = (D_geom * r0 + D_tot * (shell_size - radiusa)) / (Da + D_geom)

        assert a_R + a_r * Da / D_tot + radius1 >= \
               a_R + a_r * Db / D_tot + radius2

        assert abs(a_R + a_r * Da / D_tot + radiusa - shell_size) \
            < 1e-12 * shell_size


        if __debug__:
          log.debug('a %g, r %g, R %g r0 %g' % 
                 (shell_size, a_r, a_R, r0))
        if __debug__:
            tr = ((a_r - r0) / math.sqrt(6 * self.D_tot))**2
            if self.D_R == 0:
                tR = numpy.inf 
            else:
                tR = (a_R / math.sqrt(6*self.D_R))**2
            log.debug('tr %g, tR %g' % (tr, tR))


        assert a_r > 0
        assert a_r > r0, '%g %g' % (a_r, r0)
        assert a_R > 0 or (a_R == 0 and (D1 == 0 or D2 == 0))

        return a_R, a_r

    def draw_com_escape_or_iv_event_time_tuple(self, r0):
        """Returns a (event time, event type, reactingsingle=None) tuple.
        
        """
        dt_com = draw_time_wrapper(self.com_greens_function())
        dt_iv = draw_time_wrapper(self.iv_greens_function(r0))
        if dt_com < dt_iv:
            return dt_com, EventType.COM_ESCAPE, None
        else:
            # Note: we are not calling pair.draw_iv_event_type yet, but 
            # postpone it to the very last minute (when this event is 
            # executed in fire_pair). So IV_EVENT can still be a iv 
            # event or a iv reaction.
            return dt_iv, EventType.IV_EVENT, None

    def draw_single_reaction_time_tuple(self):
        """Return a (reaction time, event type, reactingsingle)-tuple.

        """
        dt_reaction1, event_type1 = self.single1.draw_reaction_time_tuple()
        dt_reaction2, event_type2 = self.single2.draw_reaction_time_tuple()
        if dt_reaction1 < dt_reaction2:
            return dt_reaction1, event_type1, self.single1
        else:
            return dt_reaction2, event_type2, self.single2

    def determine_next_event(self, r0):
        """Return a (event time, event type, reactingsingle)-tuple.

        """
        return min(self.draw_com_escape_or_iv_event_time_tuple(r0), 
                   self.draw_single_reaction_time_tuple()) 

    def draw_iv_event_type(self, r0):
        gf = self.iv_greens_function(r0)
        event_kind = draw_event_type_wrapper(gf, self.dt)
        if event_kind == PairEventKind.IV_REACTION:
            return EventType.IV_REACTION
        elif event_kind == PairEventKind.IV_ESCAPE:
            return EventType.IV_ESCAPE
        raise NotImplemented()

    def draw_new_positions(self, dt, r0, old_iv, event_type):
        """Calculate new positions of the pair particles using a new 
        center-of-mass, a new inter-particle vector, and an old 
        inter-particle vector.

        """
        new_com = self.draw_new_com(dt, event_type)
        new_iv = self.draw_new_iv(dt, r0, old_iv, event_type)

        D1 = self.single1.pid_particle_pair[1].D
        D2 = self.single2.pid_particle_pair[1].D

        newpos1 = new_com - new_iv * (D1 / self.D_tot)
        newpos2 = new_com + new_iv * (D2 / self.D_tot)
        return newpos1, newpos2

    def draw_new_com(self, dt, event_type):
        if event_type == EventType.COM_ESCAPE:
            r = self.a_R
        else:
            gf = self.com_greens_function()
            r = draw_r_wrapper(gf, dt, self.a_R)

        displacement = self.create_com_vector(r)

        # Add displacement to old CoM. This assumes (correctly) that 
        # r0=0 for the CoM. Compare this to 1D singles, where r0 is not  
        # necesseraly 0.
        return self.com + displacement

    def draw_new_iv(self, dt, r0, old_iv, event_type):
        gf = self.choose_pair_greens_function(r0, dt)
        if event_type == EventType.IV_ESCAPE:
            r = self.a_r
        elif event_type == EventType.IV_REACTION:
            r = self.sigma
        else:
            r = draw_r_wrapper(gf, dt, self.a_r, self.sigma)

        return self.create_interparticle_vector(gf, r, dt, r0, old_iv)

    def check(self):
        pass

    def __str__(self):
        sid = self.single1.pid_particle_pair[1].sid
        name = self.world.model.get_species_type_by_id(sid)["name"]
        if name[0] != '(':
            name = '(' + name + ')'
        return 'Pair[%s: %s, %s, %s]' % (
            self.domain_id,
            self.single1.pid_particle_pair[0],
            self.single2.pid_particle_pair[0],
            name)

class SphericalPair(Pair):
    """2 Particles inside a (spherical) shell not on any surface.

    """
    def __init__(self, domain_id, com, single1, single2, shell_id,
                 r0, shell_size, rt, surface):
        Pair.__init__(self, domain_id, com, single1, single2, shell_id,
                      r0, shell_size, rt, surface)

    def com_greens_function(self):
        # Green's function for centre of mass inside absorbing sphere.
        return GreensFunction3DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # Green's function for interparticle vector inside absorbing 
        # sphere.  This exact solution is used for drawing times.
        return GreensFunction3DRadAbs(self.D_tot, self.rt.ktot, r0,
                                              self.sigma, self.a_r)

    def create_new_shell(self, position, radius, domain_id):
        return SphericalShell(domain_id, Sphere(position, radius))

    def choose_pair_greens_function(self, r0, t):
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(6.0 * self.D_tot * t)

        if distance_from_sigma < threshold_distance:
        
            if distance_from_shell < threshold_distance:
                # near both a and sigma;
                # use GreensFunction3DRadAbs
                if __debug__:
                    log.debug('GF: normal')
                return self.iv_greens_function(r0)
            else:
                # near sigma; use GreensFunction3DRadInf
                if __debug__:
                    log.debug('GF: only sigma')
                return GreensFunction3DRadInf(self.D_tot, self.rt.ktot, r0,
                                               self.sigma)
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF: only a')
                return GreensFunction3DAbs(self.D_tot,
                                                                 r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF: free')
                return GreensFunction3D(self.D_tot, r0)

    def create_com_vector(self, r):
        return random_vector(r)

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        theta = draw_theta_wrapper(gf, r, dt)

        new_inter_particle_s = numpy.array([r, theta, 
                                         myrandom.uniform() * 2 * Pi])
        new_iv = spherical_to_cartesian(new_inter_particle_s)

        #FIXME: need better handling of angles near zero and pi.

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotation_axis = crossproduct([0,0,1], inter_particle)
        angle = vector_angle_against_z_axis(old_iv)
        if angle % numpy.pi != 0.0:
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            rotated = rotate_vector(new_iv, rotation_axis, angle)
        elif angle == 0.0:
            rotated = new_iv
        else:
            rotated = numpy.array([new_iv[0], new_iv[1], - new_iv[2]])
        return rotated

    def __str__(self):
        return 'Spherical' + Pair.__str__(self)


class PlanarSurfacePair(Pair):
    """2 Particles inside a (cylindrical) shell on a PlanarSurface. 
    (Hockey pucks).

    """
    def __init__(self, domain_id, com, single1, single2, shell_id,
                 r0, shell_size, rt, surface):
        Pair.__init__(self, domain_id, com, single1, single2, shell_id,
                      r0, shell_size, rt, surface)

    def com_greens_function(self):
        # Todo. 2D gf Abs Sym.
        return GreensFunction3DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # Todo. 2D gf Rad Abs.
        # This exact solution is used for drawing times.
        return GreensFunction3DRadAbs(self.D_tot, self.rt.ktot, r0,
                                              self.sigma, self.a_r)

    def create_new_shell(self, position, radius, domain_id):
        # The size (thickness) of a hockey puck is not more than it has 
        # to be (namely the radius of the particle), so if the particle 
        # undergoes an unbinding reaction we still have to clear the 
        # target volume and the move may be rejected (NoSpace error).
        orientation = crossproduct(self.surface.shape.unit_x,
                                   self.surface.shape.unit_y)
        size = max(self.single1.pid_particle_pair[1].radius,
                   self.single2.pid_particle_pair[1].radius)
        return CylindricalShell(domain_id, Cylinder(position, radius, 
                                                    orientation, size))

        a_R, a_r = self.determine_radii()

    def choose_pair_greens_function(self, r0, t):
        # Todo
        return self.iv_greens_function(r0)

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        return x * self.surface.shape.unit_x + y * self.surface.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        theta = draw_theta_wrapper(gf, r, dt)

        #FIXME: need better handling of angles near zero and pi?
        unit_x = self.surface.shape.unit_x
        unit_y = self.surface.shape.unit_y
        angle = vector_angle(unit_x, old_iv)
        # Todo. Test if nothing changes when theta == 0.
        new_angle = angle + theta

        new_iv = r * math.cos(new_angle) * unit_x + \
                 r * math.sin(new_angle) * unit_y

        return new_iv

    def __str__(self):
        return 'PlanarSurface' + Pair.__str__(self)


class CylindricalSurfacePair(Pair):
    """2 Particles inside a (cylindrical) shell on a CylindricalSurface.  
    (Rods).

    """
    def __init__(self, domain_id, com, single1, single2, shell_id,
                 r0, shell_size, rt, surface):
        Pair.__init__(self, domain_id, com, single1, single2, shell_id,
                      r0, shell_size, rt, surface)

    def com_greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        return GreensFunction1DAbsAbs(self.D_R, self.v_R, 0.0, -self.a_R, self.a_R)

    def iv_greens_function(self, r0):
        return GreensFunction1DRadAbs(self.D_tot, self.v_r, self.rt.ktot, r0, self.sigma, self.a_r)

    def create_new_shell(self, position, size, domain_id):
        # The radius of a rod is not more than it has to be (namely the 
        # radius of the biggest particle), so if the particle undergoes 
        # an unbinding reaction we still have to clear the target volume 
        # and the move may be rejected (NoSpace error).
        radius = max(self.single1.pid_particle_pair[1].radius,
                     self.single2.pid_particle_pair[1].radius)
        orientation = self.surface.shape.unit_z
        return CylindricalShell(domain_id,
                                Cylinder(position, radius, orientation, size))

    def choose_pair_greens_function(self, r0, t):
        # Todo
        return self.iv_greens_function(r0)

    def create_com_vector(self, r):
        return r * self.surface.shape.unit_z

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        # Note: using self.surface.shape.unit_z here might accidently 
        # interchange the particles.
        return r * normalize(old_iv)

    def get_shell_size(self):
        # Heads up.
        return self.shell_list[0][1].shape.size

    def __str__(self):
        return 'CylindricalSurface' + Pair.__str__(self)




