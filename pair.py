from _gfrd import *
from greens_function_wrapper import *

class Pair( object ):
    """There are 3 types of pairs:
        * SphericalPair
        * PlanarSurfacePair
        * CylindricalSurfacePair

    """
    # CUTOFF_FACTOR is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__(self, domain_id, CoM, single1, single2, shell_id, r0, 
                 shellSize, rt, surface):
        self.multiplicity = 2

        # Order single1 and single2 so that D1 < D2.
        if single1.pid_particle_pair[1].D <= single2.pid_particle_pair[1].D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.a_R, self.a_r = self.determineRadii(r0, shellSize)

        self.rt = rt

        self.eventID = None

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.surface = surface

        # Create shell.
        shell = self.createNewShell(CoM, shellSize, domain_id)

        self.shell_list = [(shell_id, shell), ]
        self.domain_id = domain_id

    def __del__( self ):
        if __debug__:
            log.debug( 'del %s' % str( self ) )

    def getCoM(self):
        return self.shell_list[0][1].shape.position
    CoM = property(getCoM)

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

    def get_D_tot( self ):
        return self.single1.pid_particle_pair[1].D + \
               self.single2.pid_particle_pair[1].D
    D_tot = property(get_D_tot)

    def get_D_R(self):
        return (self.single1.pid_particle_pair[1].D *
                self.single2.pid_particle_pair[1].D) / self.D_tot
    D_R = property(get_D_R)

    def get_sigma(self):
        return self.single1.pid_particle_pair[1].radius + \
               self.single2.pid_particle_pair[1].radius
    sigma = property(get_sigma)

    def initialize( self, t ):
        self.lastTime = t
        self.dt = 0
        self.eventType = None

    def determineRadii(self, r0, shellSize):
        """Determine a_r and a_R.

        Todo. Make dimension (1D/2D/3D) specific someday. Optimization only.

        """
        single1 = self.single1
        single2 = self.single2
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius

        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D

        shellSize /= SAFETY  # FIXME:

        D_tot = D1 + D2
        D_geom = math.sqrt(D1 * D2)

        assert r0 >= self.sigma, \
            '%s;  r0 %g < sigma %g' % ( self, r0, self.sigma )

        # equalize expected mean t_r and t_R.
        if ((D_geom - D2) * r0) / D_tot + shellSize +\
                math.sqrt(D2 / D1) * (radius1 - shellSize) - radius2 >= 0:
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
        a_R = (D_geom * (Db * (shellSize - radiusa) + \
                         Da * (shellSize - r0 - radiusa))) /\
              (Da * Da + Da * Db + D_geom * D_tot)

        #ar
        a_r = (D_geom * r0 + D_tot * (shellSize - radiusa)) / (Da + D_geom)

        assert a_R + a_r * Da / D_tot + radius1 >= \
               a_R + a_r * Db / D_tot + radius2

        assert abs(a_R + a_r * Da / D_tot + radiusa - shellSize ) \
            < 1e-12 * shellSize


        if __debug__:
          log.debug( 'a %g, r %g, R %g r0 %g' % 
                  ( shellSize, a_r, a_R, r0 ) )
        if __debug__:
          log.debug( 'tr %g, tR %g' % 
                     ( ( ( a_r - r0 ) / math.sqrt(6 * self.D_tot))**2,\
                           (a_R / math.sqrt( 6*self.D_R ))**2 ) )
        assert a_r > 0
        assert a_r > r0, '%g %g' % ( a_r, r0 )
        assert a_R > 0 or ( a_R == 0 and ( D1 == 0 or D2 == 0 ) )

        return a_R, a_r

    def draw_com_escape_or_iv_event_time_tuple(self, r0):
        """Returns a (event time, event type, reactingsingle=None) tuple.
        
        """
        dt_com = draw_time_wrapper(self.com_greens_function())
        dt_iv = draw_time_wrapper(self.iv_greens_function(), r0)
        if dt_com < dt_iv:
            return dt_com, EventType.COM_ESCAPE, None
        else:
            # Note: we are not calling pair.draw_iv_event_type yet, but 
            # postpone it to the very last minute (when this event is executed 
            # in firePair). So IV_EVENT can still be a iv event or a iv 
            # reaction.
            return dt_iv, EventType.IV_EVENT, None

    def draw_single_reaction_time_tuple( self ):
        """Return a (reaction time, event type, reactingsingle)-tuple.

        """
        dt_reaction1, event_type1 = self.single1.draw_reaction_time_tuple()
        dt_reaction2, event_type2 = self.single2.draw_reaction_time_tuple()
        if dt_reaction1 < dt_reaction2:
            return dt_reaction1, event_type1, self.single1
        else:
            return dt_reaction2, event_type2, self.single2

    def determineNextEvent(self, r0):
        """Return a (event time, event type, reactingsingle)-tuple.

        """
        return min(self.draw_com_escape_or_iv_event_time_tuple(r0), 
                   self.draw_single_reaction_time_tuple()) 

    def draw_iv_event_type(self, r0):
        gf = self.iv_greens_function()
        return draw_eventtype_wrapper(gf, self.dt, r0)

    def drawNewPositions(self, dt, r0, old_iv, eventType):
        """Calculate new positions of the pair particles using a new 
        center-of-mass, a new inter-particle vector, and an old inter-particle 
        vector.

        """
        new_com = self.drawNewCoM(dt, eventType)
        new_iv = self.drawNewIV(dt, r0, old_iv, eventType)

        D1 = self.single1.pid_particle_pair[1].D
        D2 = self.single2.pid_particle_pair[1].D

        newpos1 = new_com - new_iv * (D1 / self.D_tot)
        newpos2 = new_com + new_iv * (D2 / self.D_tot)
        return newpos1, newpos2

    def check( self ):
        pass

    def __str__(self):
        return 'Pair[%s: %s, %s: eventID=%s]' % (
            self.domain_id,
            self.single1.pid_particle_pair[0],
            self.single2.pid_particle_pair[0],
            self.eventID)

class SphericalPair(Pair):
    """2 Particles inside a (spherical) shell not on any surface.

    """
    def __init__(self, domain_id, CoM, single1, single2, shell_id,
                 r0, shellSize, rt, surface):
        Pair.__init__(self, domain_id, CoM, single1, single2, shell_id,
                      r0, shellSize, rt, surface)

    def com_greens_function(self):
        # Green's function for centre of mass inside absorbing sphere.
        gf = FirstPassageGreensFunction(self.D_R)
        gf.seta(self.a_R)
        return gf

    def iv_greens_function(self):
        # Green's function for interparticle vector inside absorbing sphere.  
        # This exact solution is used for drawing times.
        gf = FirstPassagePairGreensFunction(self.D_tot, self.rt.k, self.sigma)
        gf.seta(self.a_r)
        return gf

    def createNewShell(self, position, radius, domain_id):
        return SphericalShell(domain_id, Sphere(position, radius))

    def choosePairGreensFunction( self, r0, t ):
        distanceFromSigma = r0 - self.sigma
        distanceFromShell = self.a_r - r0

        thresholdDistance = Pair.CUTOFF_FACTOR * \
            math.sqrt( 6.0 * self.D_tot * t )

        if distanceFromSigma < thresholdDistance:
        
            if distanceFromShell < thresholdDistance:
                # near both a and sigma;
                # use FirstPassagePairGreensFunction
                if __debug__:
                    log.debug( 'GF: normal' )
                return self.iv_greens_function()
            else:
                # near sigma; use BasicPairGreensFunction
                if __debug__:
                    log.debug( 'GF: only sigma' )
                pgf = BasicPairGreensFunction( self.D_tot, self.rt.k, 
                                               self.sigma )
                return pgf
        else:
            if distanceFromShell < thresholdDistance:
                # near a;
                if __debug__:
                    log.debug( 'GF: only a' )
                pgf = FirstPassageNoCollisionPairGreensFunction( self.D_tot )
                pgf.seta(self.a_r)
                return pgf
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug( 'GF: free' )
                pgf = FreePairGreensFunction( self.D_tot )
                return pgf

    def drawNewCoM(self, dt, eventType):
        gf = self.com_greens_function()
        r_R = draw_displacement_wrapper(gf, dt, eventType, self.a_R)
        return self.CoM + randomVector(r_R)

    def drawNewIV(self, dt, r0, old_iv, eventType): 
        gf = self.choosePairGreensFunction(r0, dt)
        r, theta = draw_displacement_iv_wrapper(gf, r0, dt, eventType,
                                                self.a_r, self.sigma)
        newInterParticleS = numpy.array([r, theta, 
                                         myrandom.uniform() * 2 * Pi])
        new_iv = sphericalToCartesian(newInterParticleS)

        #FIXME: need better handling of angles near zero and pi.

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotationAxis = crossproduct( [ 0,0,1 ], interParticle )
        angle = vectorAngleAgainstZAxis(old_iv)
        if angle % numpy.pi != 0.0:
            rotationAxis = crossproductAgainstZAxis(old_iv)
            rotationAxis = normalize(rotationAxis)
            rotated = rotateVector(new_iv, rotationAxis, angle)
        elif angle == 0.0:
            rotated = new_iv
        else:
            rotated = numpy.array([new_iv[0], new_iv[1], - new_iv[2]])
        return rotated

    def __str__(self):
        return 'Spherical' + Pair.__str__(self)


class PlanarSurfacePair( Pair ):
    """2 Particles inside a (cylindrical) shell on a PlanarSurface. (Hockey 
    pucks).

    """
    def __init__(self, domain_id, CoM, single1, single2, shell_id,
                 r0, shellSize, rt, surface):
        Pair.__init__(self, domain_id, CoM, single1, single2, shell_id,
                      r0, shellSize, rt, surface)

    def com_greens_function(self):
        # Todo. 2D gf Abs Sym.
        gf = FirstPassageGreensFunction(self.D_R)
        gf.seta(self.a_R)
        return gf

    def iv_greens_function(self):
        # Todo. 2D gf Rad Abs.
        # This exact solution is used for drawing times.
        gf = FirstPassagePairGreensFunction(self.D_tot, self.rt.k, self.sigma)
        gf.seta(self.a_r)
        return gf

    def createNewShell(self, position, radius, domain_id):
        # The size (thickness) of a hockey puck is not more than it has to be 
        # (namely the radius of the particle), so if the particle undergoes an 
        # unbinding reaction we still have to clear the target volume and the 
        # move may be rejected (NoSpace error).
        orientation = self.surface.shape.unit_z
        size = max(self.single1.pid_particle_pair[1].radius,
                   self.single2.pid_particle_pair[1].radius)
        return CylindricalShell(position, radius, orientation, size, domain_id)

        a_R, a_r = self.determineRadii()

    def drawNewCoM(self, dt, eventType):
        gf = self.com_greens_function()
        r_R = draw_displacement_wrapper(gf, dt, eventType, self.a_R)
        x, y = randomVector2D(r_R)
        return(self.CoM + x * self.surface.shape.unit_x
                        + y * self.surface.shape.unit_y)

    def drawNewIV(self, dt, r0, old_iv, eventType): 
        # Todo.
        #gf = self.choosePairGreensFunction(r0, dt)
        gf = self.iv_greens_function()
        r, theta = draw_displacement_iv_wrapper(gf, r0, dt, eventType,
                                                self.a_r, self.sigma)
        assert r > self.sigma and r <= self.a_r

        unitX = self.surface.shape.unit_x
        unitY = self.surface.shape.unit_y
        angle = vectorAngle(unitX, old_iv)
        # Todo. Test if nothing changes when theta == 0.
        new_angle = angle + theta

        new_iv = r * math.cos(new_angle) * unitX + \
                 r * math.sin(new_angle) * unitY

        return new_iv

    def __str__(self):
        return 'PlanarSurface' + Pair.__str__(self)


class CylindricalSurfacePair( Pair ):
    """2 Particles inside a (cylindrical) shell on a CylindricalSurface. 
    (Rods).

    """
    def __init__(self, domain_id, CoM, single1, single2, shell_id,
                 r0, shellSize, rt, surface):
        Pair.__init__(self, domain_id, CoM, single1, single2, shell_id,
                      r0, shellSize, rt, surface)

    def com_greens_function(self):
        # Todo. 1D gf Abs Abs  should be -a to a.
        #gf = FirstPassageGreensFunction1D(self.D_R)
        gf = FirstPassageGreensFunction(self.D_R)
        gf.seta(self.a_R)
        return gf

    def iv_greens_function(self):
        # Todo. 1D gf Rad Abs  should be sigma to a.
        #gf = FirstPassageGreensFunction1DRad(self.D_tot, self.rt.k)
        # This exact solution is used for drawing times.
        gf = FirstPassagePairGreensFunction(self.D_tot, self.rt.k, self.sigma)
        gf.seta(self.a_r)
        return gf

    def createNewShell(self, position, size, domain_id):
        # The radius of a rod is not more than it has to be (namely the radius 
        # of the biggest particle), so if the particle undergoes an unbinding 
        # reaction we still have to clear the target volume and the move may 
        # be rejected (NoSpace error).
        radius = max(self.single1.pid_particle_pair[1].radius,
                     self.single2.pid_particle_pair[1].radius)
        orientation = self.surface.shape.unit_z
        return CylindricalShell(position, radius, orientation, size, domain_id)

    def drawNewCoM(self, dt, eventType):
        gf = self.com_greens_function()
        # Draw displacement (not absolute position).
        r_R = draw_displacement_wrapper(gf, dt, eventType, self.a_R) # Todo.
        return self.CoM + r_R * self.surface.shape.unit_z

    def drawNewIV(self, dt, r0, old_iv, eventType): 
        # Todo.
        #gf = self.choosePairGreensFunction(r0, dt)
        gf = self.iv_greens_function()

        r = draw_displacement_wrapper(gf, dt, eventType, self.a_r, r0, 
                                      self.sigma)
        assert r > self.sigma and r <= self.a_r

        # Note: using self.surface.shape.unit_z here might accidently 
        # interchange the particles.
        return r * normalize(old_iv)

    def get_shell_size(self):
        # Heads up.
        return self.shell_list[0][1].size

    def __str__(self):
        return 'CylindricalSurface' + Pair.__str__(self)



