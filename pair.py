
import numpy

from _gfrd import *

from utils import *
import myrandom

import logging

log = logging.getLogger('ecell')

class Pair( object ):
    
    # CUTOFF_FACTOR is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    def __init__(self, domain_id, single1, single2, shell_id_shell_pair, rt):
        self.multiplicity = 2

        # Order single1 and single2 so that D1 < D2.
        if single1.pid_particle_pair[1].D <= single2.pid_particle_pair[1].D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.rt = rt

        D1 = self.single1.pid_particle_pair[1].D
        D2 = self.single2.pid_particle_pair[1].D

        self.D_tot = D1 + D2
        self.D_R = ( D1 * D2 ) / self.D_tot

        self.sigma = self.single1.pid_particle_pair[1].radius + \
                     self.single2.pid_particle_pair[1].radius

        self.eventID = None

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.shell_list = [shell_id_shell_pair, ]
        self.domain_id = domain_id

    def __del__( self ):
        if __debug__:
            log.debug( 'del %s' % str( self ) )

    def getShell(self):
        return self.shell_list[0]

    def setShell(self, value):
        self.shell_list[0] = value

    shell = property(getShell, setShell)

    def initialize( self, t ):

        self.lastTime = t
        self.dt = 0
        self.eventType = None

    def getD( self ):
        return self.D_tot #FIXME: is this correct?

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
        self.a_R = (D_geom * (Db * (shellSize - radiusa) + \
                               Da * (shellSize - r0 - radiusa))) /\
                               (Da * Da + Da * Db + D_geom * D_tot)

        #ar
        self.a_r = (D_geom * r0 + D_tot * (shellSize - radiusa)) /\
            (Da + D_geom)

        assert self.a_R + self.a_r * Da / D_tot + radius1 >= \
            self.a_R + self.a_r * Db / D_tot + radius2

        assert abs( self.a_R + self.a_r * Da / D_tot + radiusa - shellSize ) \
            < 1e-12 * shellSize


        if __debug__:
          log.debug( 'a %g, r %g, R %g r0 %g' % 
                  ( shellSize, self.a_r, self.a_R, r0 ) )
        if __debug__:
          log.debug( 'tr %g, tR %g' % 
                     ( ( ( self.a_r - r0 ) / math.sqrt(6 * self.D_tot))**2,\
                           (self.a_R / math.sqrt( 6*self.D_R ))**2 ) )
        assert self.a_r > 0
        assert self.a_r > r0, '%g %g' % ( self.a_r, r0 )
        assert self.a_R > 0 or ( self.a_R == 0 and ( D1 == 0 or D2 == 0 ) )

    def determinePairEvent(self, t, r0, shellSize):
        self.lastTime = t

        self.determineRadii(r0, shellSize)

        sgf = FirstPassageGreensFunction(self.D_R)
        sgf.seta(self.a_R)


        # draw t_R
        try:
            self.t_R = self.drawTime_single( sgf )
        except Exception, e:
            raise Exception, 'sgf.drawTime() failed; %s; %s' %\
                ( str( e ), sgf.dump() )

        pgf = FirstPassagePairGreensFunction( self.D_tot, 
                                              self.rt.k, self.sigma )
        pgf.seta( self.a_r )

        # draw t_r
        try:
            self.t_r = self.drawTime_pair(pgf, r0)
        except Exception, e:
            raise Exception, \
                'pgf.drawTime() failed; %s; r0=%g, %s' % \
                ( str( e ), r0, pgf.dump() )


        # draw t_reaction
        t_reaction1 = self.single1.drawReactionTime()
        t_reaction2 = self.single2.drawReactionTime()

        if t_reaction1 < t_reaction2:
            self.t_single_reaction = t_reaction1
            self.reactingsingle = self.single1
        else:
            self.t_single_reaction = t_reaction2
            self.reactingsingle = self.single2

        self.dt = min( self.t_R, self.t_r, self.t_single_reaction )

        assert self.dt >= 0
        if __debug__:
            log.debug( 'dt %g, t_R %g, t_r %g' % 
                     ( self.dt, self.t_R, self.t_r ) )

        if self.dt == self.t_r:  # type = 0 (REACTION) or 1 (ESCAPE_r)
            try:
                self.eventType = self.drawEventType(pgf, r0, self.t_r)
            except Exception, e:
                raise Exception,\
                    'pgf.drawEventType() failed; %s; r0=%g, %s' %\
                    ( str( e ), r0, pgf.dump() )

        elif self.dt == self.t_R: # type = ESCAPE_R (2)
            self.eventType = 2
        elif self.dt == self.t_single_reaction:  # type = single reaction (3)
            self.eventType = 3 
        else:
            raise AssertionError, "Never get here"

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
                pgf = FirstPassagePairGreensFunction( self.D_tot, 
                                                      self.rt.k, self.sigma )

                return pgf
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
                return pgf
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug( 'GF: free' )
                pgf = FreePairGreensFunction( self.D_tot )
                return pgf

    def drawTime_single( self, sgf ):
        rnd = myrandom.uniform()
        return sgf.drawTime( rnd )

    def drawTime_pair( self, pgf, r0 ):
        rnd = myrandom.uniform()
        #print 'r0 = ', r0, ', rnd = ', rnd[1],\
        #    pgf.dump()
        return pgf.drawTime( rnd, r0 )

    def drawEventType( self, pgf, r0, t ):
        rnd = myrandom.uniform()
        return pgf.drawEventType( rnd, r0, t )

    def drawR_single( self, sgf, t ):
        rnd = myrandom.uniform()
        try:
            r = sgf.drawR( rnd, t )
            while r > self.a_R: # redraw; shouldn't happen often
                if __debug__:
                    log.info( 'drawR_single: redraw' )
                rnd = myrandom.uniform()
                r = sgf.drawR( rnd, t )
        except Exception, e:
            raise Exception,\
                'gf.drawR_single() failed; %s; rnd= %g, t= %g, %s' %\
                ( str( e ), rnd, t, sgf.dump() )

        return r

    def drawR_pair( self, r0, t, a ):
        '''
        Draw r for the pair inter-particle vector.
        '''
        gf = self.choosePairGreensFunction( r0, t )

        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        rnd = myrandom.uniform()
        try:
            r = gf.drawR( rnd, r0, t )
            # redraw; shouldn't happen often
            while r >= self.a_r or r <= self.sigma: 
                if __debug__:
                    log.info( 'drawR_pair: redraw' )
                #self.sim.rejectedMoves += 1  #FIXME:
                rnd = myrandom.uniform()
                r = gf.drawR( rnd, r0, t )
        except Exception, e:
            raise Exception,\
                'gf.drawR_pair() failed; %s; rnd= %g, r0= %g, t= %g, %s' %\
                ( str( e ), rnd, r0, t, gf.dump() )


        return r

    def drawTheta_pair( self, rnd, r, r0, t, a ):
        '''
        Draw theta for the pair inter-particle vector.
        '''
        gf = self.choosePairGreensFunction( r0, t )

        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        try:
            theta = gf.drawTheta( rnd, r, r0, t )
        except Exception, e:
            raise Exception,\
                'gf.drawTheta() failed; %s; rnd= %g, r= %g, r0= %g, t=%g, %s' %\
                ( str( e ), rnd, r, r0, t, gf.dump() )

        return theta

    def check( self ):
        pass

    def __repr__( self ):
        return 'Pair[%s: %s, %s: eventID=%s]' % (
            self.domain_id,
            self.single1.pid_particle_pair[0],
            self.single2.pid_particle_pair[0],
            self.eventID )

    def calculatePairPos(self, CoM, newInterParticle, oldInterParticle):
        '''
        Calculate new positions of the particles in the Pair using
        a new center-of-mass, a new inter-particle vector, and
        an old inter-particle vector.
        '''
        #FIXME: need better handling of angles near zero and pi.

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotationAxis = crossproduct( [ 0,0,1 ], interParticle )

        angle = vectorAngleAgainstZAxis(oldInterParticle)
        if angle % numpy.pi != 0.0:
            rotationAxis = crossproductAgainstZAxis(oldInterParticle)
            rotationAxis = normalize(rotationAxis)
            rotated = rotateVector(newInterParticle, rotationAxis, angle)
        elif angle == 0.0:
            rotated = newInterParticle
        else:
            rotated = numpy.array([newInterParticle[0], newInterParticle[1],
                                   - newInterParticle[2]])

        D1 = self.single1.pid_particle_pair[1].D
        D2 = self.single2.pid_particle_pair[1].D
        newpos1 = CoM - rotated * (D1 / (D1 + D2))
        newpos2 = CoM + rotated * (D2 / (D1 + D2))

        return newpos1, newpos2

