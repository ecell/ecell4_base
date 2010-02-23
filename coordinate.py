import numpy
import myrandom
from _gfrd import EventType #, FirstPassageGreensFunction1D (TODO)
from utils import *

import logging
log = logging.getLogger('ecell')


'''
                                                                     radius of
                                                                     particle
                                                                        |
                                                                        V
RCoordinate                 |---------------------------------------|-------|
                          r0=0                                      a      end
                                                                           of
                                                                          shell


                    radii of particles                               radius of
                  (or particle+surface)                              particle
                            |                                           |
                            V                                           v
RThetaCoordinates   |--------------|-------------|------------------|-------|
(radial part)       0            sigma           r0                 a      end
                                                                           of
                                                                          shell
'''


class Coordinate(object):
    """1 coordinate (2 in case of RThetaCoordinates) of a position vector.

    The position of a particle inside a shell can be specified by a vector. 
    Which coordinates are used for this vector depends on the surface the 
    particle is on (if any). For example for a particle inside a sphere, not 
    on a surface, this is (r, theta, phi). But for a particle on a 
    CylindricalSurface, moving in 1 dimension, just (z) is sufficient.

    In the files single.py and pair.py it is defined which coordinates are 
    needed to define the position of the particle inside its shell for each 
    Single and Pair. Then, for each such coordinate, that can not be randomly 
    chosen, a domain is constructed.

    There are 3 types of domains:
        * ZCoordinate: cartesian.
        * RCoordinate: radially symetric.
        * RThetaCoordinates: cylindrical/spherical.

    All Greens Functions should be called from this file only.

    """
    def __init__(self, gf, a=None):
        self.gf = gf         # Green's function.
        if not a == None:
            self.a = a       # Outer radius.

    def geta(self):
        return self._a
    def seta(self, a):
        self._a = a
        # Set a of Greens' function here.
        self.gf.seta(a)
    a = property(geta, seta)


class RCoordinate(Coordinate):
    """The initial position of the particle is always at r = 0, so theta can and 
    should be choosen at random.

    """
    def __init__(self, gf, a):
        Coordinate.__init__(self, gf, a)

    def drawTime(self):
        try:
            rnd = myrandom.uniform()
            log.debug('        *Radial drawTime. ') #+ str(self.gf))
            dt = self.gf.drawTime(rnd)
        except Exception, e:
            raise Exception('gf.drawTime() failed, %s, rnd = %g, a = %g, %s' %
                            (str(e), rnd, self.a, gf.dump()))
        return dt

    def drawEventType(self, dt):
        # Either SINGLE_ESCAPE or COM_ESCAPE, doesn't matter.
        return EventType.COM_ESCAPE

    def drawDisplacement(self, dt, eventType):
        if(eventType == EventType.COM_ESCAPE or
           eventType == EventType.SINGLE_ESCAPE):
            # Escape through this coordinate. We already know the new r.
            return self.a

        try:
            rnd = myrandom.uniform()
            log.debug('        *Radial drawR. ') #+ str(self.gf))
            r = self.gf.drawR(rnd, dt)
            while r > self.a: # redraw; shouldn't happen often
                if __debug__:
                    log.debug('        *Radial drawR: redraw')
                rnd = myrandom.uniform()
                r = self.gf.drawR(rnd, dt)
        except Exception, e:
            raise Exception('gf.drawR failed, %s, rnd = %g, dt = %g, a = %g, '
                            '%s' % (str(e), rnd, dt, self.a, gf.dump()))

        return r


class RThetaCoordinates(Coordinate):
    """Used for PairGreensFunctions (3D as well as 2D).

    """
    def __init__(self, gf, sigma, r0, a):
        Coordinate.__init__(self, gf, a)    # Greens function holds D, sigma, k.
        self.sigma = sigma                  # Inner radius.
        self.r0 = r0                        # Starting position.

    def drawTime(self):
        try:
            rnd = myrandom.uniform()
            log.debug('        *Radial2D drawTime. ') #+ str(self.gf))
            dt = self.gf.drawTime(rnd, self.r0)
        except Exception, e:
            raise Exception('gf.drawTime() failed, %s, rnd = %g, sigma = %g, '
                            'r0 = %g, a = %g, %s' %
                            (str(e), rnd, self.sigma, self.r0, self.a,
                             gf.dump()))
        return dt

    def drawEventType(self, dt):
        try:
            rnd = myrandom.uniform()
            log.debug('        *Radial2D drawEventType. ') #+ str(self.gf))
            eventType = self.gf.drawEventType(rnd, self.r0, dt)
        except Exception, e:
            raise Exception('gf.drawEventType() failed, %s, sigma = %g,'
                            'r0 = %g, a = %g, dt = %g, %s' %
                            (str(e), self.sigma, self.r0, self.a, dt,
                             gf.dump()))
        return eventType     # (PAIR_REACTION or IV_ESCAPE)

    def drawDisplacement(self, gf, dt, eventType):
        if eventType == EventType.PAIR_REACTION:
            r = self.sigma
        elif eventType == EventType.IV_ESCAPE:
            r = self.a
        else:
            r = self.drawR_pair(gf, dt)
        theta = self.drawTheta_pair(gf, r, dt)
        return r, theta

    def drawR_pair(self, gf, dt):
        """Draw r for the pair inter-particle vector.

        """
        try:
            rnd = myrandom.uniform()
            log.debug('        *Radial2D drawR_pair. ') #+ str(self.gf))
            r = gf.drawR(rnd, self.r0, dt)
            # redraw; shouldn't happen often
            while r >= self.a or r <= self.sigma: 
                if __debug__:
                    log.info('    drawR_pair: redraw')
                #self.sim.rejectedMoves += 1  #FIXME:
                rnd = myrandom.uniform()
                r = gf.drawR(rnd, self.r0, dt)
        except Exception, e:
            raise Exception('gf.drawR_pair() failed, %s, rnd = %g, sigma = %g, '
                            'r0 = %g, a = %g, dt = %g, %s' %
                            (str(e), rnd, self.sigma, self.r0, self.a, dt,
                             gf.dump()))
        return r

    def drawTheta_pair(self, gf, r, dt):
        """Draw theta for the pair inter-particle vector.

        """
        try:
            rnd = myrandom.uniform()
            log.debug('        *Radial2D drawTheta_pair. ')#+ str(self.gf))
            theta = gf.drawTheta(rnd, r, self.r0, dt)
        except Exception, e:
            raise Exception('gf.drawTheta() failed, %s, rnd = %g, r = %g, '
                            'sigma = %g, r0 = %g, a = %g, dt = %g' %
                            (str(e), rnd, r, self.sigma, self.r0, 
                             self.a, dt))#, gf.dump()))

        # Heads up. For cylinders theta should be between [-pi, pi]. For 
        # spheres it doesn't matter.
        return myrandom.choice(-1, 1) * theta

