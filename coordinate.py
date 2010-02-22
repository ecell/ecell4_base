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

    def drawDisplacement(self, dt):
        try:
            rnd = myrandom.uniform()
            log.debug('        *Radial drawR. ') #+ str(self.gf))
            r = self.gf.drawR(rnd, dt)
        except Exception, e:
            raise Exception('gf.drawR failed, %s, rnd = %g, dt = %g, a = %g, '
                            '%s' % (str(e), rnd, dt, self.a, gf.dump()))

        return r


