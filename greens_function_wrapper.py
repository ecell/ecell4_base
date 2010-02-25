import numpy
import myrandom
from _gfrd import EventType
from utils import *

import logging
log = logging.getLogger('ecell')


def draw_time_wrapper(gf, r0=None):
    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        *drawTime. ') #+ str(gf))
    try:
        if r0 == None:
            # Todo. Let gf handle this.
            dt = gf.drawTime(rnd)
        else:
            dt = gf.drawTime(rnd, r0)
    except Exception, e:
        raise Exception('gf.drawTime() failed, '
                        '%s, rnd = %g, %s' %
                        (str(e), rnd, gf.dump()))
    return dt

def draw_eventtype_wrapper(gf, dt, r0):
    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        *drawEventType. ') #+ str(gf))
    try:
        eventType = gf.drawEventType(rnd, r0, dt)
    except Exception, e:
        raise Exception('gf.drawEventType() failed, '
                        '%s, rnd = %g, r0 = %g, dt = %g, %s' %
                        (str(e), rnd, r0, dt, gf.dump()))
    return eventType

def draw_displacement_wrapper(gf, dt, eventType, a):
    if(eventType == EventType.COM_ESCAPE or
       eventType == EventType.SINGLE_ESCAPE):
        # Escape through this coordinate. We already know the new r.
        # Todo. Let gf handle this.
        return a

    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        *drawR. ') #+ str(gf))
    try:
        r = gf.drawR(rnd, dt)
        while r > a: # redraw; shouldn't happen often
            if __debug__:
                log.debug('        *drawR: redraw')
            rnd = myrandom.uniform()
            r = gf.drawR(rnd, dt)
    except Exception, e:
        raise Exception('gf.drawR() failed, '
                        '%s, rnd = %g, dt = %g, %s' %
                        (str(e), rnd, dt, gf.dump()))

    return r

def draw_displacement_iv_wrapper(gf, r0, dt, eventType, a, sigma):
    def draw_r_wrapper(gf, r0, dt, a, sigma):
        """Draw r for the inter-particle vector.

        """
        rnd = myrandom.uniform()

        if __debug__:
            log.debug('        *drawR_pair. ') #+ str(gf))
        try:
            r = gf.drawR(rnd, r0, dt)
            # redraw; shouldn't happen often
            while r >= a or r <= sigma: 
                if __debug__:
                    log.info('    drawR_pair: redraw')
                #self.sim.rejectedMoves += 1  #FIXME:
                rnd = myrandom.uniform()
                r = gf.drawR(rnd, r0, dt)
        except Exception, e:
            raise Exception('gf.drawR_pair() failed, '
                            '%s, rnd = %g, r0 = %g, dt = %g, %s' %
                            (str(e), rnd, r0, dt, gf.dump()))
        return r

    def draw_theta_wrapper(gf, r, r0, dt):
        """Draw theta for the inter-particle vector.

        """
        rnd = myrandom.uniform()

        if __debug__:
            log.debug('        *drawTheta. ')#+ str(gf))
        try:
            theta = gf.drawTheta(rnd, r, r0, dt)
        except Exception, e:
            raise Exception('gf.drawTheta() failed, '
                            '%s, rnd = %g, r = %g, r0 = %g, dt = %g' %
                            (str(e), rnd, r, r0, dt))#, gf.dump()))

        # Heads up. For cylinders theta should be between [-pi, pi]. For 
        # spheres it doesn't matter.
        return myrandom.choice(-1, 1) * theta

    if eventType == EventType.IV_REACTION:
        # Todo. Let gf handle this.
        r = sigma
    elif eventType == EventType.IV_ESCAPE:
        # Todo. Let gf handle this.
        r = a
    else:
        r = draw_r_wrapper(gf, r0, dt, a, sigma)
    theta = draw_theta_wrapper(gf, r, r0, dt)
    return r, theta

