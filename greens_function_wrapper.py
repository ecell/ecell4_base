import numpy
import myrandom
from _greens_functions import EventType
from utils import *

import logging
log = logging.getLogger('ecell')


def drawTime_wrapper(gf):
    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        *drawTime. ' + gf.__class__.__name__)
    try:
        # Todo. Let gf handle this.
        dt = gf.drawTime(rnd)
    except Exception, e:
        raise Exception('gf.drawTime() failed, '
                        '%s, rnd = %g, %s' %
                        (str(e), rnd, gf.dump()))
    return dt

def draw_eventtype_wrapper(gf, dt):
    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        *drawEventType. ' + gf.__class__.__name__)
    try:
        event_type = gf.drawEventType(rnd, dt)
    except Exception, e:
        raise Exception('gf.drawEventType() failed, '
                        '%s, rnd = %g, dt = %g, %s' %
                        (str(e), rnd, dt, gf.dump()))
    return event_type

# Todo. Returns r, not displacement.
def draw_displacement_wrapper(gf, dt, event_type, a, sigma=None):
    if(((event_type == EventType.COM_ESCAPE or
         event_type == EventType.SINGLE_ESCAPE) and sigma == None) or
       (event_type == EventType.IV_ESCAPE and sigma != None)):
        # Escape through this coordinate. We already know the new r.
        # Todo. Let gf handle this.
        return a
    elif event_type == EventType.IV_REACTION and sigma != None:
        # Todo. Let gf handle this.
        return sigma

    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        *drawR. ' + gf.__class__.__name__)
    try:
        r = gf.drawR(rnd, dt)
        while r > a or r <= sigma: # redraw; shouldn't happen often
            if __debug__:
                log.debug('        *drawR: redraw')
            rnd = myrandom.uniform()
            r = gf.drawR(rnd, dt)
    except Exception, e:
        raise Exception('gf.drawR() failed, '
                        '%s, rnd = %g, dt = %g, %s' %
                        (str(e), rnd, dt, gf.dump()))

    return r

# Todo. Returns (r,theta), not displacement.
def draw_displacement_iv_wrapper(gf, dt, event_type, a, sigma):
    def drawTheta_wrapper(gf, r, dt):
        """Draw theta for the inter-particle vector.

        """
        rnd = myrandom.uniform()

        if __debug__:
            log.debug('        *drawTheta. ' + gf.__class__.__name__)
        try:
            theta = gf.drawTheta(rnd, r, dt)
        except Exception, e:
            raise Exception('gf.drawTheta() failed, '
                            '%s, rnd = %g, r = %g, dt = %g' %
                            (str(e), rnd, r, dt))#, gf.dump()))

        # Heads up. For cylinders theta should be between [-pi, pi]. For 
        # spheres it doesn't matter.
        return myrandom.choice(-1, 1) * theta

    r = draw_displacement_wrapper(gf, dt, event_type, a, sigma)
    theta = drawTheta_wrapper(gf, r, dt)
    return r, theta

