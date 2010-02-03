import math
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


