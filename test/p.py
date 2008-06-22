#!/usr/bin/env python


import sys

import numpy
from matplotlib.pylab import *


import _gfrd

N_A = 6.0221367e23


sigma = 1e-8
#r0 = sigma 
D = 1e-12
#kf = 1000 * sigma * D
#kf=1e-8
#kf=1e-10
kf=1e-10
#a = 1e-7
a = sigma*5
#r0 = a * (1.0-1e-7)
r0 = sigma * 3
#r0 = a * 0.999
#r0 = (a-sigma) * 0.5 + sigma

tau = sigma*sigma / D
#T = tau * .1
#T = 1e-300
T = 1e-2

rmin = sigma


def plot_p_survival( gf, T ):

    N = 1000

    x = numpy.mgrid[0:T:T/N]
    parray1 = numpy.array( [ gf.p_survival( t ) for t in x ] )
    plot( x, parray1, '-', label='psurvival' )




def plot_p_survival_i( gf ):

    N = 1000

    x = range( N )
    parray1 = numpy.array( [ gf.p_survival_i_exp( i, T, r0 ) for i in x ] )
    print len(parray1[:-1]), len( parray1[1:])
    parray2 = parray1[:-1:2] + parray1[1::2]
    parray2 = parray2[:-1:2] + parray2[1::2]
    plot( range(len(parray2)), parray2, '-', label='psurvival_i' )
    plot( range(len(parray2)), parray2.cumsum(), '-', label='psurvival_i_sum' )
    #plot( range(N), parray1, '.', label='psurvival_i' )

def plot_p_survival_alpha( gf ):

    N = 1000

    alpha = numpy.array(range( N )) * 1e-3
    parray1 = numpy.array( [ gf.p_survival_i_alpha( al, T, r0 ) for al in alpha ] )
    print len(parray1[:-1]), len( parray1[1:])
    parray2 = parray1[:-1:2] + parray1[1::2]
    parray2 = parray2[:-1:2] + parray2[1::2]
    #plot( range(len(parray2)), parray2, '-', label='psurvival_i' )
    #plot( range(len(parray2)), parray2.cumsum(), '-', label='psurvival_i_sum' )
    plot( alpha, parray1, '.', label='psurvival_i' )



def plot_p_leaveas( gf, t ):

    N = 100000

    tmax = 1e-3
    tmin = 1e-10


    ttick = ( tmax - tmin ) / N
    tarray = numpy.mgrid[tmin:tmax:ttick]

    parray1 = array( [ 1 - gf.p_survival( t, r0 ) for t in tarray ] )
    semilogx( tarray , parray1, '-', label='psurvival' )

    gf2 = _gfrd.BasicPairGreensFunction( D, kf, sigma )
    parray2 = array( [ 1 - gf2.p_survival( t, r0 ) for t in tarray ] )
    semilogx( tarray , parray2, '-', label='psurvival basic' )

#     parray2 = array( [ gf.p_leavea( t, r0 )  for t in tarray ] )
#     parray2 = 1 - parray2# / gf.p_leavea( 0, r0 )
#     semilogx( tarray , parray2, '-', label='pleavea' )

#     parray3 = array( [ gf.p_leaves( t, r0 )  for t in tarray ] )
#     parray3 = 1 - parray3# / gf.p_leaves( 0, r0 )
#     semilogx( tarray , parray3, '-', label='pleaves' )

#     semilogx( tarray , parray2 + parray3 - 1, '-', label='s+a' )

    #semilogx( tarray , parray2 + parray3, '-', label='a+s' )


def plot_leaveas( gf, t ):

    N = 3000

    #tmax = 2.4e-5
    #tmin = 1.1e-5
    tmax = 2.5e-2
    tmin = 2.2e-8

    ttick = ( tmax - tmin ) / N
    tarray = numpy.mgrid[tmin:tmax:ttick]

    #parray1 = array( [ 1 - gf.p_survival( t, r0 ) for t in tarray ] )
    #semilogx( tarray , parray1, '-', label='psurvival' )

    parray2 = array( [ gf.leavea( t, r0 ) * 4 * numpy.pi * a * a
                       for t in tarray ] )
    parray3 = array( [ gf.leaves( t, r0 ) * 4 * numpy.pi * sigma * sigma
                       for t in tarray ] )
    parray4 = array( [ gf.dp_survival( t, r0 )  for t in tarray ] )

    #semilogx( tarray, parray2 / (parray2+parray3), '-', label='leavea' )
    semilogx( tarray, parray2, '-', label='leavea' )
    semilogx( tarray, parray3, '-', label='leaves' )
    #semilogx( tarray, parray3 / (parray2+parray3), '-', label='leaves' )
    #semilogx( tarray, parray4 / gf.dp_survival(0,r0) , '-', label='dp_survival' )
    #semilogx( tarray, (parray2 + parray3)/(parray2[0]+parray3[0]) , '-', label='a+s' )

    #semilogx( tarray , parray2, '-', label='a' )
    #semilogx( tarray , parray3, '-', label='s' )




def plot_p_int_r( gf, t ):

    N = 10000

    rmax = min( a, r0 + 4 * math.sqrt( 6 * D * t ) )
    #rmax = a
    #rmin = max( sigma, r0 - 4 * math.sqrt( 6 * D * t ) )
    #rmin = max( 0, r0 - 4 * math.sqrt( 6 * D * t ) )
    rmin = 0

    rtick = ( rmax - rmin ) / N
    rarray = numpy.mgrid[rmin:rmax:rtick]

    #surv = gf.p_survival( t, r0 )
    surv = gf.p_survival( t )
    print surv
    #parray = array( [ gf.p_int_r( r, t, r0 ) for r in rarray ] ) / surv
    parray = array( [ gf.p_int_r( r, t ) for r in rarray ] ) / surv

    plot( rarray / sigma , parray, '-', label='f' )


def plot_ip_theta( gf, r, t ):

    N = 300

    thetamax = 0
    thetamin = numpy.pi

    thetatick = ( thetamax - thetamin ) / N
    thetaarray = numpy.mgrid[thetamin:thetamax:thetatick]

    p0 = gf.p_0( t, r, r0 ) * 2
    parray = array( [ gf.ip_theta( theta, r, r0, t ) 
                      for theta in thetaarray ] ) 
    parray /= p0

    plot( thetaarray, parray, '-', label='f' )

def p( r, t ):
    surv = gf.p_survival( t, r0 )
    return gf.p_int_r( r, t, r0 ) / surv

if __name__ == '__main__':

    #gf = _gfrd.FirstPassagePairGreensFunction( D, kf, sigma )
    #gf = _gfrd.FirstPassageNoCollisionPairGreensFunction( D )
    gf = _gfrd.FirstPassageGreensFunction( D )
    gf.seta( a )
     
    #plot_p_int_r( gf, T )
    plot_p_int_r( gf, 1e-6 )
    #plot_p_survival( gf, T )

    #plot_ip_theta( gf, r0, T )
    #plot_p_leaveas( gf, r0 )
    #plot_leaveas( gf, r0 )
    #plot_p_survival_i( gf )
    #plot_p_survival_alpha( gf )

    #xlabel( 'r / sigma' )
    #ylabel( 'p_irr' )
    legend()
    show()

