#!/usr/bin/env python

import math
import numpy



Pi = numpy.pi

'''
g(r,t) = [ 1 / (2 a2 r ) ] Sum_j j sin( Pi j r )
              exp( - D ( ( Pi j )2 t ) / a2 )
'''

def fpp1(r,t,D,a):
    a1 = ( 1.0 / (2.0 * a*a * r ) )

    a2 = 0.0
    j = 1
    while True:
        term = j * math.sin( Pi * j * r ) *\
               math.exp( - ( D * t * (Pi * j)**2 ) / ( a * a ) )
        a2 += term
        #print term
        if abs( a2 * 1e-18 ) > abs( term ) or term == 0:
            break

        j += 1

    print a1,a2


    return a1 * a2 * 4.0 * Pi * r * r


from scipy import integrate

print 'fpp', fpp1( 1,1e-2,1,1.0 )

print integrate.quad( fpp1, 0.0, 1.2, args=(1e-1,1.0,1.2) )
