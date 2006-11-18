#!/usr/bin/env python


import math
import random
import numpy

import scipy
from scipy import optimize
from scipy import integrate
from scipy import special

N_A = 6.0221367e23


def p1( D, dt ):
    ro = math.sqrt( 2.0 * D * dt )

    displacement = numpy.random.normal( 0.0, ro, 3 )

    return displacement


def p2_R( D1, D2, dt ):
    ro = math.sqrt( 2.0 * ( D1 + D2 ) * dt )

    displacement = numpy.random.normal( 0.0, ro, 3 )

    return displacement


def erfc_scaled( x ):
    '''
    Calculates exp( x^2 ) * erfc( x )

    See asymptotic expansion here:
    http://en.wikipedia.org/wiki/Error_function
    
    '''

    if x < 25.0:  
        return numpy.exp( x * x ) * special.erfc( x )
    else: 
        M_1_SQRTPI = 1.0 / math.sqrt( scipy.pi )

        x2sq = 4.0 * x * x

        # up to second term in the expansion.
        # abs err ~= 9e-8 at x == 20, 3e-8 at x == 25

        # the third term below and beyond doesn't have a big
        # contribution for large x.
        # - ( 2 / ( x2sq * x2sq * x2sq ) )       

        return M_1_SQRTPI * ( 1.0 / x ) * \
               ( 1.0 + ( - 2.0 / x2sq ) + ( 4.0 / ( x2sq * x2sq ) ) )

    


def W( a, b ):
    '''
    exp( 2 a b + b^2 ) erfc( a + b )
    '''
    
    #print 'W a b', a, b

    return numpy.exp( - a * a ) * erfc_scaled( a + b )


def p_irr( r, t, radius, r0, D, k ):

    '''
    dr = r - radius
    alpha = ( 1 + (kr/kD)) sqrt(D) / radius
    kD = 4 pi radius D

    (
    1 / ( sqrt( 4 pi t ) ) *

    ( exp( - ( dr^2 ) / 4Dt ) + exp( - (dr - 2 radius)^2 / 4Dt ) )

    - alpha W( ( dr - 2radius)/ sqrt(4Dt), alpha sqrt(t) )
    )

    / 4 pi r r0 sqrt(D)

    '''

    sqrtD = numpy.sqrt(D)
    Dt4 = 4.0 * D * t 


    dr = r - r0
    kD = 4.0 * scipy.pi * radius * D
    alpha = ( ( 1.0 + ( k / kD ) ) * sqrtD ) / radius



    a1 = math.exp( - ( dr * dr ) / Dt4 ) +\
         math.exp( - ( ( r + r0 - 2.0 * radius ) ** 2 ) / Dt4 )

    #print math.exp( - ( dr * dr ) / Dt4 ),\
    #     math.exp( - ( ( r + r0 - 2.0 * radius ) ** 2 ) / Dt4 )

    a1 /= numpy.sqrt( 4.0 * scipy.pi * t )
    #print a1
    a2 = alpha * W( ( r + r0 - 2.0 * radius ) / numpy.sqrt( Dt4 ),
                    alpha * numpy.sqrt(t) )

    #print a2
    #a3 = math.sqrt( 4.0 * scipy.pi * t )
    #print a1, a2, a3

    factor = 4.0 * scipy.pi * r * r0 * sqrtD

    #print factor

    
    return  (a1/factor - a2/factor) * r * r

    

def p_survival2( t, radius, r0, D, k ):

    kD = 4.0 * scipy.pi * radius * D
    alpha = ( ( 1.0 + k / kD ) * math.sqrt( D ) ) / radius

    #print 'kD alpha', kD, alpha

    r0__radius = ( r0 - radius )
    if r0__radius <= 0.0:
        r0__radius = 0.0

    sqrt__4_D_t = math.sqrt( 4.0 * D * t )

    a1 = ( radius / r0 )
    a2 = ( k / ( k + kD ) )
    #    print 'k kD', k, kD
    
    a3 = special.erfc( r0__radius / sqrt__4_D_t )

    a4 = W( ( r0__radius / sqrt__4_D_t ), alpha * math.sqrt( t ) )

    #    print 'kd alpha', kD, alpha

    #print a1, a2, a3, a4

    return 1.0 - ( a1 * a2 * ( a3 - a4 ) )

def p2( r, t, radius, r0, D, k ):
    '''
  p := ( ( 1 / Sqrt[ 4 Pi t ] ) ( Exp[ - ( ( r - r0 ) ^ 2 / ( 4 D t ) ) ] + Exp[ - ( r + r0 - 2 radius ) ^ 2 / ( 4 D t ) ] ) - ( ( 1 + k / kD ) Sqrt[D] / radius ) Exp[ 2 ap bp + bp ^ 2 ] Erfc[ ap + bp ] ) / ( 4 Pi r r0 Sqrt[D] )


num1 = ( Exp[ - ( ( r - r0 ) ^ 2 / ( 4 D t ) ) ] + Exp[ - ( r + r0 - 2 radius ) ^ 2 / ( 4 D t ) ] ) / ( Sqrt[ 4 Pi t ] ) 

num2 = - ( ( 1 + k / kD ) Sqrt[D] / radius ) Exp[ 2 ap bp + bp ^ 2 ] Erfc[ ap + bp ]

      den = ( 4 Pi r r0 Sqrt[D] )

    '''

    kD = 4.0 * scipy.pi * radius * D
    alpha = ( ( 1.0 + k / kD ) * math.sqrt( D ) ) / radius
    
    Dt = D * t

    tSqrt = math.sqrt( t )
    DSqrt = math.sqrt( D )
    PiSqrt = math.sqrt( scipy.pi )
    DtSqrt = math.sqrt( Dt )
    
    a = ( r + r0 ) / ( 2.0 * DtSqrt ) - ( 2.0 * radius )/ ( 2 * DtSqrt )
    b = alpha * tSqrt
    print 'a', alpha
    print - ( ( ( r - r0 ) / ( 2 * DtSqrt ) ) ** 2 )
    print ( math.sqrt( 4 * scipy.pi * t ) ) 
    num1 = ( math.exp( - ( ( ( r - r0 ) / ( 2 * DtSqrt ) ) ** 2 ) ) \
             + math.exp( - ( ( r + r0 - 2 * radius ) ** 2 ) / ( 4 * Dt ) ) ) \
             / ( math.sqrt( 4 * scipy.pi * t ) )

    num2 = alpha * W( a, b )

#    den = 4.0 * scipy.pi * r * r0 * DSqrt
#    den = DSqrt
    
    print num1, num2
#    res = ( num1 - num2 ) / den
    res = ( num1 - num2 )
    return 



def q2( t, radius, r0, D, k ):
    '''
(Sqrt[D]*k*(-(kD/(E^((r0 - radius)^2/(4*D*t))*Sqrt[Pi]*Sqrt[t])) +
   (Sqrt[D]*E^((D*(k + kD)^2*t)/(kD^2*radius^2) +
       (Sqrt[D]*(k + kD)*(r0 - radius)*Sqrt[t])/(kD*radius*Sqrt[D*t]))*(k + kD)*
     Erfc[(Sqrt[D]*(k + kD)*Sqrt[t])/(kD*radius) +
       (r0 - radius)/(2*Sqrt[D*t])])/radius))/(kD^2*r0)



C = Sqrt[D]*k

num1 = -(kD/(E^((r0 - radius)^2/(4*D*t))*Sqrt[Pi]*Sqrt[t]))


num2a = Sqrt[D]

num2b = E^((D*(k + kD)^2*t)/(kD^2*radius^2) +
            (Sqrt[D]*(k + kD)*(r0 - radius)*Sqrt[t])/(kD*radius*Sqrt[D*t]))

num2c = (k + kD) *
         Erfc[(Sqrt[D]*(k + kD)*Sqrt[t])/(kD*radius) +
             (r0 - radius)/(2*Sqrt[D*t])]

num2 = ( num2a*b*c ) / radius

den =  kD^2*r0

    '''


    kD = 4.0 * scipy.pi * radius * D

    Dt = D * t
    
    tSqrt = math.sqrt( t )
    DSqrt = math.sqrt( D )
    PiSqrt = math.sqrt( scipy.pi )
    DtSqrt = math.sqrt( Dt )

    kD_k = kD + k
    r0__radius = r0 - radius
    kDradius = kD * radius

#    print 'kd,k', kD, k

    C = DSqrt * k
    den = ( kD ** 2 ) * r0

#    print 'c den', C,den

    num1 = - ( kD / ( math.exp( ( r0__radius ** 2 ) / ( 4 * Dt ) ) * PiSqrt * tSqrt ) )

    num2a = DSqrt

#    print ( D * ( kD_k ** 2 ) * t )
#    print ( ( kD ** 2 ) * ( radius ** 2 ) )

#    print ( DSqrt * kD_k * r0__radius * tSqrt ) / \
#                      ( kDradius * DtSqrt )

    num2b = math.exp( ( D * ( kD_k ** 2 ) * t ) / ( ( kD ** 2 ) * ( radius ** 2 ) ) + \
                      ( DSqrt * kD_k * r0__radius * tSqrt ) / \
                      ( kDradius * DtSqrt ) )

    num2c = kD_k * \
            special.erfc( ( DSqrt* kD_k *tSqrt)/(kDradius) + r0__radius / ( 2 * DtSqrt ) )

    print num1
    print num2a, num2b, num2c
    print C
    print den

    num2 = ( num2a * num2b * num2c ) / radius 


    r = C * ( num1 + num2 ) / den

    return r


if __name__ == '__main__':

    def s( f, u, radius, r0, D, k ):
        print ( f, u, radius, r0, D, k )
        print u - p_survival2( f, radius, r0, D, k )
        return u - p_survival2( f, radius, r0, D, k )
    
    def t( f, u, radius, r0, D, k ):
        print f,u
        return u - f


    def t_react2( radius, r0, D, k ):
        dtMax = 0.001
        u = 0.8
        
        infp = p_survival2_inf( radius, r0, D, k )
        print 'sur', infp
        if u <= infp:
            return inf
        
        maxtp = p_survival2( dtMax, radius, r0, D, k )
        print maxtp
        if u <= maxtp:
            return inf
        
        res = optimize.brenth( s, 1e-300, 1, args=( u, radius, r0, D, k) )
        return res


    radius = 5e-8
    r0 = 5.0e-8 + 1e-12
    D = 1e-11
    k = 1e9/ (1e3 * N_A)
    dt = 1e-10
    
#    print p_irr( r0, 1e-20, radius, r0, D,  k / (1e3* N_A)   )

    print p_survival2( 1e-1, radius, r0, D, k   )

#    print p2( r, dt, radius, r0, D, k )

#    print integrate.quad( p_irr, r0, r0, args=( 1e-12, radius, r0, D, k ) )


    #print t_react2( radius, r0, D, k / N_A )
#print q2( dt, radius, r0, D, k / N_A )



#print p_survival2_inf( radius, r0, D, k / N_A )

#print p_survival2( 1e-12, 1e300, 1e-11, 1e-10, 1e5 )





