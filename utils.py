
import math
import numpy
import scipy


Pi = scipy.pi
Pi2 = scipy.pi * 2.0
PiSqrt = math.sqrt( scipy.pi )

N_A = 6.0221367e23
INF = numpy.Inf

NOWHERE = numpy.array( ( INF, INF, INF ) )


def MsTom3s( rate ):
    return rate / ( 1000 * N_A )


'''
Transpose the position pos so that it can be used with another 
position vector basis.


Both pos and basis must be within the cyclic boundary.
'''

def cyclicTranspose( pos, basis, fsize ):
    halfsize = fsize * 0.5

    diff = basis - pos

    reloc = numpy.greater( diff, halfsize ) * fsize - \
        numpy.less( diff, - halfsize ) * fsize

    return pos + reloc


def distanceSq_Simple( position1, position2, fsize = 0 ):
    diff = position1 - position2
    return numpy.dot( diff, diff )

def distance( position1, position2, fsize = 0 ):
    return math.sqrt( distanceSq_Simple( position1, position2 ) )

def distanceSqArray_Simple( position1, positions, fsize = 0 ):
    tmp = positions - position1
    return ( tmp * tmp ).sum(1)


def distanceSq_Cyclic( position1, position2, fsize ):

    diff = numpy.abs( position2 - position1 )
    diff -= numpy.greater( diff, fsize * 0.5 ) * fsize # transpose

    return numpy.dot( diff, diff )


def distanceSqArray_Cyclic( position1, positions, fsize ):

    diff = numpy.abs( positions - position1 )
    diff -= numpy.greater( diff, fsize * 0.5 ) * fsize # transpose

    return ( diff * diff ).sum(1)



def cartesianToSpherical( c ):
    # x, y, z = c
    r = math.sqrt( ( c ** 2 ).sum() )
    theta = math.acos( c[2] / r )
    phi = math.atan2( c[1], c[0] )
    if phi < 0.0:  # atan2 returns [- PI, PI]
        phi += 2.0 * Pi
    return numpy.array( [ r, theta, phi ] )


def sphericalToCartesian( s ):
    r, theta, phi = s
    sintheta = math.sin( theta )
    return numpy.array( [ r * math.cos( phi ) * sintheta,
                          r * math.sin( phi ) * sintheta,
                          r * math.cos( theta ) ] )


def randomUnitVectorS():
    s = numpy.array( [ 1.0, numpy.random.uniform( 0, Pi ),
                       numpy.random.uniform( 0, Pi2 ) ] )
    return s


def randomUnitVector():
    return sphericalToCartesian( randomUnitVectorS() )


def length( a ):
    return math.sqrt( numpy.dot( a, a ) )

def normalize( a ):
    return a / length( a )


def vectorAngle( a, b ):
    cosangle = numpy.dot( a, b ) / ( length( a ) * length( b ) )
    return math.acos( cosangle )

def vectorAngleAgainstZAxis( b ):
    cosangle = b[2] / length( b )
    return math.acos( cosangle )

def crossproduct( a, b ):
    M = numpy.array( [ [    0.0, - a[2],   a[1] ],
                       [   a[2],    0.0, - a[0] ],
                       [ - a[1],   a[0],    0.0 ] ] )
    return numpy.dot( M, b )

def crossproductAgainstZAxis( a ):
    return numpy.array( [ - a[1], a[0], 0.0 ] )


'''
v: vector to rotate
r: normalized rotation axis
alpha: rotation angle in radian
'''
def rotateVector( v, r, alpha ):
    cosalpha = math.cos( alpha )
    sinalpha = math.sin( alpha )
    cosalphac = 1.0 - cosalpha

    M = numpy.array( [ [ cosalpha + cosalphac * r[0] * r[0],
                         cosalphac * r[0] * r[1] - r[2] * sinalpha,
                         cosalphac * r[0] * r[2] + r[1] * sinalpha ],
                       [ cosalphac * r[0] * r[1] + r[2] * sinalpha,
                         cosalpha + cosalphac * r[1] * r[1],
                         cosalphac * r[1] * r[2] - r[0] * sinalpha ],
                       [ cosalphac * r[0] * r[2] - r[1] * sinalpha,
                         cosalphac * r[1] * r[2] + r[0] * sinalpha,
                         cosalpha + cosalphac * r[2] * r[2] ] ] )

    return numpy.dot( M,v )
