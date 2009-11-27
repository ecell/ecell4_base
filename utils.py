
import math
import numpy
from numpy.random import uniform
import scipy

import _gfrd

Pi = scipy.pi
Pi2 = scipy.pi * 2.0
PiSqrt = math.sqrt( scipy.pi )

N_A = 6.0221367e23
INF = numpy.inf

ZEROPOS = numpy.array( [ 0., 0., 0. ] )
NOWHERE = numpy.array( ( INF, INF, INF ) )

def Mtom3( rate ):
    return rate / ( 1000 * N_A )

def meanArrivalTime( r, D ):
    return ( r * r ) / ( 6.0 * D )

def uniq( l ):
    nset = {}
    map( nset.__setitem__, l, [] )
    return nset.keys()

def _cyclic_transpose(pos1, pos2, world_size):
    '''
    Transpose the position pos1 so that it can be used with another 
    position pos2.

    pos1 is transposed into one of mirror images of the cyclic boundary
    condition so that the distance between pos1 and pos2 is smallest.

    Both of given pos1 and pos2 must be within the cyclic boundary.  However,
    note that the returned transposed pos1 may not be within the cyclic boundary.
    '''
    return _gfrd.cyclic_transpose(pos1, pos2, world_size)

def cyclic_transpose( pos1, pos2, fsize ):
    halfsize = fsize * 0.5

    diff = pos2 - pos1
  
    reloc = numpy.greater( diff, halfsize ) * fsize - \
        numpy.less( diff, - halfsize ) * fsize

    return pos1 + reloc

def distanceSq_Simple( position1, position2, fsize = None ):
    return _gfrd.distanceSq( position1, position2 )

def distance_Simple( position1, position2, fsize = 0 ):
    return _gfrd.distance( position1, position2 )

def distanceSqArray_Simple( position1, positions, fsize = None ):
    return numpy.square( positions - position1 ).sum( 1 )

def distanceArray_Simple( position1, positions, fsize = None ):
    return numpy.sqrt( distanceSqArray_Simple( position1, positions ) )

def distanceSq_Cyclic( position1, position2, fsize ):
    return _gfrd.distanceSq_Cyclic( position1, position2, fsize )

def distance_Cyclic( position1, position2, fsize ):
    return _gfrd.distance_Cyclic( position1, position2, fsize )

def distanceSqArray_Cyclic( position1, positions, fsize ):
    diff = numpy.abs( positions - position1 )
    diff -= numpy.greater( diff, fsize * 0.5 ) * fsize # transpose
    return numpy.square( diff ).sum( 1 )

def distanceArray_Cyclic( position1, positions, fsize = 0 ):
    return numpy.sqrt( distanceSqArray_Cyclic( position1, positions, fsize ) )

def cartesianToSpherical( c ):
    # x, y, z = c
    r = length( c )
    theta = math.acos( c[2] / r )
    phi = math.atan2( c[1], c[0] )
    if phi < 0.0:  # atan2 returns [- PI, PI]
        phi += 2.0 * Pi
    return numpy.array( [ r, theta, phi ] )

def sphericalToCartesian( s ):
    #FIXME: it's possible that the below is a source of some bias.
    r, theta, phi = s
    sintheta = math.sin( theta )
    return numpy.array( [ r * math.cos( phi ) * sintheta,
                          r * math.sin( phi ) * sintheta,
                          r * math.cos( theta ) ] )

def randomUnitVectorS():
    s = numpy.array( [ 1.0, uniform( 0, Pi ),
                       uniform( 0, Pi2 ) ] )
    return s

def randomUnitVector():
    v = uniform( -1, 1, 3 )
    return _gfrd.normalize(v, 1)

def randomVector( r ):
    v = uniform( -1, 1, 3 )
    return _gfrd.normalize(v, r)

def length( a ):
    return _gfrd.length( a )

def normalize( a, l=1 ):
    return _gfrd.normalize( a, l )

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

def rotateVector( v, r, alpha ):
    '''
    v: vector to rotate
    r: normalized rotation axis
    alpha: rotation angle in radian
    '''
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

def _calculate_pair_CoM( pos1, pos2, D1, D2, world_size ):
    return _gfrd.calculate_pair_CoM(pos1, pos2, D1, D2, world_size);

def calculate_pair_CoM( pos1, pos2, D1, D2, world_size ):
    '''
    Just a free func ver of Pair.getCoM().
    '''
    pos2t = cyclic_transpose( pos2, pos1, world_size )
    retval = ( ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 ) ) % world_size
    #cpp_result = _calculate_pair_CoM( pos1, pos2, D1, D2, world_size )
    #assert retval == cpp_result, "%g:%g" % ( retval, cpp_result )
    return retval

def _apply_boundary(pos, world_size):
    return _gfrd.apply_boundary(pos, world_size)

def apply_boundary(pos, world_size):
    return pos % world_size

def permutate(seq):
    """permutate a sequence and return a list of the permutations"""
    if not seq:
        return [seq] # is an empty sequence
    else:
        temp = []

        for k in range(len(seq)):
            part = seq[:k] + seq[k+1:]
            for m in permutate(part):
                temp.append(seq[k:k+1] + m)
        return temp

def k_D( D, sigma ):
    return 4.0 * numpy.pi * D * sigma

def k_a( kon, kD ):
    if kon > kD:
        raise RuntimeError, 'kon > kD.'
    ka = 1 / ( ( 1 / kon ) - ( 1 / kD ) )
    return ka

def k_d( koff, kon, kD ):
    return k_a( kon, kD ) * koff / kon

def k_on( ka, kD ):
    kon = 1 / ( ( 1 / kD ) + ( 1 / ka ) )  # m^3/s
    return kon

def C2N( c, V ):
    return c * V * N_A  # round() here?

