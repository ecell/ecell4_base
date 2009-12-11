
from utils import *
import math
import numpy
import myrandom as rng

class Surface:
    def __init__( self ):
        pass

    '''
    This method computes the absolute distance between pos and
    a closest point on the surface.
    '''
    def distance( self, pos ):
        return abs( self.signedDistance( pos ) )

    '''
    Same as distance(), but the return value is positive
    or negative depending on in which side the position is.
    
    When this surface defines a closed region in space,
    negative value means that the position is inside.
    '''
    def signedDistance( self, pos ):
        pass

    '''
    Returns a random position equidistributed within
    the region in space defined by negative signed distance.
    See also signedDistance().
    '''
    def randomPosition( self ):
        pass


'''
Planar boundary given by a x + b y + c z + d = 0

'''
class PlanarSurface( Surface ):

    '''
    params: [ a, b, c, d ] as in
    a x + b y + c z + d = 0
    where ( a, b, c ) is a normal vector of the plane.

    '''
    def __init__( self, params ):

        self.setParams( params )

    def signedDistance( self, pos ):
        return numpy.dot( self.hessianNormal,\
                          numpy.array( [ pos[0], pos[1], pos[2], 1.0 ] ) )

    def setParams( self, params ):
        self.params = numpy.array( params )

        scaling = length( self.params[:3] )

        # Hessian normal form [ nx, ny, nz, p ].
        self.hessianNormal = self.params / scaling


    def getParams( self ):
        return self.params

    def randomPosition( self ):
        raise 'not supported.  specified space has infinite volume.'
    

'''

'''
class CuboidalSurface( Surface ):

    '''
    origin -- = [ x0, y0, z0 ] is the origin.
    size -- = [ sx, sy, sz ] is the vector from the origin to the diagonal
              point.
    '''

    def __init__( self, origin, size ):
        self.setParams( origin, size )

    def signedDistance( self, pos ):
        dists = numpy.concatenate( ( self.origin - pos,
                                     self.origin+self.size-pos ) )
        absdists = numpy.abs( dists )
        i = numpy.argmin( absdists )

        return dists[i]

    def randomPosition( self ):
        return numpy.array( [ rng.uniform( self.origin[0], self.size[0] ),
                              rng.uniform( self.origin[1], self.size[1] ),
                              rng.uniform( self.origin[2], self.size[2] ) ]
                            )

    
    def setParams( self, origin, size ):

        self.origin = numpy.array( origin )
        self.size = numpy.array( size )
        self.midpoint = ( self.size - self.origin ) * 0.5

    def getParams( self ):
        return self.params


class SphericalSurface( Surface ):

    '''
    origin -- = [ x0, y0, z0 ] is the origin and
    radius -- the radius of the sphere.

    ( x - x0 )^2 + ( y - y0 )^2 + ( z - z0 )^2 = r^2
    '''
    def __init__( self, origin, radius ):

        self.setParams( origin, radius )

    def signedDistance( self, pos ):
        
        return math.sqrt( ( ( pos - self.origin ) ** 2 ).sum() ) - self.radius

    def randomPosition( self ):
        pos = randomUnitVectorS()
        pos[0] *= self.radius
        return sphericalToCartesian( pos )

    def setParams( self, origin, radius ):
        self.origin = numpy.array( origin )
        self.radius = radius

    def getParams( self ):
        return self.params
