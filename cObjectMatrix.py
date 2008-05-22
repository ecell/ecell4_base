#!/usr/bin/env python


import numpy

import object_matrix

from utils import *


class ObjectMatrix( object ):

    def __init__( self ):

        self.worldSize = 1.0
        self.setMatrixSize( 3 )
        self.initialize()


    def setWorldSize( self, size ):
        self.worldSize = size
        self.initialize()


    def setMatrixSize( self, size ):
        if size < 3:
            raise RuntimeError,\
                'Size of distance cell matrix must be at least 3'
        self.matrixSize = size
        self.initialize()

    def getSize( self ):
        return len( self.impl )

    size = property( getSize )

    def getCellSize( self ):
        return self.impl.cell_size

    cellSize = property( getCellSize )

    def initialize( self ):

        if not isinstance( self.worldSize, float ):
            raise NotImplementedError,\
                'ObjectMatrix does not support non-cubic world.' + \
                'Use SimpleObjectMatrix.'

        self.impl = object_matrix.ObjectContainer( self.worldSize, 
                                                   self.matrixSize )


    def clear( self ):
        self.initialize()


    def add( self, key, pos, radius ):

        assert radius <= self.impl.cell_size
        #assert key not in self.objectCellMap

        self.impl[ key ] = object_matrix.Sphere( pos, radius )


    def remove( self, key ):

        assert key in self.impl
        del self.impl[ key ]


    def update( self, key, pos, radius ):

        assert key in self.impl
        self.impl[ key ] = object_matrix.Sphere( pos, radius )


    def get( self, key ):
        sphere = self.impl[ key ]
        return numpy.array( [ sphere.x, sphere.y, sphere.z ] ), sphere.radius


    def getNeighborsCyclicNoSort( self, pos ):

        return self.impl.all_neighbors_array_cyclic( pos )

    def getNeighborsWithinRadiusNoSort( self, pos, radius ):

        return self.impl.neighbors_array_cyclic(\
            object_matrix.Sphere( pos, 
                                  radius ) )

    def getNeighborsCyclic( self, pos, n=None ):

        neighbors, distances = self.impl.all_neighbors_array_cyclic( pos )
        topargs = distances.argsort()[:n]

        return neighbors.take( topargs ), distances.take( topargs )


    def getNeighborsWithinRadius( self, pos, radius ):

        neighbors, distances = \
            self.impl.neighbors_array_cyclic( object_matrix.Sphere( pos, 
                                                                    radius ) )
        topargs = distances.argsort()

        return neighbors.take( topargs ), distances.take( topargs )


    def getNeighbors( self, pos, n=None ):
        return self.getNeighborsCyclic( pos, n )


    def getNeighborsNoSort( self, pos ):
        return self.getNeighborsCyclicNoSort( pos )


    '''
    def getNeighbors( self, pos, n=None ):

        neighbors, distances = self.impl.all_neighbors_array( pos )
        if n:
            topargs = distances.argsort()[:n]
        else:
            topargs = distances.argsort()
        distances = distances.take( topargs )
        neighbors = [ neighbors[arg].id for arg in topargs ]
 
        return neighbors, distancesq
        '''

    def check( self ):
        pass

