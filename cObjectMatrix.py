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
        print radius, self.impl.cell_size, self.impl.world_size

        assert radius <= self.impl.cell_size
        #assert key not in self.objectCellMap

        self.impl[ key ] = object_matrix.Sphere( pos, radius )


    def remove( self, key ):
        print 'del', key
        assert key in self.impl
        del self.impl[ key ]


    def update( self, key, pos, radius ):
        print 'update', key
        assert key in self.impl
        del self.impl[ key ]
        self.impl[ key ] = object_matrix.Sphere( pos, radius )


    def get( self, key ):
        return self.impl[ key ]


    def getNeighborsCyclic( self, pos, n=None ):

        neighbors, distances = self.impl.all_neighbors_array_cyclic( pos )
        if n:
            topargs = distances.argsort()[:n]
        else:
            topargs = distances.argsort()
        distances = distances.take( topargs )
        neighbors = [ neighbors[arg].id for arg in topargs ]
        print neighbors, distances
        return neighbors, distances


    def getNeighborsWithinRadius( self, pos, radius ):

        neighbors, distances = \
            self.impl.neighbors_array_cyclic( object_matrix.Sphere( pos, 
                                                                    radius ) )
        topargs = distances.argsort()
        distances = distances.take( topargs )
        neighbors = [ neighbors[arg].id for arg in topargs ]
 
        return neighbors, distances

    def getNeighbors( self, pos, n=None ):
        return self.getNeighborsCyclic( pos, n )

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

    #def check( self ):

