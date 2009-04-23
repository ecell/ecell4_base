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

        return self.impl.size()

    size = property( getSize )


    def getCellSize( self ):
        return self.impl.cell_size

    cellSize = property( getCellSize )


    def initialize( self ):

        self.impl = object_matrix.ObjectContainer( self.worldSize, 
                                                   self.matrixSize )



    def clear( self ):
        self.initialize()


    def add( self, key, pos, radius ):

        assert radius < self.cellSize * .5
        assert not self.impl.contains( key )

        self.impl.update( key, pos, radius )


    def remove( self, key ):

        assert self.impl.contains( key )

        self.impl.erase( key )



    def update( self, key, pos, radius ):

        assert radius < self.cellSize * .5
        assert self.impl.contains( key )

        self.impl.update( key, pos, radius )


    def get( self, key ):

        return self.impl.get( key )


    def getNeighborsCyclicNoSort( self, pos ):

        return self.impl.all_neighbors_array_cyclic( pos )


    def getNeighborsWithinRadiusNoSort( self, pos, radius ):

        assert radius < self.cellSize * .5

        return self.impl.neighbors_array_cyclic( pos, radius )


    def getNeighborsCyclic( self, pos, n=None ):

        neighbors, distances = self.impl.all_neighbors_array_cyclic( pos )
        topargs = distances.argsort()[:n]

        return neighbors.take( topargs ), distances.take( topargs )


    def getNeighborsWithinRadius( self, pos, radius ):

        assert radius < self.cellSize * .5

        neighbors, distances = \
            self.impl.neighbors_array_cyclic( pos, radius )

        topargs = distances.argsort()

        return neighbors.take( topargs ), distances.take( topargs )


    def getNeighbors( self, pos, n=None ):

        return self.getNeighborsCyclic( pos, n )


    def getNeighborsNoSort( self, pos ):

        return self.getNeighborsCyclicNoSort( pos )


    def check( self ):

        pass

