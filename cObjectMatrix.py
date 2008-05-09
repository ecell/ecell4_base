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

    def initialize( self ):

        if not isinstance( self.worldSize, float ):
            raise NotImplementedError,\
                'ObjectMatrix does not support non-cubic world.' + \
                'Use SimpleObjectMatrix.'

        self.impl = object_matrix.ObjectContainer( self.worldSize, 
                                                    self.matrixSize )



    def getSize( self ):
        return len( self.impl )

    size = property( getSize )

    def clear( self ):
        self.initialize()


    def add( self, key, pos, radius ):

        assert radius <= self.impl.cell_size
        #assert key not in self.objectCellMap

        self.impl[ key ] = object_matrix.Sphere( pos, radius )

    #def remove( self, key ):


    #def update( self, key, pos, radius ):

    def get( self, key ):
        return self.impl[ key ]


    def getNeighborsCyclic( self, pos, n=None, dummy=None ):

        n, d = self.impl.all_neighbors_array_cyclic( pos )
        return [ i.id for i in n ], d


    def getNeighborsWithinRadius( self, pos, radius ):

        n,d = self.impl.neighbors_array( pos, radius )
        return [ i.id for i in n ], d

    def getNeighbors( self, pos, n=None, dummy=None ):
        n,d = self.impl.all_neighbors_array( pos )
        return [ i.id for i in n ], d

    #def check( self ):

