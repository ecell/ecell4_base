#!/usr/bin/env python


import numpy

from utils import *


class ObjectMatrixCell( object ):

    def __init__( self ):

        self.clear()


    def clear( self ):

        self.size = 0
        self.objList = []

        self.positions = numpy.array( [], numpy.floating )
        self.positions.shape = ( 0, 3 )

        self.radii = numpy.array( [], numpy.floating )

        #self.distanceMatrix = numpy.array( [[]], numpy.floating )
        
        
    def add( self, obj ):

        assert obj not in self.objList

        self.size += 1
        self.__resizeArrays( self.size )

        self.objList.append( obj )
        self.positions[ -1 ] = obj.pos
        self.radii[ -1 ] = obj.radius

    def remove( self, obj ):

        self.size -= 1
        i = self.objList.index( obj )

        if i != self.size:
            self.objList[ i ] = self.objList.pop()
            self.positions[ i ] = self.positions[ -1 ]
            self.radii[ i ] = self.radii[ -1 ]
        else:
            self.objList.pop()

        self.__resizeArrays( self.size )


    def update( self, obj ):
        i = self.objList.index( obj )
        self.positions[ i ] = obj.pos
        self.radii[ i ] = obj.radius

    def __resizeArrays( self, newsize ):

        self.positions.resize( ( newsize, 3 ) )
        self.radii.resize( newsize )

    def check( self ):
 
        assert self.size == len( self.objList ) 
        assert self.size == len( self.positions )
        assert self.size == len( self.radii )

        for i, obj in enumerate( self.objList ):
            if ( obj.pos - self.positions[i] ).sum() != 0:
                raise RuntimeError, 'objMatrix positions consistency failed'
            if obj.radius != self.radii[i]:
                raise RuntimeError, 'objMatrix radii consistency failed'




class SimpleObjectMatrix( object ):

    def __init__( self ):

        self.matrix = ObjectMatrixCell()
        self._distanceSq = distanceSq_Cyclic
        self._distanceSqArray = distanceSqArray_Cyclic


    def setWorldSize( self, size ):
        self.worldSize = size

    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1,\
                                                 positions ) )

    def getSize( self ):
        return self.matrix.size

    def getObjList( self ):
        return self.matrix.objList

    def getPositions( self ):
        return self.matrix.positions

    def getRadii( self ):
        return self.matrix.radii

    size = property( getSize )
    objList = property( getObjList )
    positions = property( getPositions )
    radii = property( getRadii )

    def clear( self ):
        self.matrix.clear()
        
    def add( self, obj ):

        self.matrix.add( obj )

    def remove( self, obj ):

        self.matrix.remove( obj )

    def update( self, obj ):

        self.matrix.update( obj )


    def getNeighbors( self, pos, n=None ):

        objMatrix = self.matrix

        size = objMatrix.size

        if not n:
            n = size

        distances = self.distanceArray( objMatrix.positions, pos ) -\
           objMatrix.radii

        topargs = distances.argsort()[:n]
        distances = distances.take( topargs )
        neighbors = [ objMatrix.objList[arg] for arg in topargs ]

        return neighbors, distances


    def check( self ):

        self.matrix.check()





class ObjectMatrix( object ):

    def __init__( self ):

        self.worldSize = 1.0

        self.objCellMap = {}

        self._distanceSq = distanceSq_Cyclic
        self._distanceSqArray = distanceSqArray_Cyclic

        self.setMatrixSize( 1 )

        self.initialize()

    def setWorldSize( self, size ):
        self.worldSize = size
        self.initialize()

    def setMatrixSize( self, size ):
        self.matrixSize = size

        self.cellMatrix = [ [ [ ObjectMatrixCell(), ] * size ] * size ] * size

        self.initialize()

    def initialize( self ):

        self.cellSize = self.worldSize / self.matrixSize

    def hashPos( self, pos ):
        #print pos
        #assert pos.max() < self.worldSize and pos.min() >= 0.0
        return ( pos % self.worldSize / self.cellSize ).astype( numpy.int )

    def cellByPos( self, pos ):
        i = self.hashPos( pos )
        return self.cellMatrix[ i[0] ][ i[1] ][ i[2] ]

    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1,\
                                                 positions ) )

    def getSize( self ):
        size = 0
        for i in self.cellMatrix:
            for j in i:
                for k in j:
                    size += k.size
        return size


    size = property( getSize )

    def clear( self ):

        for i in self.cellMatrix:
            for j in i:
                for k in j:
                    k.clear()

        self.objCellMap.clear()


    def add( self, obj ):

        assert obj.radius <= self.cellSize
        matrix = self.cellByPos( obj.pos )
        matrix.add( obj )

        self.objCellMap[ obj ] = matrix


    def remove( self, obj ):

        matrix = self.objCellMap[ obj ]
        matrix.remove( obj )

        #assert matrix == self.cellByPos( obj.pos )

        del self.objCellMap[ obj ]


    def update( self, obj ):

        matrix = self.objCellMap[ obj ]
        newMatrix = self.cellByPos( obj.pos )

        if matrix is newMatrix:
            matrix.update( obj )
        else:
            matrix.remove( obj )
            newMatrix.add( obj )
            self.objCellMap[ obj ] = newMatrix


    def getNeighbors( self, pos, n=None ):

        idx = self.hashPos( pos )

        objMatrix = self.cellMatrix[ idx[0] ][ idx[1] ][ idx[2] ]

        size = objMatrix.size

        if not n:
            n = size

        distances = self.distanceArray( objMatrix.positions, pos ) -\
           objMatrix.radii

        topargs = distances.argsort()[:n]
        distances = distances.take( topargs )
        neighbors = [ objMatrix.objList[arg] for arg in topargs ]

        return neighbors, distances


    def check( self ):

        for i in self.cellMatrix:
            for j in i:
                for k in j:
                    k.check()

        if len( self.objCellMap ) != self.size:
            raise RuntimeError, 'len( self.objCellMap ) != self.size'


