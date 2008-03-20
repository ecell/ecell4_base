#!/usr/bin/env python


import numpy

from utils import *


class ObjectMatrixCell( object ):

    def __init__( self ):

        self.clear()


    def clear( self ):

        self.size = 0
        self.keyList = []

        self.positions = numpy.array( [], numpy.floating )
        self.positions.shape = ( 0, 3 )

        self.radii = numpy.array( [], numpy.floating )

        #self.distanceMatrix = numpy.array( [[]], numpy.floating )
        
        
    def add( self, key, pos, radius ):

        assert key not in self.keyList

        self.size += 1
        self.__resizeArrays( self.size )

        self.keyList.append( key )
        self.positions[ -1 ] = pos
        self.radii[ -1 ] = radius

    def remove( self, key ):

        self.size -= 1
        i = self.keyList.index( key )
        #i = self.keyList[ key ]

        if i != self.size:
            self.keyList[ i ] = self.keyList.pop()
            self.positions[ i ] = self.positions[ -1 ]
            self.radii[ i ] = self.radii[ -1 ]
        else:
            self.keyList.pop()

        self.__resizeArrays( self.size )


    def update( self, key, pos, radius ):
        i = self.keyList.index( key )
        self.positions[ i ] = pos
        self.radii[ i ] = radius

    def __resizeArrays( self, newsize ):

        self.positions.resize( ( newsize, 3 ) )
        self.radii.resize( newsize )

    def check( self ):
 
        assert self.size == len( self.keyList ) 
        assert self.size == len( self.positions )
        assert self.size == len( self.radii )
        '''
        for i, key in enumerate( self.keyList ):
            if ( key.pos - self.positions[i] ).sum() != 0:
                raise RuntimeError, 'keyMatrix positions consistency failed'
            if key.radius != self.radii[i]:
                raise RuntimeError, 'keyMatrix radii consistency failed'
        '''



class SimpleObjectMatrix( object ):

    def __init__( self ):

        self.matrix = ObjectMatrixCell()
        self._distanceSq = distanceSq_Cyclic
        self._distanceSqArray = distanceSqArray_Cyclic


    def setWorldSize( self, size ):
        self.worldSize = size
        self.cellSize = self.worldSize

    def setMatrixSize( self, size ):
        print 'SimpleObjectMatrix.setMatrixSize() ignored.'

    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1,\
                                                 positions ) )

    def getSize( self ):
        return self.matrix.size

    def getKeyList( self ):
        return self.matrix.keyList

    def getPositions( self ):
        return self.matrix.positions

    def getRadii( self ):
        return self.matrix.radii

    size = property( getSize )
    keyList = property( getKeyList )
    positions = property( getPositions )
    radii = property( getRadii )

    def clear( self ):
        self.matrix.clear()
        
    def add( self, key, pos, radius ):

        self.matrix.add( key, pos, radius )

    def remove( self, key ):

        self.matrix.remove( key )

    def update( self, key, pos, radius ):

        self.matrix.update( key, pos, radius )


    def getNeighbors( self, pos, n=None ):

        keyMatrix = self.matrix

        size = keyMatrix.size

        if not n:
            n = size

        distances = self.distanceArray( keyMatrix.positions, pos ) -\
           keyMatrix.radii

        topargs = distances.argsort()[:n]
        distances = distances.take( topargs )
        neighbors = [ keyMatrix.keyList[arg] for arg in topargs ]

        return neighbors, distances


    def check( self ):

        self.matrix.check()





class ObjectMatrix( object ):

    def __init__( self ):

        self.worldSize = 1.0

        self.objectCellMap = {}

        self._distanceSq = distanceSq_Cyclic
        self._distanceSqArray = distanceSqArray_Cyclic

        self.setMatrixSize( 3 )

        self.initialize()

        self.TRANSPOSES = self.allTransposes()

    def allTransposes( self ):
        transposes = []
        a = [ -1, 0, 1 ]
        for i in a:
            for j in a:
                for k in a:
                    transposes.append( [ i, j, k ] )
        return numpy.array( transposes )




    def setWorldSize( self, size ):
        self.worldSize = size
        self.initialize()

    def setMatrixSize( self, size ):
        if size < 3:
            raise RuntimeError,\
                'Size of distance cell matrix must be at least 3'


        self.matrixSize = size

        self.cellMatrix = \
            [ [ [ ObjectMatrixCell() \
                      for i in range( size ) ] \
                    for j in range( size ) ] \
                  for k in range( size ) ]

        self.initialize()

    def initialize( self ):

        self.cellSize = self.worldSize / self.matrixSize

    def hashPos( self, pos ):
        return ( ( pos % self.worldSize ) / self.cellSize ).astype( numpy.int )

    def cellByPos( self, pos ):
        i = self.hashPos( pos )
        return self.cellMatrix[ i[0] ][ i[1] ][ i[2] ]

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

        self.objectCellMap.clear()


    def add( self, key, pos, radius ):

        assert radius <= self.cellSize
        matrix = self.cellByPos( pos )
        matrix.add( key, pos, radius )

        self.objectCellMap[ key ] = matrix


    def remove( self, key ):

        matrix = self.objectCellMap[ key ]
        matrix.remove( key )

        del self.objectCellMap[ key ]


    def update( self, key, pos, radius ):

        matrix = self.objectCellMap[ key ]
        newMatrix = self.cellByPos( pos )

        if matrix is newMatrix:
            matrix.update( key, pos, radius )
        else:
            matrix.remove( key )
            newMatrix.add( key, pos, radius )
            self.objectCellMap[ key ] = newMatrix


    def getNeighborsCyclic( self, pos, n=None, dummy=None ):

        centeridx = self.hashPos( pos )

        transposes = self.TRANSPOSES + centeridx

        positions = [0,] * 27
        radii = [0,] * 27
        neighbors = []

        for i, idx in enumerate( transposes ):

            idxp = idx % self.matrixSize
            keyMatrix = self.cellMatrix[ idxp[0] ][ idxp[1] ][ idxp[2] ]
            
            # offset the positions; no need to use the cyclic distance later.
            offsetp = idxp - idx
            if offsetp[0] == offsetp[1] == offsetp[2] == 0: 
                positions[i] = keyMatrix.positions
            else:
                offset = self.cellSize * offsetp
                positions[i] = keyMatrix.positions - offset

            radii[i] = keyMatrix.radii
            neighbors += keyMatrix.keyList

        positions = numpy.concatenate( positions )
        radii = numpy.concatenate( radii )

        if len( positions ) == 0:
            return [dummy(),], [numpy.inf,]

        distances = distanceArray_Simple( positions, pos ) - radii

        if not n:
            n = len( distances )

        topargs = distances.argsort()[:n]
        distances = distances.take( topargs )
        neighbors = [ neighbors[arg] for arg in topargs ]

        return neighbors, distances


    def getNeighbors( self, pos, n=None, dummy=None ):
        return self.getNeighborsCyclic( pos, n, dummy )


    def check( self ):

        for i in self.cellMatrix:
            for j in i:
                for k in j:
                    k.check()

        if len( self.objectCellMap ) != self.size:
            raise RuntimeError, 'len( self.objectCellMap ) != self.size'


