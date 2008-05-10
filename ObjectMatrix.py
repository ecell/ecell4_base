#!/usr/bin/env python


import numpy

from utils import *


class ObjectMatrixCell( object ):

    def __init__( self ):

        self.clear()


    def clear( self ):

        self.size = 0
        self.keyList = []
        self.indexMap = {}

        self.positions = numpy.array( [], numpy.floating )
        self.positions.shape = ( 0, 3 )

        self.radii = numpy.array( [], numpy.floating )

        #self.distanceMatrix = numpy.array( [[]], numpy.floating )
        
        
    def add( self, key, pos, radius ):

        assert key not in self.keyList
        assert key not in self.indexMap

        self.keyList.append( key )
        self.indexMap[ key ] = self.size

        self.size += 1
        self.__resizeArrays( self.size )

        self.positions[ -1 ] = pos
        self.radii[ -1 ] = radius

    def remove( self, key ):

        self.size -= 1
        i = self.indexMap[key]

        if i != self.size:
            k = self.keyList.pop()
            self.keyList[ i ] = k
            self.indexMap[k] = i
            self.positions[ i ] = self.positions[ -1 ]
            self.radii[ i ] = self.radii[ -1 ]
        else:
            self.keyList.pop()

        del self.indexMap[key]

        self.__resizeArrays( self.size )


    def update( self, key, pos, radius ):

        i = self.indexMap[key]
        self.positions[ i ] = pos
        self.radii[ i ] = radius

    def get( self, key ):

        i = self.indexMap[key]
        return self.positions[ i ], self.radii[ i ]


    def __resizeArrays( self, newsize ):

        self.positions.resize( ( newsize, 3 ) )
        self.radii.resize( newsize )

    def check( self ):
 
        assert self.size == len( self.keyList ) 
        assert self.size == len( self.positions )
        assert self.size == len( self.radii )



class SimpleObjectMatrix( object ):

    def __init__( self ):

        self.matrix = ObjectMatrixCell()
        self._distanceSq = distanceSq_Cyclic
        self._distanceSqArray = distanceSqArray_Cyclic


    def setWorldSize( self, size ):
        self.worldSize = size
        if isinstance( size, float ):
            self.cellSize = self.worldSize
        else:
            self.cellSize = numpy.max( self.worldSize )

    def setMatrixSize( self, size ):
        print 'SimpleObjectMatrix.setMatrixSize() ignored. (', size , ')'

    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1, positions ) )

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

    def get( self, key ):
        return self.matrix.get( key )


    def getNeighbors( self, pos, n=None, dummy=None ):

        matrix = self.matrix

        if not n:
            n = matrix.size

        if len( matrix.positions ) == 0:
            return [], []

        distances = self.distanceArray( matrix.positions, pos ) - matrix.radii

        topargs = distances.argsort()[:n]
        distances = distances.take( topargs )
        neighbors = [ matrix.keyList[arg] for arg in topargs ]

        return neighbors, distances


    def getNeighborsWithinRadius( self, pos, radius ):

        matrix = self.matrix

        if len( matrix.positions ) == 0:
            return [], []

        distances = self.distanceArray( matrix.positions, pos ) - matrix.radii

        args = ( distances < radius ).nonzero()[0]
        distances = distances.take( args )
        neighbors = [ matrix.keyList[arg] for arg in args ]

        args = distances.argsort()

        distances = distances.take( args )
        neighbors = [ neighbors[arg] for arg in args ]

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

    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1, positions ) )
        

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
        #numpy.array(\
        self.cellMatrix =    \
            [ [ [ ObjectMatrixCell() \
                      for i in range( size ) ] \
                    for j in range( size ) ] \
                  for k in range( size ) ] #)

        self.initialize()

    def initialize( self ):

        if not isinstance( self.worldSize, float ):
            raise NotImplementedError,\
                'ObjectMatrix does not support non-cubic world.' + \
                'Use SimpleObjectMatrix.'

        self.cellSize = self.worldSize / self.matrixSize

    def hashPos( self, pos ):
        return ( ( pos % self.worldSize ) / self.cellSize ).astype( numpy.int )

    def cellByPos( self, pos ):
        i = self.hashPos( pos )
        return self.cellMatrix[ i[0] ][ i[1] ][ i[2] ]

    def getKeyList( self ):
        return self.objectCellMap.keys()

    def getSize( self ):
        return len( self.objectCellMap )

    keyList = property( getKeyList )
    size = property( getSize )

    def clear( self ):

        for i in self.cellMatrix:
            for j in i:
                for k in j:
                    k.clear()

        self.objectCellMap.clear()


    def add( self, key, pos, radius ):

        assert radius <= self.cellSize
        assert key not in self.objectCellMap
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

    def get( self, key ):
        matrix = self.objectCellMap[ key ]
        return matrix.get( key )


    def getNeighborsCyclic( self, pos, n=None, dummy=None ):

        centeridx = self.hashPos( pos )

        positions = []
        radii = []
        neighbors = []

        transposes = self.TRANSPOSES + centeridx
        for idx in transposes:

            idxp = idx % self.matrixSize
            #matrix = self.cellMatrix[idxp[0],idxp[1],idxp[2]]
            matrix = self.cellMatrix[idxp[0]][idxp[1]][idxp[2]]
            if matrix.size == 0:
                continue
            
            # offset the positions; no need to use the cyclic distance later.
            offsetp = idxp - idx
            if offsetp[0] == offsetp[1] == offsetp[2] == 0: 
                positions += [matrix.positions,]
            else:
                offset = self.cellSize * offsetp
                positions += [matrix.positions - offset,]

            radii += [ matrix.radii, ]
            neighbors += matrix.keyList

        if len( neighbors ) == 0:
            return [dummy,], [INF,]

        positions = numpy.concatenate( positions )
        radii = numpy.concatenate( radii )
        if len( neighbors ) == 0:
            return [dummy,], [numpy.inf,]

        distances = distanceArray_Simple( positions, pos ) - radii

        if not n:
            n = len( distances )

        topargs = distances.argsort()[:n]
        distances = distances.take( topargs )
        neighbors = [ neighbors[arg] for arg in topargs ]

        return neighbors, distances


    def getNeighborsWithinRadius( self, pos, radius ):

        centeridx = self.hashPos( pos )

        distances = []
        neighbors = []

        transposes = ( self.TRANSPOSES + centeridx ) % self.matrixSize
        for idx in transposes:

            matrix = self.cellMatrix[idx[0]][idx[1]][idx[2]]
            #matrix = self.cellMatrix[idx[0],idx[1],idx[2]]
            if matrix.size == 0:
                continue
            dists = self.distanceArray( matrix.positions, pos ) - matrix.radii
            args = ( dists < radius ).nonzero()[0]

            if len( args ) != 0:
                distances += [ dists.take( args ), ]
                neighbors += [ matrix.keyList[arg] for arg in args ]

        if len( neighbors ) == 0:
            return [], []

        distances = numpy.concatenate( distances )

        # sort again
        args = distances.argsort()
        distances = distances.take( args )
        neighbors = [ neighbors[arg] for arg in args ]

        return neighbors, distances


    def getNeighbors( self, pos, n=None, dummy=None ):
        return self.getNeighborsCyclic( pos, n, dummy )


    def check( self ):

        for i in self.cellMatrix:
            for j in i:
                for k in j:
                    k.check()

        size = 0
        for i in self.cellMatrix:
            for j in i:
                for k in j:
                    size += k.size

        if size != self.size:
            raise RuntimeError, 'size consistency broken.'


