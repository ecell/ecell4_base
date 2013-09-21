# coding; utf-8
#
# from ecell3-spatiocyte/SpatiocyteStepper.cpp
#

from math import sqrt
import struct

HCP_LATTICE = 0
CUBIC_LATTICE = 1
NVR = 0.5# the Normalized Voxel Radius
HCPl = NVR / sqrt(2)
HCPx = NVR * sqrt(8.0/3)
HCPy = NVR * sqrt(3)

latticeType = HCP_LATTICE

def coord2point(coord, row_size, layer_size):
    """
    input : coord, row_size, layer_size
    output: a tuple (x, y, z)
    """
    (grow, glayer, gcol) = coord2global(coord, row_size, layer_size)
    if latticeType == HCP_LATTICE:
        y = (gcol % 2) * HCPl + HCPy * glayer
        z = grow * 2 * NVR + ((glayer + gcol) % 2) * NVR
        x = gcol * HCPx
    elif latticeType == CUBIC_LATTICE:
        y = glayer * 2 * NVR
        z = grow * 2 * NVR
        x = gcol * 2 * NVR
    else:
        return
    point = (x, y, z)
    return point

def coord2global(coord, row_size, layer_size):
    """
    input : coord, row_size, layer_size
    output: a tuple (grow, glayer, gcol)
    """
    gcol = coord / (row_size * layer_size)
    glayer = (coord % (row_size * layer_size)) / row_size
    grow = (coord % (row_size * layer_size)) % row_size
    point = (grow, glayer, gcol)
    return point

class SpatiocyteLogReader:

    def __init__(self, logfile):
        self.logfile = open(logfile, 'rb')
        self.readInitialization()
        self.readCompVacant()
        self.headerSeek = self.tell()
        self.logfile.seek(0,2)
        self.footerSeek = self.tell()
        self.logfile.seek(self.headerSeek)

    def close(self):
        self.logfile.close()

    def tell(self):
        return self.logfile.tell()

    def isEnd(self):
        return self.tell() == self.footerSeek

    def readInitialization(self):
        '''
        corresponding to VisualizationLogProcess::initializeLog()
        '''
        self.readHeader()
        self.readLatticeSpecies()
        #self.readPolymerSpecies()
        #self.readOffLatticeSpecies()


    def readHeader(self):
        data = {}
        data['aLatticeType'] = struct.unpack('<I', self.logfile.read(4))[0]
        data['theMeanCount'] = struct.unpack('I', self.logfile.read(4))[0]
        data['aStartCoord'] = struct.unpack('I', self.logfile.read(4))[0]

        data['aRowSize'] = struct.unpack('I', self.logfile.read(4))[0]
        data['aLayerSize'] = struct.unpack('I', self.logfile.read(4))[0]
        data['aColSize'] = struct.unpack('I', self.logfile.read(4))[0]

        data['aRealRowSize'] = struct.unpack('d', self.logfile.read(8))[0]
        data['aRealLayerSize'] = struct.unpack('d', self.logfile.read(8))[0]
        data['aRealColSize'] = struct.unpack('d', self.logfile.read(8))[0]

        data['theLatticeSpSize'] = struct.unpack('I', self.logfile.read(4))[0]
        data['thePolymerSize'] = struct.unpack('I', self.logfile.read(4))[0]
        data['aResersvedSize'] = struct.unpack('I', self.logfile.read(4))[0]
        data['theOffLatticeSpSize'] = struct.unpack('I', self.logfile.read(4))[0]
        data['theLogMarker'] = struct.unpack('I', self.logfile.read(4))[0]
        data['aVoxelRadius'] = struct.unpack('d', self.logfile.read(8))[0]
        '''
        header_format = '<IIIIIIdddIIIIId'
        header_titles = ['aLatticeType', 'theMeanCount', 'aStartCoord',
                'aRowSize', 'aLayerSize', 'aColSize', 'aRealRowSize',
                'aRealLayerSize', 'aRealColSize', 'theLatticeSpSize',
                'thePolymerSize', 'aResersvedSize', 'theOffLatticeSpSize',
                'theLogMarker', 'aVoxelRadius']
        data = struct.unpack(header_format,f.read(4*19))
        '''
        self.header = data

    def getHeader(self):
        return self.header


    def readLatticeSpecies(self):
        species = []

        for i in range(self.header['theLatticeSpSize']):
            aStringSize = struct.unpack('I', self.logfile.read(4))[0]
            aString = struct.unpack(str(aStringSize) + 's',
                    self.logfile.read(aStringSize))[0]
            aRadius = struct.unpack('d', self.logfile.read(8))[0]
            print (aString, aRadius)
            species.append((aString, aRadius))

        self.header['latticeSpecies'] = species


    def readPolymerSpecies(self):
        species = []

        for i in range(self.header['thePolymerSize']):
            aRadius = struct.unpack('d', self.logfile.read(8))[0]
            species.append(aRadius)

        self.header['polymerSpecies'] = species


    def readOffLatticeSpecies(self):
        species = []
        for i in range(self.header['theOffLatticeSpSize']):
            aStringSize = struct.unpack('I', self.logfile.read(4)[0])
            aString = struct.unpack(str(aStringSize) + 's',
                    self.logfile.read(aStringSize))[0]
            aRadius = struct.unpack('d', logfile.read(8))[0]
            species.append((aString, aRadius))

        self.header['offLatticeSpecies'] = species


    def readCompVacant(self):
        '''
        corresponding to VisualizationLogProces::logCompVacant()
        '''
        data = {}

        aCurrentTime = struct.unpack('d', self.logfile.read(8))[0]
        i = 0
        data['Coords'] = {}
        for index in range(self.header['theLatticeSpSize']+1):
            i = struct.unpack('I', self.logfile.read(4))[0]
            if i == self.header['theLogMarker']:
                break
            aSize = struct.unpack('i', self.logfile.read(4))[0]
            data['Coords'][i] = []
            for j in range(aSize):
                aCoord = struct.unpack('I', self.logfile.read(4))
                data['Coords'][i].append(aCoord)
        data['Points'] = {}
        for index in range(self.header['theOffLatticeSpSize']+1):
            i = struct.unpack('I', self.logfile.read(4))[0]
            if i == self.header['theLogMarker']:
                break
            aSize = struct.unpack('i', self.logfile.read(4))[0]
            data['Points'][i] = []
            for j in range(aSize):
                (x, y, z) = struct.unpack('ddd', self.logfile.read(8*3))
                data['Points'][i].append({'x':x, 'y':y, 'z':z})

        self.header['compVacant'] = data


    def readSpecies(self):
        '''
        corresponding to VisualizationLogProcess::logSpecies()
        '''
        data = {}

        aCurrentTime = struct.unpack('d', self.logfile.read(8))[0]
        data['aCurrentTime'] = aCurrentTime

        data['Molecules'] = []
        for i in range(self.header['theLatticeSpSize']):
            molecules = self.readMolecules()
            data['Molecules'].append(molecules)

        data['SourceMolecules'] = []
        for i in range(self.header['thePolymerSize']):
            molecules = self.readSourceMolecules()
            data['SourceMolecules'].append(molecules)

        data['TargetMolecules'] = []
        for i in range(self.header['thePolymerSize']):
            molecules = self.readTargetMolecules()
            data['TargetMolecules'].append(molecules)

        data['SharedMolecules'] = []
        for i in range(self.header['thePolymerSize']):
            molecule = self.readSharedMolecules()
            data['SharedMolecules'].append(molecule)

        theLogMarker0 = struct.unpack('I', self.logfile.read(4))[0]
        if theLogMarker0 != self.header['theLogMarker']:
            print '[ERROR]\tthe log marker is different!'
            sys.exit()

        data['Polymers'] = []
        for i in range(self.header['thePolymerSize']):
            polymer = self.readPolymers()
            data['Polymers'].append(polymer)

        data['OffLattice'] = []
        for i in range(self.header['theOffLatticeSpSize']):
            offLattice = self.readOffLattice()
            data['OffLattice'].append(offLattice)

        theLogMarker1 = struct.unpack('I', self.logfile.read(4))[0]
        if theLogMarker1 != self.header['theLogMarker']:
            print '[ERROR]\tthe log marker is different!'
            sys.exit()

        return data


    def skipSpecies(self):

        self.logfile.seek(8,1)

        for i in range(self.header['theLatticeSpSize']):
            self.skipMolecules()

        for i in range(self.header['thePolymerSize']):
            self.skipSourceMolecules()

        for i in range(self.header['thePolymerSize']):
            self.skipTargetMolecules()

        for i in range(self.header['thePolymerSize']):
            self.skipSharedMolecules()

        theLogMarker0 = struct.unpack('I', self.logfile.read(4))[0]
        if theLogMarker0 != self.header['theLogMarker']:
            print '[ERROR]\tthe log marker is different!'
            sys.exit()

        for i in range(self.header['thePolymerSize']):
            polymer = self.skipPolymers()

        for i in range(self.header['theOffLatticeSpSize']):
            offLattice = self.skipOffLattice()


        theLogMarker1 = struct.unpack('I', self.logfile.read(4))[0]
        if theLogMarker1 != self.header['theLogMarker']:
            print '[ERROR]\tthe log marker is different!'
            sys.exit()

    def skipSpeciesTo(self, index):
        currentSeek = self.tell()
        self.logfile.seek(self.headerSeek)

        for i in range(index):
            self.skipSpecies()
            if self.isEnd():
                self.logfile.seek(currentSeek)
                print index," is out of bound."
                return

        return self.readSpecies()


    def readMolecules(self):
        '''
        read aSpecies->getCoord(i) i(0:aSpecies->size())
        '''
        molecules = {}
        (index, size) = struct.unpack('ii', self.logfile.read(8))
        molecules['index'] = index
        molecules['Coords'] = []
        for i in range(size):
            aCoord = struct.unpack('I', self.logfile.read(4))[0]
            molecules['Coords'].append(aCoord)
        return molecules


    def skipMolecules(self):
        self.logfile.seek(4,1)
        size = struct.unpack('i', self.logfile.read(4))[0]
        self.logfile.seek(4*size,1)


    def readSourceMolecules(self):
        '''
        read aSpecies->getSourceCoords()
        '''
        data = {}
        (aSourceIndex, aSize) = struct.unpack('ii', self.logfile.read(8))
        data['index'] = aSourceIndex
        data['Coords'] = []
        for i in range(aSize):
            aCoord = struct.unpack('I', self.logfile.read(4))[0]
            data['Coords'].append(aCoord)
        return data


    def skipSourceMolecules(self):
        self.logfile.seek(4,1)
        size = struct.unpack('i', self.logfile.read(4))[0]
        self.logfile.seek(4*size)


    def readTargetMolecules(self):
        '''
        read aSpecies->getTargetCoords()
        '''
        data = {}
        (aTargetIndex, aSize) = struct.unpack('ii', self.logfile.read(8))
        data['index'] = aTargetIndex
        data['Coords'] = []
        for i in range(aSize):
            aCoord = struct.unpack('I', self.logfile.read(4))[0]
            data['Coords'].append(aCoord)
        return data


    def readTargetMolecules(self):
        self.logfile.seek(4,1)
        size = struct.unpack('i', self.logfile.read(4))[0]
        self.logfile.seek(4*size)


    def readSharedMolecules(self):
        '''
        read aSpecies->getSharedCoords()
        '''
        data = {}
        (aSharedIndex, aSize) = struct.unpack('ii', self.logfile.read(8))
        data['index'] = aSharedIndex
        data['Coords'] = []
        for i in range(aSize):
            aCoord = struct.unpack('I', self.logfile.read(4))[0]
            data['Coords'].append(aCoord)
        return data

    def skipSharedMolecules(self):
        self.logfile.seek(4,1)
        size = struct.unpack('i', self.logfile.read(4))[0]
        self.logfile.seek(4*size)


    def readPolymers(self):
        '''
        read aSpecies->getPoint(i) i(0:aSpecies->size())
        '''
        data = {}
        (anIndex, aSize) = struct.unpack('ii', self.logfile.read(8))
        data['index'] = anIndex
        data['Points'] = []
        for i in range(aSize):
            (x, y, z) = struct.unpack('ddd', self.logfile.read(8*3))
            data['Points'].append({'x':x, 'y':y, 'z':z})
        return data

    def skipPolymers(self):
        self.logfile.seek(4,1)
        size = struct.unpack('i', self.logfile.read(4))[0]
        self.logfile.seek(24*size)


    def readOffLattice(self):
        '''
        read aSpecies->getPoint(i)
        or  aSpecies->getMultiscaleStructurePoint(i) i(0:aSpecies->size())
        '''
        data = {}
        (anIndex, aSize) = struct.unpack('ii', self.logfile.read(8))
        data['index'] = anIndex
        data['Points'] = []
        for i in range(aSize):
            (x, y, z) = struct.unpack('ddd', self.logfile.read(8*3))
            data['Points'].append({'x':x, 'y':y, 'z':z})
        return data

    def skipOffLattice(self):
        self.logfile.seek(4,1)
        size = struct.unpack('i', self.logfile.read(4))[0]
        self.logfile.seek(24*size)

