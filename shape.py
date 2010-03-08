import numpy
from utils import *
import math

class Shape(object):
    def __init__(self):
        pass

    def distanceTo(self, pos):
        """Note: cyclicTranspose 'pos' before you use this.

        """
        return abs(self.signedDistanceTo(pos))


class Cylinder(Shape):
    def __init__(self, origin, radius, orientation, size):
        Shape.__init__(self)
        # Origin is the centre of the cylinder. This has the slight advantage 
        # over other options that we can make use of symmetry sometimes.
        self.origin = numpy.array(origin)
        self.radius = radius
        self.unitZ = normalize(numpy.array(orientation))
        # Size is the half length of the cylinder.
        self.size = size                            
        self.vectorZ = self.unitZ * size # Extra.

        # Select basis vector in which self.unitZ is smallest.
        _, basisVector = min(zip(abs(self.unitZ), 
                                 [[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
        # Find (any) 2 vectors perpendicular to self.unitZ.
        self.unitX = normalize(numpy.cross(self.unitZ, basisVector))
        self.unitY = normalize(numpy.cross(self.unitZ, self.unitX))

    def getSize(self):
        return self._size
    def setSize(self, size):
        self._size = size
    size = property(getSize, setSize)

    def signedDistanceTo(self, pos):
        """Note: cyclicTranspose 'pos' before you use this.

        """
        r, z, = self.toInternal(pos)
        dz = abs(z) - self.size
        dr = r - self.radius
        if dz > 0:
            # pos is (either) to the right or to the left of the cylinder.
            if dr > 0:
                # Compute distance to edge.
                distance = math.sqrt(dz * dz + dr * dr)
            else:
                distance = dz
        else:
            if dr > 0:
                # pos is somewhere 'parallel' to the cylinder.
                distance = dr
            else:
                # Inside cylinder, dz and dr are negative.
                distance = max(dz, dr)
        return distance

    def toInternal(self, pos):
        """Return the (z, r) components of pos in a coordinate system 
        defined by the vectors unitR and unitZ, where unitR is choosen such 
        that unitR and unitZ define a plane in which pos lies.

        """
        posVector = pos - self.origin

        z = numpy.dot(posVector, self.unitZ) # Can be <0.
        posVectorZ = z * self.unitZ

        posVectorR = posVector - posVectorZ
        r = length(posVectorR)       # Always >= 0.
          
        return r, z

    def projectedPoint(self, pos):
        """Returns:
        1. the position of the projected point of 'pos' onto the main axis of 
        the cylinder.
        2. the distance (always positive) between that point and 'pos'.

        """
        r, z = self.toInternal(pos)
        return self.origin + z * self.unitZ, r

    def __str__(self):
        return("Cylinder: " + str(self.origin) + " " + str(self.radius) +
               " " + str(self.unitZ) + " " + str(self.size))


class Box(Shape):
    def __init__(self, origin, vectorX, vectorY, vectorZ, Lx, Ly, Lz):
        Shape.__init__(self)
        self.origin = numpy.array(origin)

        self.unitX = normalize(numpy.array(vectorX))
        self.unitY = normalize(numpy.array(vectorY))
        self.unitZ = normalize(numpy.array(vectorZ))

        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz

        # vectorX, vectorY and vectorZ didn't have to be normalized before.
        self.vectorX = self.unitX * Lx
        self.vectorY = self.unitY * Ly
        self.vectorZ = self.unitZ * Lz

    def signedDistanceTo(self, pos):
        """Note: cyclicTranspose 'pos' before you use this.

        """
        x, y, z = self.toInternal(pos)
        dx = abs(x) - self.Lx
        dy = abs(y) - self.Ly
        dz = abs(z) - self.Lz

        if dx > 0:
            if dy > 0:
                if dz > 0:
                    distance = sqrt(dx * dx + dy * dy + dz * dz)
                else:
                    distance = sqrt(dx * dx + dy * dy)
            else:
                if dz > 0:
                    distance = sqrt(dx * dx + dz * dz)
                else:
                    distance = dx
        else:
            if dy > 0:
                if dz > 0:
                    distance = sqrt(dy * dy + dz * dz)
                else:
                    distance = dy
            else:
                if dz > 0:
                    distance = dz
                else:
                    # Inside box. Pick negative distance closest to 0.
                    distance = max(max(dx, dy), dz)

        return distance

    def toInternal(self, pos):
        """First compute the (x, y, z) components of pos in a coordinate 
        system defined by the vectors unitX, unitY, unitZ of the box.

        """
        posVector = pos - self.origin
        x = numpy.dot(posVector, self.unitX)
        y = numpy.dot(posVector, self.unitY)
        z = numpy.dot(posVector, self.unitZ)

        return x, y, z

    def projectedPoint(self, pos):
        """Returns:
        1. the position of the projected point of 'pos' onto the xy-plane of 
        the box.
        2. the distance (positive or negative) between that point and 'pos'.

        Note: cyclicTranspose 'pos' if you are going to use the 'z' value.

        """
        x, y, z = self.toInternal(pos)
        return self.origin + x * self.unitX + y * self.unitY, z

    def __str__(self):
        return("Box: " + str(self.origin) + " " + str(self.vectorX) + " " +
               str(self.vectorY) + " " + str(self.vectorZ))


