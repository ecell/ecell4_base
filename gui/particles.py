import sys
import numpy

import matplotlib
from matplotlib.axes3d import Axes3D

import pylab

import random

class Particles:
    def __init__( self ):
        self.pos = numpy.array([[]])
        self.pos.shape = ( 0, 3 )
        self.radii = numpy.array([])
    


def loadParticles( filename ):

    file = open( filename )

    particlePools = {}

    for line in file.readlines():
        #print line
        if line[0] == '#':
            continue

        id, x, y, z, r = line.split()

        if not id in particlePools:
            particlePools[id] = Particles()

        pool = particlePools[id]
        
        pool.pos = numpy.append( pool.pos, 
                                 [[ float(x), float(y), float(z) ]],
                                 axis=0 )

        pool.radii = numpy.append( pool.radii, float( r ) )


    file.close()

    return particlePools
                        

def plotParticles( particlePools ):
    
    ax = Axes3D(pylab.figure())

    for id in particlePools.keys():
        particles = particlePools[id]
        xyz = particles.pos.transpose()
        radius = particles.radii[0]
        print id, xyz, radius

        ax.scatter3D( xyz[0], xyz[1], xyz[2], c='r', s=1 )
        #s=radius*1e8)

    pylab.show()


def test_scatter():

    ax = Axes3D(pylab.figure())
    #
    #
    n = 100
    for c,zl,zh in [('r',-50,-25),('b',-30,-5)]:
        xs,ys,zs = zip(*
                       [(random.randrange(23,32),
                         random.randrange(100),
                         random.randrange(zl,zh)
                         ) for i in range(n)])
        print xs, ys, zs
        ax.scatter3D(xs,ys,zs, c=c)
    #
    ax.set_xlabel('------------ X Label --------------------')
    ax.set_ylabel('------------ Y Label --------------------')
    ax.set_zlabel('------------ Z Label --------------------')

    pylab.show()

if __name__ == '__main__':

    filename = sys.argv[1]

    particlePools = loadParticles( filename )
    #print id, pos, r
    #test_scatter()
    plotParticles( particlePools )
