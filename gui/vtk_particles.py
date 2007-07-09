#!/usr/bin/python


import sys

import vtk
import time

import numpy
import scipy.io


zoom = 1e5
size = 1e-5  * zoom
#size = 2.1544e-7  * zoom

radii = [ 1e-7, 3.2e-9, 4.02e-9 ]


colors = [ (1, .2, .2 ),  ( .2,.2,1 ), ( .8, .8, .3 )]


def renderParticles( ren, positions, n ):


    for pos in positions:
        
        renderParticle( pos, radii[n], colors[n] )



def renderParticle( pos, radius, color ):
        vpos = pos * zoom
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter( vpos[0],vpos[1],vpos[2] )
        sphere.SetRadius( radius * zoom )
        sphereMapper = vtk.vtkPolyDataMapper()
        sphereMapper.SetInputConnection( sphere.GetOutputPort() )

        sphereActor = vtk.vtkActor()
        sphereActor.SetMapper( sphereMapper )
        sphereActor.GetProperty().SetColor( color )
    
        ren.AddActor( sphereActor )




ren = vtk.vtkRenderer()
ren.SetBackground( 1, 1, 1 )

cube = vtk.vtkCubeSource()
cube.SetBounds(0,size,0,size,0,size)
cubeMapper = vtk.vtkPolyDataMapper()
cubeMapper.SetInputConnection( cube.GetOutputPort() )
cubeActor = vtk.vtkActor()
cubeActor.SetMapper( cubeMapper )
cubeActor.GetProperty().SetOpacity( 0.2 )

ren.AddActor( cubeActor )


i = 0
for filename in sys.argv[1:]:
    file = open( filename )
    print file
    try:
        data = scipy.io.read_array( file )
    except:
        i += 1
        continue

    renderParticles( ren, data, i )
    i += 1



renWin = vtk.vtkRenderWindow()
renWin.AddRenderer( ren )
renWin.SetSize( 600, 600 )


iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

style = vtk.vtkInteractorStyleTrackballCamera()
iren.SetInteractorStyle(style)

iren.Initialize()
iren.Start()
