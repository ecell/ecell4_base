#!/usr/bin/python


import sys
import os
import subprocess

import vtk
import time

import numpy
import scipy.io


zoom = 1e5
size = 1e-5  * zoom
#size = 2.1544e-7  * zoom

radii = [ 1e-7, 3.2e-9, 4.02e-9 ]


colors = [ (1, .2, .2 ),  ( .2,.2,1 ), ( .8, .8, .3 )]


def addParticles( ren, positions, n ):


    for pos in positions:
        
        addParticle( ren, pos, radii[n], colors[n] )



def addParticle( ren, pos, radius, color ):
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


def writeFrame( infiles, outfile ):


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
    
    
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer( ren )
    renWin.SetSize( 400, 400 )
    #iren = vtk.vtkRenderWindowInteractor()
    #iren.SetRenderWindow(renWin)
    
    #style = vtk.vtkInteractorStyleTrackballCamera()
    #iren.SetInteractorStyle(style)

    #iren.Initialize()

    i = 0
    for filename in infiles:
        file = open( filename )
        print file
        try:
            data = scipy.io.read_array( file )
        except:
            i += 1
            continue
        
        addParticles( ren, data, 0 ) #i )
        i += 1

    ren.ResetCamera(0,size,0,size,0,size)
    camera = ren.GetActiveCamera()
    camera.Zoom(1.3)
    #camera.Dolly(1.5)
    #camera.SetDistance(.1)

    renWin.Render()
        
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput( renWin )
        
    wr = vtk.vtkPNGWriter()
    wr.SetInputConnection( w2if.GetOutputPort() )
    wr.SetFileName( outfile )
    wr.Write()
        




def main():

    tmpdir='tmpimg'
    outfile='out.mpeg'

    if not os.access( tmpdir, os.W_OK ):
        os.mkdir( tmpdir )

    i = 0
    for filename in sys.argv[1:]:
        
        outfile = os.path.join( tmpdir, '%04d.png' % i )
        print outfile
        writeFrame( (filename,), outfile )
        i += 1

    subprocess.call( ['ffmpeg', '-y', '-i',  '%s/%%04d.png' % tmpdir, '-r', '24',
                      outfile ] )

    #os.rmdir(tmpdir)

if __name__ == '__main__':
    main()

