import vtk
import ecell4
import os.path
import argparse
import sys

colors = [
    ('Red', [0.8, 0.1, 0.1]), 
    ('Green', [0.27, 0.8, 0.21]), 
    ('Blue', [0.24, 0.41, 0.7]), 
    ('Yellow', [1.0, 0.5, 0.0]), 
    ('Orange', [1.0, 0.37, 0.05]), 
    ('White', [1, 1, 1]), 
    ('Magenta', [0.72, 0.29, 1.0]), 
    ('Cyan', [0.1, 1.0, 0.6]), 
    ('Black', [0.1, 0.1, 0.1]), 
    ('Grey', [0.46, 0.46, 0.46]), 
    ('LightBlue', [0.32, 0.42, 1]), 
    ('LightGrey', [0.6, 0.6, 0.6]), 
    ('BrightGreen', [0.4, 1.0, 0.14]), 
    ('BrightYellowGreen', [0.64, 1.0, 0.05]), 
    ('BrightYellow', [1.0, 0.67, 0.0]), 
    ('WhiteGray', [0.9, 0.9, 0.9]), 
    ('WhiteMagenta', [0.8, 0.48, 1.0]), 
    ('WhiteYellow', [1.0, 0.75, 0.17]), 
    ('WhitePurple', [0.67, 0.6, 1.0]), 
    ('DarkRed', [0.46, 0.1, 0.1]), 
    ('DarkGreen', [0.1, 0.5, 0.1]), 
    ('DarkBlue', [0.1, 0.2, 0.5]), 
    ('DarkOrange', [0.845, 0.179, 0.102])]

class vtkTimerCallback():

    def __init__(self, filenames, saveimage=False, txt=None):
         self.timer_count = 0
         self.filenames = filenames
         self.timer_id = None
         self.txt = txt
         self.saveimage = saveimage

    def filename(self):
        return self.filenames[min(self.timer_count, len(self.filenames) - 1)]

    def stop(self, iren):
        if self.timer_id is not None:
            iren.DestroyTimer(self.timer_id)

    def update(self, renWin):
        if (len(self.filenames) <= self.timer_count
            or not os.path.isfile(self.filenames[self.timer_count])):
            return False
        filename = self.filename()

        renWin.SetWindowName(filename)
        ren = renWin.GetRenderers().GetFirstRenderer()
        remove_all_actors(ren)
        w = load_world(filename)
        add_actors(w, ren, source)
        if self.txt is not None:
            self.txt.SetInput("t={0:g}".format(w.t()))

        # camera = ren.GetActiveCamera()
        # print(camera.GetPosition(), camera.GetFocalPoint())
        # # (0.03282041493100417, -0.04556241788024047, 1.2086963582474413)
        # # (0.0, 0.0, 0.0)

        renWin.Render()

        if self.saveimage:
            screenshot(renWin, os.path.splitext(filename)[0] + '.png')

        self.timer_count += 1
        return True

    def execute(self, iren, event=None):
        renWin = iren.GetRenderWindow()
        if not self.update(renWin):
            self.stop(iren)
            return

def remove_all_actors(ren):
    actors = ren.GetActors()
    while actors.GetNumberOfItems() > 0:
        actor = actors.GetLastActor()
        ren.RemoveActor(actor)
        actors.RemoveItem(actor)

def load_world(filename):
    return ecell4.meso.MesoscopicWorld(filename)

def add_actors(w, ren, source):
    ren.SetBackground(0.2, 0.3, 0.4)

    L = max(w.edge_lengths())
    shift = w.edge_lengths() / L * 0.5

    for cnt, sp in enumerate(w.list_species()):
        particles = w.list_particles(sp)
        points = vtk.vtkPoints()
        points.SetNumberOfPoints(len(particles))
        for i, (pid, p) in enumerate(particles):
            pos = tuple(p.position() / L - shift)
            points.SetPoint(i, pos[0], pos[1], pos[2])

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)

        atoms = vtk.vtkGlyph3D()
        atoms.SetInput(polydata)
        atoms.SetSource(source.GetOutput())
        atoms.SetScaleFactor(1.0)
        atoms.SetScaleModeToScaleByScalar()

        # mapper
        mapper = vtk.vtkPolyDataMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.SetInput(atoms.GetOutput())
        else:
            mapper.SetInputConnection(atoms.GetOutputPort())
         
        # actor
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(*colors[cnt % len(colors)][1])
        actor.SetMapper(mapper)
         
        # assign actor to the renderer
        ren.AddActor(actor)

def screenshot(renWin, filename):
    # ffmpeg -r 15 -i input%03d.png -qscale 0 output.mp4
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(renWin)
    w2if.Update()
 
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(filename)
    writer.SetInput(w2if.GetOutput())
    writer.Write()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Visualizing HDF5s with VTK.')
    parser.add_argument(
        'filenames', metavar='filename', type=str, nargs='+',
        help='HDF5 filenames for visualizing')
    parser.add_argument(
        '--radius', '-r', metavar='radius', type=float, default=0.002,
        help='a default radius rescaled')
    parser.add_argument(
        '--resolution', '-R', metavar='res', type=int, default=None,
        help='a resolution for a sphere')
    parser.add_argument(
        '--with-nolabel', action='store_true', dest='nolabel',
        help='whether displaying a label for the time or not')
    parser.add_argument(
        '--save-image', action='store_true', dest='saveimage',
        help='whether saving a screenshot for each view')
    parser.add_argument(
        '--offscreen', action='store_true', dest='offscreen',
        help='enable offscreen rendering')
    parser.add_argument(
        '--interval', metavar='dt', type=int, default=100,
        help='an interval to switch and update HDF5 files in milliseconds')
    args = parser.parse_args()

    filenames = args.filenames
    filenames.sort()
    res = args.resolution
    nolabel = args.nolabel
    dt = args.interval
    radius = args.radius
    saveimage = args.saveimage
    offscreen = args.offscreen
    if offscreen:
        saveimage = True

    # create a rendering window and renderer
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize(768, 432)
    if offscreen:
        renWin.SetOffScreenRendering(True)
     
    ren = vtk.vtkRenderer()
    renWin.AddRenderer(ren)

    # create source
    if res is None:
        source = vtk.vtkPointSource()
        source.SetRadius(radius)
    else:
        source = vtk.vtkSphereSource()
        source.SetThetaResolution(res)
        source.SetPhiResolution(res)
        source.SetRadius(radius)

    if not nolabel:
        txt = vtk.vtkTextActor()
        txt.SetInput("t={0:g}".format(0.0))
        txtprop = txt.GetTextProperty()
        txtprop.SetFontFamilyToArial()
        txtprop.SetFontSize(18)
        txtprop.SetColor(1, 1, 1)
        txt.SetDisplayPosition(20, 30)
        ren.AddActor(txt)

    if not offscreen:
        # create a renderwindowinteractor
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
         
        # enable user interface interactor
        iren.Initialize()

    # Sign up to receive TimerEvent
    cb = vtkTimerCallback(filenames, saveimage, None if nolabel else txt)

    if not offscreen:
        cb.execute(iren)
        if len(filenames) > 1:
            iren.AddObserver('TimerEvent', cb.execute)
            cb.timer_id = iren.CreateRepeatingTimer(dt)
        #start the interaction and timer
        iren.Start()
    else:
        while cb.update(renWin):
            print(cb.filename())
