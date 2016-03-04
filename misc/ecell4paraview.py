from paraview.simple import *

meso = GetActiveSource()

renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
viewLayout1 = GetLayout()

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=meso)
programmableFilter1.Script = """
sidmap = None  # {3: 3, 4: 4}  # A serial ID mapper
r0 = 0.05  # Default radius when the input is zero

import numpy as np

inputs0 = inputs[0]

if sidmap is not None:
    sid = inputs0.RowData['sid'].copy()
    mask = np.logical_or.reduce([sid == key for key in sidmap.keys()])
    for key, value in sidmap.items():
        sid[sid == key] = value
    output.RowData.append(sid[mask], 'sid')
else:
    mask = np.ones_like(inputs0.RowData['sid'], dtype=bool)
    output.RowData.append(inputs0.RowData['sid'][mask], 'sid')

for key in ('x', 'y', 'z'):
    output.RowData.append(inputs0.RowData[key][mask], key)

if r0 is not None:
    r = inputs0.RowData['r'][mask].copy()
    r = np.where(r <= 0, r0, r)
    output.RowData.append(r, 'r')
else:
    output.RowData.append(inputs0.RowData['r'][mask], 'r')
"""
programmableFilter1.RequestInformationScript = ''
programmableFilter1.RequestUpdateExtentScript = ''
programmableFilter1.PythonPath = ''

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=programmableFilter1)
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# create a new 'Glyph'
glyph1 = Glyph(Input=tableToPoints1, GlyphType='Sphere')
glyph1.Scalars = ['POINTS', 'r']
glyph1.ScaleMode = 'scalar'
glyph1.Vectors = ['POINTS', 'None']
glyph1.ScaleFactor = 1.0
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'All Points'

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(glyph1)

# get color transfer function/color map for 'r'
sidLUT = GetColorTransferFunction('sid')

# show data in view
glyph1Display_1 = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display_1.ColorArrayName = ['POINTS', 'sid']
glyph1Display_1.LookupTable = sidLUT
glyph1Display_1.GlyphType = 'Arrow'
glyph1Display_1.SetScaleArray = ['POINTS', 'sid']
glyph1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display_1.OpacityArray = ['POINTS', 'sid']
glyph1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display_1.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

sidPWF = GetOpacityTransferFunction('sid')
sidLUT.NumberOfTableValues = 32
sidLUT.ColorSpace = 'HSV'
sidLUT.RescaleTransferFunction(0.0, 32.0)
sidPWF.RescaleTransferFunction(0.0, 32.0)

# Properties modified on renderView1
renderView1.UseGradientBackground = 1

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-5.776684363101421, 8.094366607688107, 5.7143244602859875]
renderView1.CameraFocalPoint = [2.3135256972163916, 0.550464017316699, 0.5488972440361977]
renderView1.CameraViewUp = [0.237542764231946, 0.7092276746338985, -0.6637541266873143]
renderView1.CameraParallelScale = 3.203095814674732

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
