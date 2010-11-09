# http://www.paraview.org/Wiki/Python_Programmable_Filter
#
# This script generates a helix double.
# This is intended as the script of a 'Programmable Source'.
#
# Note: it won't run standalone.

import math
 
num_pts = 160               # Points along each helix. 80. 
length = 1.0                # Length of each helix. 8.
rounds = 6.0                # Number of times around. 3.
phase_shift = math.pi/1.5   # Phase shift between helixes.
scaleyz = 0.5
if 'HELIX_RADIUS_SCALE_FACTOR' in vars():
    scaleyz *= HELIX_RADIUS_SCALE_FACTOR

# Get a vtk.PolyData object for the output
#pdo = vtk.vtkPolyData()
pdo = self.GetPolyDataOutput()
 
# This will store the points for the helix.
new_pts = vtk.vtkPoints()
for i in range(0, num_pts):
   # Generate points for first helix.
   x = i*length/num_pts - length/2
   y = scaleyz * math.sin(- i*rounds*2*math.pi/num_pts)
   z = scaleyz * math.cos(- i*rounds*2*math.pi/num_pts)
   # Stupid Paraview wants a normal vector to the cylinder to orient
   # it. So x and y swapped.
   new_pts.InsertPoint(i, y, x, z)
 
   # Generate Points for second helix. Add a phase offset to y and z.
   y = scaleyz * math.sin(- i*rounds*2*math.pi/num_pts+phase_shift)
   z = scaleyz * math.cos(- i*rounds*2*math.pi/num_pts+phase_shift)
   # Offset helix 2 pts by 'num_pts' to keep separate from helix 1 Pts.
   # Stupid Paraview wants a normal vector to the cylinder to orient
   # it. So x and y swapped.
   new_pts.InsertPoint(i+num_pts, y,x,z)
 
# Add the points to the vtk_poly_data object.
pdo.SetPoints(new_pts)
 
# Make two vtkPolyLine objects to hold curve construction data.
aPolyLine1 = vtk.vtkPolyLine()
aPolyLine2 = vtk.vtkPolyLine()
 
# Indicate the number of points along the line.
aPolyLine1.GetPointIds().SetNumberOfIds(num_pts)
aPolyLine2.GetPointIds().SetNumberOfIds(num_pts)
for i in range(0,num_pts):
   # First Helix - use the first set of points.
   aPolyLine1.GetPointIds().SetId(i, i)
   # Second Helix - use the second set of points.
   # (Offset the point reference by 'num_pts').
   aPolyLine2.GetPointIds().SetId(i,i+num_pts)
 
# Allocate the number of 'cells' that will be added. Two 'cells' for the helix 
# curves, and one 'cell' for every 3rd point along the helixes.
links = range(0,num_pts,3)
pdo.Allocate(2+len(links), 1)
 
# Add the poly line 'cell' to the vtk_poly_data object.
pdo.InsertNextCell(aPolyLine1.GetCellType(), aPolyLine1.GetPointIds())
pdo.InsertNextCell(aPolyLine2.GetCellType(), aPolyLine2.GetPointIds())
 
for i in links:
   # Add a line connecting the two helixes.
   aLine = vtk.vtkLine()
   aLine.GetPointIds().SetId(0, i)
   aLine.GetPointIds().SetId(1, i+num_pts)
   pdo.InsertNextCell(aLine.GetCellType(), aLine.GetPointIds())

