from object_matrix import *

c = ObjectContainer(1.0, 10)
assert c.matrix_size == 10
assert c.cell_size == 0.1
assert len(c) == 0
c[0] = Sphere(0.5, 0.3, 0.2, 0.1)
assert len(c) == 1
assert c[1] == None
c[1] = Sphere(0.0, 0.3, 0.9, 0.1)
assert len(c) == 2
assert c[2] == None
assert isinstance(c[0], SphereRef)
for i in c.iterneighbors(Sphere(0.45, 0.23, 0.13, 0.09)):
    print i
for i in c.iterneighbors_cyclic(Sphere(0.9, 0.23, 0.0, 0.1)):
    print i
