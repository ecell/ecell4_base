from types cimport Real
from core cimport Cpp_PlanarSurface

cdef extern from "ecell4/core/PlanarSurface.hpp" namespace "ecell4":
    Cpp_PlanarSurface create_x_plane(Real)
    Cpp_PlanarSurface create_y_plane(Real)
    Cpp_PlanarSurface create_z_plane(Real)
