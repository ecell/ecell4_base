#ifndef ECELL4_POLYGON_HDF5_WRITER_HPP
#define ECELL4_POLYGON_HDF5_WRITER_HPP

#include <cstring>
#include <memory>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Triangle.hpp"

namespace ecell4
{

class Polygon; // forward declaration

void save_polygon_hdf5(const Polygon&,  H5::Group*);
void load_polygon_hdf5(const H5::Group&,  Polygon*);

} // ecell4
#endif// ECELL4_POLYGON_HDF5_WRITER_HPP
