#ifndef ECELL4_POLYGON_HDF5_WRITER_HPP
#define ECELL4_POLYGON_HDF5_WRITER_HPP

#include <cstring>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Triangle.hpp"

namespace ecell4
{

class Polygon; // forward declaration

struct PolygonHDF5Traits
{
    struct h5_triangle_struct {
        double p1x;
        double p1y;
        double p1z;
        double p2x;
        double p2y;
        double p2z;
        double p3x;
        double p3y;
        double p3z;
    };

    static H5::CompType get_triangle_comp_type()
    {
        H5::CompType h5_triangle_comp_type(sizeof(h5_triangle_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(h5_triangle_comp_type.getId(), #member,\
                HOFFSET(h5_triangle_struct, member), type.getId())
        INSERT_MEMBER(p1x, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p1y, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p1z, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p2x, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p2y, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p2z, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p3x, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p3y, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(p3z, H5::PredType::NATIVE_DOUBLE);
#undef INSERT_MEMBER
        return h5_triangle_comp_type;
    }
};

void save_triangles_polygon(const Polygon&, H5::Group*);
void load_triangles_polygon(const H5::Group&, Polygon*);

} // ecell4
#endif// ECELL4_POLYGON_HDF5_WRITER_HPP
