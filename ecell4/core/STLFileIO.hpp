#ifndef ECELL4_STL_FILE_READER
#define ECELL4_STL_FILE_READER
#include <ecell4/core/Triangle.hpp>
#include <string>
#include <vector>
#include <cstdint>

namespace ecell4
{

class Polygon; // forward declaration

enum class STLFormat : std::uint8_t
{
    Ascii  = 0,
    Binary = 1,
};

std::vector<Triangle> read_stl_format(const std::string& fname, const STLFormat);
void                 write_stl_format(const std::string& fname, const STLFormat,
                                      const std::vector<Triangle>& triangles);

Polygon read_polygon(const std::string& filename, const STLFormat,
                     const Real3& edge_lengths);
void   write_polygon(const std::string& filename, const STLFormat,
                     const Polygon&);

}// ecell4
#endif /* ECELL4_STL_FILE_READER */
