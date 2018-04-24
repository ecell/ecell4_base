#ifndef ECELL4_STL_FILE_READER
#define ECELL4_STL_FILE_READER
#include <ecell4/core/Triangle.hpp>
#include <string>
#include <vector>
#include <fstream>

namespace ecell4
{

struct STLFormat {enum Kind {Ascii, Binary};};

std::vector<Triangle>
read_stl_format(const std::string& filename, const STLFormat::Kind);

void
write_stl_format(const std::string& filename,
                 const std::vector<Triangle>& triangles, const STLFormat::Kind);

}// ecell4
#endif /* ECELL4_STL_FILE_READER */
