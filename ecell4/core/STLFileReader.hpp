#ifndef ECELL4_STL_FILE_READER
#define ECELL4_STL_FILE_READER
#include <ecell4/core/Real3.hpp>
#include <string>
#include <vector>
#include <fstream>

namespace ecell4
{

struct StlTriangle
{
    StlTriangle(){}
    StlTriangle(const Real3& n, const boost::array<Real3, 3>& vtx)
        : normal(n), vertices(vtx)
    {}
    Real3 normal;
    boost::array<Real3, 3> vertices;
};

class StlFileReader
{
  public:
    typedef StlTriangle triangle_type;

    enum FileType
    {
        Ascii,
        Binary,
    };

  public:

    StlFileReader(){}
    ~StlFileReader(){}

    std::vector<triangle_type>
    read(const std::string& filename, const FileType t) const;

    void
    dump(const std::string& filename, const std::vector<triangle_type>& tri) const;

  private:

    struct endsolid_exception{};

    std::vector<triangle_type> read_ascii(const std::string& filename) const;
    Real3         read_ascii_vertex(const std::string& line) const;
    Real3         read_ascii_normal(const std::string& line) const;
    triangle_type read_ascii_triangle(std::ifstream& ifs) const;

    std::vector<triangle_type> read_binary(const std::string& filename) const;
    Real3         read_binary_vector(std::ifstream& ifs) const;
    triangle_type read_binary_triangle(std::ifstream& ifs) const;
};


}// ecell4
#endif /* ECELL4_STL_FILE_READER */
