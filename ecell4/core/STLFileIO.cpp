#include <ecell4/core/STLFileIO.hpp>
#include <boost/cstdint.hpp>
#include <boost/format.hpp>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace ecell4
{

struct endsolid_appeared{};

static Real3 read_ascii_stl_vertex(const std::string& line)
{
    std::istringstream iss(line);

    std::string prefix;
    iss >> prefix;
    if(prefix != "vertex")
    {
        throw std::runtime_error("syntax error: missing vertex line");
    }

    Real x, y, z;
    iss >> x >> y >> z;
    return Real3(x, y, z);
}

static Real3 read_ascii_stl_normal(const std::string& line)
{
    std::istringstream iss(line);

    std::string facet, normal;
    iss >> facet >> normal;
    if(facet != "facet" || normal != "normal")
    {
        throw std::runtime_error("syntax error: missing `facet normal`");
    }

    Real x, y, z;
    iss >> x >> y >> z;
    return Real3(x, y, z);
}

static Triangle read_ascii_stl_triangle(std::ifstream& ifs)
{
    boost::array<Real3, 3> vs;
    bool normal_read = false;
    std::size_t vertex_index = 0;
    while(!ifs.eof())
    {
        std::string line;
        std::getline(ifs, line);
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if(prefix == "facet")
        {
            if(normal_read)
            {
                throw std::runtime_error("syntax error: duplicated `normal`");
            }
            normal_read = true;
            // XXX ignore normal written in the file
            const volatile Real3 normal = read_ascii_stl_normal(line);
        }
        else if(prefix == "outer")
        {
            ; // outer loop. ignore.
        }
        else if(prefix == "vertex")
        {
            if(vertex_index > 2)
            {
                throw NotSupported("STL contains more than 3 vertices");
            }
            vs.at(vertex_index) = read_ascii_stl_vertex(line);
            ++vertex_index;
        }
        else if(prefix == "endloop")
        {
            ;
        }
        else if(prefix == "endfacet")
        {
            return Triangle(vs);
        }
        else if(prefix == "endsolid")
        {
            throw endsolid_appeared();
        }
        else
        {
            // comment line? do nothing.
        }
        ifs.peek();
    }
    throw std::runtime_error("invalid syntax");
}

static std::vector<Triangle> read_ascii_stl(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());
    if(!ifs.good())
    {
        throw std::runtime_error("file open error: " + filename);
    }

    while(!ifs.eof())
    {
        std::string line;
        std::getline(ifs, line);
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        if(prefix == "solid")
        {
//             std::cerr << "found solid." << std::endl;
//             std::cerr << line << std::endl;
            break;
        }
        ifs.peek();
    }
    if(ifs.eof())
    {
        throw std::runtime_error("could not find solid line");
    }

    std::vector<Triangle> retval;
    while(!ifs.eof())
    {
        try
        {
            retval.push_back(read_ascii_stl_triangle(ifs));
        }
        catch(endsolid_appeared& esl)
        {
            break;
        }
        ifs.peek();
    }
    return retval;
}

static Real3 read_binary_stl_vector(std::ifstream& ifs)
{
    float x, y, z;
    ifs.read(reinterpret_cast<char*>(&x), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&y), sizeof(float));
    ifs.read(reinterpret_cast<char*>(&z), sizeof(float));
    return Real3(x, y, z);
}

static Triangle read_binary_stl_triangle(std::ifstream& ifs)
{
    // ignore normal vector written in the file
    const volatile Real3 normal = read_binary_stl_vector(ifs);
    boost::array<Real3, 3> vs;
    vs[0] = read_binary_stl_vector(ifs);
    vs[1] = read_binary_stl_vector(ifs);
    vs[2] = read_binary_stl_vector(ifs);
    ifs.ignore(2);
    return Triangle(vs);
}

static std::vector<Triangle>
read_binary_stl(const std::string& filename)
{
    std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
    if(!ifs.good())
    {
        throw std::runtime_error("file open error: " + filename);
    }

    ifs.seekg(0, ifs.end);
    const std::size_t size_of_file = ifs.tellg();
    ifs.seekg(0, ifs.beg);

    char ch_header[81];
    ifs.read(ch_header, 80);
    ch_header[80] = '\0';
//     std::cerr << "header   : " << ch_header << std::endl;

    boost::uint32_t num_triangle = 0;
    ifs.read(reinterpret_cast<char*>(&num_triangle), 4);
//     std::cerr << "# of face: " << num_triangle << std::endl;

    if(50 * num_triangle + 84 != size_of_file)
    {
        throw std::runtime_error((boost::format("ecell4::read_binary_stl: "
            "invalid filesize: %1% != %2% triagnles * 50 + header(84)") %
            size_of_file % num_triangle).str());
    }

    std::vector<Triangle> retval(num_triangle);
    for(boost::uint32_t i=0; i < num_triangle; ++i)
    {
        retval.at(i) = read_binary_stl_triangle(ifs);
    }
    return retval;
}

std::vector<Triangle>
read_stl_format(const std::string& filename, const STLFormat::Kind kind)
{
    switch(kind)
    {
        case STLFormat::Ascii:  return read_ascii_stl(filename);
        case STLFormat::Binary: return read_binary_stl(filename);
        default: throw std::invalid_argument("read_stl_format: unknown format");
    }
}

static void write_binary_stl(
    const std::string& filename, const std::vector<Triangle>& tri)
{
    std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
    if(!ofs.good())
    {
        throw std::runtime_error(
            "ecell4::write_stl_format: file open error: " + filename);
    }

    const std::string header("this file is generated by ecell4::write_"
                             "stl_format.                             ");
    assert(header.size() == 80);
    ofs.write(header.c_str(), 80);

    const boost::uint32_t num_triangle = tri.size();
    ofs.write(reinterpret_cast<const char*>(&num_triangle), 4);

    for(typename std::vector<Triangle>::const_iterator
            iter = tri.begin(); iter != tri.end(); ++iter)
    {
        const float nx(iter->normal()[0]);
        const float ny(iter->normal()[1]);
        const float nz(iter->normal()[2]);
        const float v0x(iter->vertex_at(0)[0]);
        const float v0y(iter->vertex_at(0)[1]);
        const float v0z(iter->vertex_at(0)[2]);
        const float v1x(iter->vertex_at(1)[0]);
        const float v1y(iter->vertex_at(1)[1]);
        const float v1z(iter->vertex_at(1)[2]);
        const float v2x(iter->vertex_at(2)[0]);
        const float v2y(iter->vertex_at(2)[1]);
        const float v2z(iter->vertex_at(2)[2]);

        ofs.write(reinterpret_cast<const char*>(&nx),  4);
        ofs.write(reinterpret_cast<const char*>(&ny),  4);
        ofs.write(reinterpret_cast<const char*>(&nz),  4);
        ofs.write(reinterpret_cast<const char*>(&v0x), 4);
        ofs.write(reinterpret_cast<const char*>(&v0y), 4);
        ofs.write(reinterpret_cast<const char*>(&v0z), 4);
        ofs.write(reinterpret_cast<const char*>(&v1x), 4);
        ofs.write(reinterpret_cast<const char*>(&v1y), 4);
        ofs.write(reinterpret_cast<const char*>(&v1z), 4);
        ofs.write(reinterpret_cast<const char*>(&v2x), 4);
        ofs.write(reinterpret_cast<const char*>(&v2y), 4);
        ofs.write(reinterpret_cast<const char*>(&v2z), 4);

        const boost::uint16_t attr(0);
        ofs.write(reinterpret_cast<const char*>(&attr), 2);
    }
    ofs.close();
    return;
}

static void write_ascii_stl(
    const std::string& filename, const std::vector<Triangle>& tri)
{
    std::ofstream ofs(filename.c_str());
    if(!ofs.good())
    {
        throw std::runtime_error(
            "ecell4::write_stl_format: file open error: " + filename);
    }

    ofs << "solid ecell4\n";
    ofs << std::setprecision(16);
    for(typename std::vector<Triangle>::const_iterator
            iter = tri.begin(); iter != tri.end(); ++iter)
    {
        const Real3  n  = iter->normal();
        const Real3& v0 = iter->vertex_at(0);
        const Real3& v1 = iter->vertex_at(1);
        const Real3& v2 = iter->vertex_at(2);

        ofs << "facet normal " << n[0] << ' ' << n[1] << ' ' <<  n[1] << '\n';
        ofs << "  outer loop\n";
        ofs << "    vertex " << v0[0] << ' ' << v0[1] << ' ' << v0[2] << '\n';
        ofs << "    vertex " << v1[0] << ' ' << v1[1] << ' ' << v1[2] << '\n';
        ofs << "    vertex " << v2[0] << ' ' << v2[1] << ' ' << v2[2] << '\n';
        ofs << "  endloop\n";
        ofs << "endfacet\n";
    }
    ofs << "endsolid ecell4\n";
    ofs.close();
    return;
}

void write_stl_format(const std::string& filename,
        const std::vector<Triangle>& tri, const STLFormat::Kind kind)
{
    switch(kind)
    {
        case STLFormat::Ascii:  return write_ascii_stl (filename, tri);
        case STLFormat::Binary: return write_binary_stl(filename, tri);
        default: throw std::invalid_argument("write_stl_format: unknown format");
    }
}

}// ecell4
