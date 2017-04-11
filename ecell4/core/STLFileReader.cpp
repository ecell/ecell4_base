#include "STLFileReader.hpp"
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace ecell4
{

std::vector<StlTriangle>
StlFileReader::read(const std::string& filename,
                    const StlFileReader::FileType type) const
{
    switch(type)
    {
        case Ascii:
            return this->read_ascii(filename);
        case Binary:
            return this->read_binary(filename);
        default:
            throw std::invalid_argument("stl unknown type");
    }
}

std::vector<StlTriangle>
StlFileReader::read_ascii(const std::string& filename) const
{
    std::ifstream ifs(filename.c_str());
    if(!ifs.good())
        throw std::runtime_error("file open error");

    while(!ifs.eof())
    {
        std::string line;
        std::getline(ifs, line);
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;
        if(prefix == "solid")
        {
            std::cerr << "found solid." << std::endl;
            std::cerr << line << std::endl;
            break;
        }
    }
    if(ifs.eof())
        throw std::runtime_error("could not find solid line");

    std::vector<StlTriangle> retval;
    while(!ifs.eof())
    {
        try
        {
            retval.push_back(this->read_ascii_triangle(ifs));
        }
        catch(endsolid_exception& esl)
        {
            break;
        }
    }
    return retval;
}

StlTriangle
StlFileReader::read_ascii_triangle(std::ifstream& ifs) const
{
    StlTriangle retval;
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
                throw std::runtime_error("invalid syntax");
            normal_read = true;
            retval.normal = this->read_ascii_normal(line);
        }
        else if(prefix == "outer")
        {
            ; // outer loop
        }
        else if(prefix == "vertex")
        {
            if(vertex_index > 2)
                throw std::runtime_error("invalid syntax");
            retval.vertices.at(vertex_index) = this->read_ascii_vertex(line);
            ++vertex_index;
        }
        else if(prefix == "endloop") 
        {
            ;
        }
        else if(prefix == "endfacet")
        {
            return retval;
        }
        else if(prefix == "endsolid")
        {
            throw endsolid_exception();
        }
        else
        {
            continue; // comment line?
        }
    }
    throw std::runtime_error("invalid syntax");
}


Real3 StlFileReader::read_ascii_vertex(const std::string& line) const
{
    std::istringstream iss(line);
    std::string prefix;
    iss >> prefix;
    if(prefix != "vertex") throw std::invalid_argument("not vertex line");
    Real x, y, z;
    iss >> x >> y >> z;
    return Real3(x, y, z);
}

Real3 StlFileReader::read_ascii_normal(const std::string& line) const
{
    std::istringstream iss(line);
    std::string facet, normal;
    iss >> facet >> normal;
    if(facet != "facet" || normal != "normal")
        throw std::invalid_argument("not vertex line");
    Real x, y, z;
    iss >> x >> y >> z;
    return Real3(x, y, z);
}

std::vector<StlTriangle>
StlFileReader::read_binary(const std::string& filename) const
{
    std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
    if(not ifs.good())
        throw std::runtime_error("file open error");

    ifs.seekg(0, ifs.end);
    const std::size_t size_of_file = ifs.tellg();
    ifs.seekg(0, ifs.beg);

    char ch_header[81];
    ifs.read(ch_header, 80);
    ch_header[80] = '\0';
    const std::string header(ch_header);
    std::cerr << "header   : " << header << std::endl;

    char ch_numTriangle[sizeof(unsigned int)];
    ifs.read(ch_numTriangle, sizeof(unsigned int));
    const std::size_t num_Triangle = *reinterpret_cast<unsigned int*>(ch_numTriangle);
    std::cerr << "# of face: " << num_Triangle << std::endl;

    if(50 * num_Triangle + 80 + sizeof(unsigned int) != size_of_file)
    {
        std::cerr << "file size must be 50 * number of triangle + 80 + 4" << std::endl;
        std::cerr << " = " << 50 * num_Triangle + 80 + sizeof(unsigned int) << std::endl;
        throw std::runtime_error("invalid filesize");
    }

    std::vector<StlTriangle> retval(num_Triangle);
    for(std::size_t i=0; i < num_Triangle; ++i)
    {
        retval.at(i) = this->read_binary_triangle(ifs);
    }
    return retval;
}

Real3 StlFileReader::read_binary_vector(std::ifstream& ifs) const
{
    char float0[sizeof(float)];
    char float1[sizeof(float)];
    char float2[sizeof(float)];

    ifs.read(float0, sizeof(float));
    ifs.read(float1, sizeof(float));
    ifs.read(float2, sizeof(float));

    const float x = *reinterpret_cast<float*>(float0);
    const float y = *reinterpret_cast<float*>(float1);
    const float z = *reinterpret_cast<float*>(float2);

    return Real3(x, y, z);
}

StlTriangle
StlFileReader::read_binary_triangle(std::ifstream& ifs) const
{
    const Real3 normal = read_binary_vector(ifs);
    boost::array<Real3, 3> vertices;
    vertices[0] = this->read_binary_vector(ifs);
    vertices[1] = this->read_binary_vector(ifs);
    vertices[2] = this->read_binary_vector(ifs);
    ifs.ignore(2);

    return StlTriangle(normal, vertices);
}

void StlFileReader::dump(const std::string& filename,
        const std::vector<triangle_type>& tri) const
{
    std::ofstream ofs(filename.c_str());
    if(!ofs.good())
    {
        std::cerr << "file " << filename << " open error" << std::endl;
        return;
    }

    ofs << "solid dumped" << std::endl;
    ofs << std::setprecision(16);
    for(typename std::vector<triangle_type>::const_iterator
        iter = tri.begin(); iter != tri.end(); ++iter)
    {
        ofs << "facet normal " << iter->normal[0] << " "
            << iter->normal[1] << " " << iter->normal[1] << std::endl;
        ofs << "outer loop" << std::endl;
        ofs << "vertex " << iter->vertices.at(0)[0] << " "
            << iter->vertices.at(0)[1] << " "
            << iter->vertices.at(0)[2] << std::endl;
        ofs << "vertex " << iter->vertices.at(1)[0] << " "
            << iter->vertices.at(1)[1] << " "
            << iter->vertices.at(1)[2] << std::endl;
        ofs << "vertex " << iter->vertices.at(2)[0] << " "
            << iter->vertices.at(2)[1] << " "
            << iter->vertices.at(2)[2] << std::endl;
        ofs << "endloop" << std::endl;
        ofs << "endfacet" << std::endl;
    }
    ofs << "endsolid dumped" << std::endl;
    ofs.close();
    return;
}

}// ecell4
