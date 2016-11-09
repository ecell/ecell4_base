#ifndef GFRD_POLYGON_STL_FILE_READER
#define GFRD_POLYGON_STL_FILE_READER
#include <boost/array.hpp>
#include <stdexcept>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "Real3Type.hpp"

template<typename coordT>
struct StlTriangle
{
    StlTriangle(){}
    StlTriangle(const coordT& n, const boost::array<coordT, 3>& vtx)
        : normal(n), vertices(vtx)
    {}
    coordT normal;
    boost::array<coordT, 3> vertices;
};

template<typename coordT>
class StlFileReader
{
  public:
    typedef StlTriangle<coordT> triangle_type;
    enum FileType
    {
        Ascii,
        Binary,
    };

    StlFileReader(){}
    ~StlFileReader(){}

    std::vector<triangle_type>
    read(const std::string& filename, const FileType t) const;

  private:

    struct endsolid_exception{};

    std::vector<triangle_type> read_ascii(const std::string& filename) const;
    coordT        read_ascii_vertex(const std::string& line) const;
    coordT        read_ascii_normal(const std::string& line) const;
    triangle_type read_ascii_triangle(std::ifstream& ifs) const;

    std::vector<triangle_type> read_binary(const std::string& filename) const;
    coordT        read_binary_vector(std::ifstream& ifs) const;
    triangle_type read_binary_triangle(std::ifstream& ifs) const;
};

template<typename coordT>
std::vector<StlTriangle<coordT> >
StlFileReader<coordT>::read(
    const std::string& filename, const StlFileReader<coordT>::FileType type) const
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

template<typename coordT>
std::vector<StlTriangle<coordT> >
StlFileReader<coordT>::read_ascii(const std::string& filename) const
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

    std::vector<StlTriangle<coordT> > retval;
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

template<typename coordT>
StlTriangle<coordT>
StlFileReader<coordT>::read_ascii_triangle(std::ifstream& ifs) const
{
    StlTriangle<coordT> retval;
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


template<typename coordT>
coordT StlFileReader<coordT>::read_ascii_vertex(const std::string& line) const
{
    typedef typename element_type_of<coordT>::type valueT;
    std::istringstream iss(line);
    std::string prefix;
    iss >> prefix;
    if(prefix != "vertex") throw std::invalid_argument("not vertex line");
    valueT x, y, z;
    iss >> x >> y >> z;
    return coordT(x, y, z);
}

template<typename coordT>
coordT StlFileReader<coordT>::read_ascii_normal(const std::string& line) const
{
    typedef typename element_type_of<coordT>::type valueT;
    std::istringstream iss(line);
    std::string facet, normal;
    iss >> facet >> normal;
    if(facet != "facet" || normal != "normal")
        throw std::invalid_argument("not vertex line");
    valueT x, y, z;
    iss >> x >> y >> z;
    return coordT(x, y, z);
}

template<typename coordT>
std::vector<StlTriangle<coordT> >
StlFileReader<coordT>::read_binary(const std::string& filename) const
{
    std::ifstream ifs(filename.c_str(), std::ios::in | std::ios::binary);
    if(!ifs.good())
        throw std::runtime_error("file open error");

    ifs.seekg(0, ifs.end);
    const std::size_t size_of_file = ifs.tellg();
    ifs.seekg(0, ifs.beg);

    char *ch_header = new char[81];
    ifs.read(ch_header, 80);
    ch_header[80] = '\0';
    const std::string header(ch_header);
    delete [] ch_header;
    std::cerr << "header   : " << header << std::endl;

    char *ch_numTriangle = new char [sizeof(unsigned int)];
    ifs.read(ch_numTriangle, sizeof(unsigned int));
    const std::size_t num_Triangle = *reinterpret_cast<unsigned int*>(ch_numTriangle);
    delete [] ch_numTriangle;
    std::cerr << "# of face: " << num_Triangle << std::endl;

    if(50 * num_Triangle + 80 + sizeof(unsigned int) != size_of_file)
    {
        std::cerr << "file size must be 50 * number of triangle + 80 + 4" << std::endl;
        std::cerr << " = " << 50 * num_Triangle + 80 + sizeof(unsigned int) << std::endl;
        throw std::runtime_error("invalid filesize");
    }

    std::vector<StlTriangle<coordT> > retval(num_Triangle);
    for(std::size_t i=0; i < num_Triangle; ++i)
    {
        retval.at(i) = this->read_binary_triangle(ifs);
    }
    return retval;
}

template<typename coordT>
coordT StlFileReader<coordT>::read_binary_vector(std::ifstream& ifs) const
{
    char *float0 = new char [sizeof(float)];
    char *float1 = new char [sizeof(float)];
    char *float2 = new char [sizeof(float)];

    ifs.read(float0, sizeof(float));
    ifs.read(float1, sizeof(float));
    ifs.read(float2, sizeof(float));

    const float x = *reinterpret_cast<float*>(float0);
    const float y = *reinterpret_cast<float*>(float1);
    const float z = *reinterpret_cast<float*>(float2);

    delete [] float0;
    delete [] float1;
    delete [] float2;

    return coordT(x, y, z);
}

template<typename coordT>
StlTriangle<coordT>
StlFileReader<coordT>::read_binary_triangle(std::ifstream& ifs) const
{
    const coordT normal = read_binary_vector(ifs);
    boost::array<coordT, 3> vertices;
    vertices[0] = this->read_binary_vector(ifs);
    vertices[1] = this->read_binary_vector(ifs);
    vertices[2] = this->read_binary_vector(ifs);

    char *unused_data = new char [2];
    ifs.read(unused_data, 2);
    delete [] unused_data;

    return StlTriangle<coordT>(normal, vertices);
}

#endif /* GFRD_POLYGON_STL_FILE_READER */
