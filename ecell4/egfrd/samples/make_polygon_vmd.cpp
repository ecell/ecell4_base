#include <ecell4/egfrd/StlFileReader.hpp>
#include <ecell4/core/Real3.hpp>

int main(int argc, char **argv)
{
    if(argc != 3) return 1;
    std::string filetype(argv[1]);
    std::string filename(argv[2]);

    StlFileReader<ecell4::Real3> reader;
    std::vector<StlTriangle<ecell4::Real3> > triangles;
    if(filetype == "-asc")
        triangles = reader.read(filename, StlFileReader<ecell4::Real3>::Ascii);
    else if(filetype == "-bin")
        triangles = reader.read(filename, StlFileReader<ecell4::Real3>::Binary);
    else
        return 1;

    std::cout << "draw material Transpalent" << std::endl;
    for(std::vector<StlTriangle<ecell4::Real3> >::const_iterator
        iter = triangles.begin(); iter != triangles.end(); ++iter)
    {
        std::cout << "draw triangle {"
                  << iter->vertices[0][0] << " "
                  << iter->vertices[0][1] << " "
                  << iter->vertices[0][2] << "} {"
                  << iter->vertices[1][0] << " "
                  << iter->vertices[1][1] << " "
                  << iter->vertices[1][2] << "} {"
                  << iter->vertices[2][0] << " "
                  << iter->vertices[2][1] << " "
                  << iter->vertices[2][2] << "}" << std::endl;
    }

    return 0;
}
