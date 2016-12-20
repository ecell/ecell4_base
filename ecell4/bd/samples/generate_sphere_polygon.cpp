#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <ecell4/core/Real3.hpp>

using ecell4::Real3;

void write_stl(const Real3& a, const Real3& b, const Real3& c)
{
    const Real3 normal_ = cross_product(b-a, c-b);
    const Real3 normal = normal_ / length(normal_);
    std::cout << "facet normal " << normal << std::endl;
    std::cout << "outer loop" << std::endl;
    std::cout << "vertex " << a << std::endl;
    std::cout << "vertex " << b << std::endl;
    std::cout << "vertex " << c << std::endl;
    std::cout << "endloop" << std::endl;
    std::cout << "endfacet" << std::endl;
    return ;
}

int main(int argc, char **argv)
{
    if(argc != 3)
    {
        std::cerr << "Usage: ./" << argv[0]
                  << " <resolution> <radius>" << std::endl;;
        return 1;
    }

    try
    {
    const std::size_t N = std::atoi(argv[1]);
    const double radius = std::atof(argv[2]);
    std::cout << "solid sphere_radius_" << radius << std::endl;

    // +z hemisphere
    for(std::size_t q = 0; q < 4; ++q)
    {
        const double buf = q * M_PI / 2;
        std::vector<Real3> vertices;
        vertices.push_back(Real3(0., 0., radius));
        for(std::size_t i=1; i<=N; ++i)
        {
            const double phi  = i * M_PI / (2 * N);
            const double z    = radius * std::cos(phi);
            const double lrad = radius * std::sin(phi);

            for(std::size_t j=0; j<=i; ++j)
            {
                const double local_angle = j * M_PI / (2 * i) + buf;
                const double x = lrad * std::cos(local_angle);
                const double y = lrad * std::sin(local_angle);
                vertices.push_back(Real3(x, y, z));
            }
        }

        for(std::size_t i=1; i<N; ++i)
        {
            const std::size_t start =  i    * (i+1) / 2;
            const std::size_t goal  = (i+1) * (i+2) / 2;
            for(std::size_t j = start; j < goal-1; ++j)
            {
                write_stl(vertices.at(j), vertices.at(j+1), vertices.at(j-i));
                write_stl(vertices.at(j), vertices.at(j+i+2), vertices.at(j+1));
            }
        }
        // edge
        for(std::size_t i=N*(N+1)/2; i< (N+1)*(N+2)/2-1; ++i)
        {
            write_stl(vertices.at(i), vertices.at(i+1), vertices.at(i-N));
        }

    }

    // -z hemisphere
    for(std::size_t q = 0; q < 4; ++q)
    {
        const double buf = q * M_PI / 2;
        std::vector<Real3> vertices;
        vertices.push_back(Real3(0., 0., -radius));
        for(std::size_t i=1; i<=N; ++i)
        {
            const double phi  = i * M_PI / (2 * N);
            const double z    = -radius * std::cos(phi);
            const double lrad =  radius * std::sin(phi);

            for(std::size_t j=0; j<=i; ++j)
            {
                const double local_angle = j * M_PI / (2 * i) + buf;
                const double x = lrad * std::cos(-local_angle);
                const double y = lrad * std::sin(-local_angle);
                vertices.push_back(Real3(x, y, z));
            }
        }
        for(std::size_t i=1; i<N; ++i)
        {
            const std::size_t start =  i    * (i+1) / 2;
            const std::size_t goal  = (i+1) * (i+2) / 2;
            for(std::size_t j = start; j < goal-1; ++j)
            {
                write_stl(vertices.at(j), vertices.at(j+1), vertices.at(j-i));
                write_stl(vertices.at(j), vertices.at(j+i+2), vertices.at(j+1));
            }
        }
        // edge
        for(std::size_t i=N*(N+1)/2; i< (N+1)*(N+2)/2-1; ++i)
        {
            write_stl(vertices.at(i), vertices.at(i+1), vertices.at(i-N));
        }
    }
    }
    catch(std::exception& except)
    {
        std::cerr << "error occured. what() = " << except.what() << std::endl;
        return 1;
    }

    return 0;
}
