#define BOOST_TEST_MODULE "filters_test"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/test/included/unit_test.hpp>
#include "MatrixSpace.hpp"
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "filters.hpp"

template<typename Toc_>
struct collector
{
    void operator()(typename Toc_::iterator i,
            typename Toc_::position_type::value_type dist)
    {
        std::cout << (*i).second << ", " << dist << std::endl;
    }
};


BOOST_AUTO_TEST_CASE(spheres)
{
    typedef MatrixSpace<Sphere<double>, int> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1.0, 10);

    oc.update(std::make_pair(0, oc_type::mapped_type(pos(0.2, 0.6, 0.4), 0.15)));
    oc.update(std::make_pair(1, oc_type::mapped_type(pos(0.2, 0.7, 0.5), 0.05)));
    oc.update(std::make_pair(2, oc_type::mapped_type(pos(0.9, 0.1, 0.4), 0.07)));
    oc.update(std::make_pair(3, oc_type::mapped_type(pos(0.9, 0.95, 0.4), 0.1)));

    collector<oc_type> col;
    oc_type::const_iterator f(oc.find(1));
    take_neighbor(oc, col, (*f).second);
    // Output should be.
    //{(0.2, 0.6, 0.4), 0.15}, -0.00857864
    //{(0.2, 0.7, 0.5), 0.05}, -0.05

}

BOOST_AUTO_TEST_CASE(cylinders)
{
    typedef double length_type;
    typedef int key_type;
    typedef Cylinder<length_type> mapped_type;
    typedef Sphere<length_type> sphere_type;
    typedef MatrixSpace<mapped_type, key_type> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1.0, 10);

    oc.update(std::make_pair(0, mapped_type(pos(0.2, 0.7, 0.7), 0.15, pos(0,0,1), 0.15)));
    oc.update(std::make_pair(1, mapped_type(pos(0.2, 0.7, 0.5), 0.50, pos(0,0,1), 0.05)));
    oc.update(std::make_pair(2, mapped_type(pos(0.9, 0.1, 0.4), 0.07, pos(0,0,1), 0.07)));
    oc.update(std::make_pair(3, mapped_type(pos(0.9, 0.95, 0.4), 0.1, pos(0,0,1), 0.1)));

    std::cout << std::endl << std::endl;
    typedef collector<oc_type> col_type;
    col_type col;

    // Collect all cylinders from oc who (partly) lie within the radius of 
    // sphere s.
    sphere_type s = sphere_type(pos(0.2, 0.7, 0.5), 0.05);
    take_neighbor(oc, col, s);
    // Distance to cylinder 0 should be 0.05.
    // Note that the distance to cylinder 1 is not correct because the point 
    // lies inside the cylinder. It should compute the distance to the caps, 
    // but it computes the distance to the radial edge.

    // Output should be:
    //{(0.2, 0.7, 0.5), 0.5, (0, 0, 1), 0.05}, -0.05
    //{(0.2, 0.7, 0.7), 0.15, (0, 0, 1), 0.15}, 0.05
}
