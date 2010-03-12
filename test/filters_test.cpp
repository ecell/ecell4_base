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
    typedef typename Toc_::position_type::value_type distance_type;
    collector(): result() {}

    void operator()(typename Toc_::iterator i, distance_type dist)
    {
        result.insert(std::make_pair((*i).first, dist));
    }
    std::map<typename Toc_::key_type, distance_type> result;
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

    // Collect all spheres from oc who (partly) lie within the radius of 
    // sphere 1.
    collector<oc_type> col;
    oc_type::const_iterator f(oc.find(1));
    take_neighbor(oc, col, (*f).second);

    BOOST_CHECK_EQUAL(col.result.size(), 2);
    // Sphere 0 overlaps with sphere 1.
    BOOST_CHECK(std::fabs(std::sqrt(0.1 * 0.1 * 2) - 0.15 - col.result[0]) <
                1e-8);
    // Distance to shell of sphere 1.
    BOOST_CHECK(std::fabs(-0.05 - col.result[1]) < 1e-8);

    BOOST_CHECK(col.result.end() == col.result.find(2));
    BOOST_CHECK(col.result.end() == col.result.find(3));
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
    BOOST_CHECK(std::fabs(0.05 - col.result[0]) < 1e-8);
    // pos lies within cylinder 1.
    BOOST_CHECK(std::fabs(-0.05 - col.result[1]) < 1e-8);

    BOOST_CHECK(col.result.end() == col.result.find(2));
    BOOST_CHECK(col.result.end() == col.result.find(3));
}
