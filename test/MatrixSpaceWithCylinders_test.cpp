#define BOOST_TEST_MODULE "MatrixSpaceWithCylinders_test"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <functional>
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "Cylinder.hpp"
#include "MatrixSpace.hpp"

BOOST_AUTO_TEST_CASE(insert)
{
    typedef double length_type;
    typedef int key_type;
    typedef Cylinder<length_type> mapped_type;
    typedef MatrixSpace<mapped_type, key_type> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1000.0, 10);
    BOOST_CHECK_CLOSE(100., oc.cell_size(), 0.001);
    {
        std::pair<oc_type::iterator, bool> ir(
                oc.update(std::make_pair(
                    0, oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000))));
        BOOST_CHECK_EQUAL(true, ir.second);  // Normal insert.
        BOOST_CHECK(oc.end() != oc.find(0)); // Key 0 exists.
        BOOST_CHECK(oc.end() == oc.find(1)); // Key 1 doesn't exist.
    }
    {
        // Update.
        std::pair<oc_type::iterator, bool> ir(
                oc.update(std::make_pair(
                    0, oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000))));
        BOOST_CHECK_EQUAL(false, ir.second); // False: this was an update.
        // ir.first is an iterator to the value you inserted. So accessing 
        // it's second element should return you the object.
        BOOST_CHECK_EQUAL(oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000),
                (*ir.first).second);
        BOOST_CHECK(oc.end() != oc.find(0));
        BOOST_CHECK(oc.end() == oc.find(1));
    }
    {
        // Another update.
        std::pair<oc_type::iterator, bool> ir(
                oc.update(std::make_pair(
                    0, oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000))));
        BOOST_CHECK_EQUAL(false, ir.second);
        BOOST_CHECK_EQUAL(oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000),
                (*ir.first).second);
        BOOST_CHECK(oc.end() != oc.find(0));
        BOOST_CHECK(oc.end() == oc.find(1));
    }
}

template<typename Toc_>
struct collector
{
    void operator()(typename Toc_::iterator i)
    {
        result.insert((*i).first);
    }

    std::set<typename Toc_::key_type> result;
};

template<typename Toc_>
struct collector2
{
    void operator()(typename Toc_::iterator i,
            const typename Toc_::position_type& pos_off)
    {
        result.insert((*i).first);
    }
    std::set<typename Toc_::key_type> result;
};

BOOST_AUTO_TEST_CASE(each_neighbor)
{
    typedef double length_type;
    typedef int key_type;
    typedef Cylinder<length_type> mapped_type;
    typedef MatrixSpace<mapped_type, key_type> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1000., 10);
    BOOST_CHECK_CLOSE(100., oc.cell_size(), 0.001);
 
    // Insert value 0.
    oc.update(std::make_pair(0, oc_type::mapped_type(pos(500, 500, 0), 25, pos(0,0,1), 50)));
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() == oc.find(1));

    {
        collector<oc_type> col;
        // Should return value 0.
        oc.each_neighbor(oc.index(pos(500, 500, 100)), col);
        BOOST_CHECK_EQUAL(col.result.size(), 1);
        BOOST_CHECK(col.result.find(0) != col.result.end());
    }

    {
        collector<oc_type> col;
        // No periodic boundary condition. Should return no values.
        // Behaviour is unspecified for values at the boundary or out of the 
        // MatrixSpace (x,y,z >= 1000).
        oc.each_neighbor(oc.index(pos(500, 500, 900)), col);
        BOOST_CHECK_EQUAL(col.result.size(), 0);
    }

    {
        collector2<oc_type> col2;
        // Periodic boundary condition. Should return element 0 after applying 
        // periodic boundary condition in z (add 1000 to z coordinate of the 
        // origin of the cylinder to be in the same neighbourhood as reference 
        // point), so: (0,0,1000).
        oc.each_neighbor_cyclic(oc.index(pos(500, 500, 900)), col2);
        BOOST_CHECK_EQUAL(col2.result.size(), 1);
        BOOST_CHECK(col2.result.find(0) != col2.result.end());
    }

    // Insert value 1.
    oc.update(std::make_pair(1, oc_type::mapped_type(pos(500, 500, 900), 25, pos(0,0,1), 50)));
    {
        BOOST_CHECK(oc.end() != oc.find(0));
        BOOST_CHECK(oc.end() != oc.find(1));
        BOOST_CHECK(oc.end() == oc.find(2));
    }

    {
        collector2<oc_type> col2;
        // Periodic boundary condition. Should return element 0 (0, 0, 0) and 
        // element 1 (0,0,-1000).
        oc.each_neighbor_cyclic(oc.index(pos(500, 500, 0)), col2);
        BOOST_CHECK_EQUAL(col2.result.size(), 2);
        BOOST_CHECK(col2.result.find(0) != col2.result.end());
        BOOST_CHECK(col2.result.find(1) != col2.result.end());
    }
}

