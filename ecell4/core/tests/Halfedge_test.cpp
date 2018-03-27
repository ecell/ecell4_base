#define BOOST_TEST_MODULE "Polygon_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/serialization/strong_typedef.hpp>
#include <boost/algorithm/cxx11/is_permutation.hpp>
#include <boost/assign.hpp>
#include <ecell4/core/HalfEdgeMesh.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Triangle.hpp>
#include <utility>

using ecell4::Real;
using ecell4::Real3;
using ecell4::Triangle;
typedef ecell4::HalfEdgePolygon Polygon;

//! test data 1: tetrahedron
// below, the normal vector towords the depth of your display.
//
//          _4
//    3__--- /
//   /|\    /
//  / | \  /
// /__|__\/
//4  1|  /2
//    | /
//    |/
//    4
//
// p1 = {0, 0, 0}
// p2 = {1, 0, 0}
// p3 = {0, 1, 0}
// p4 = {0, 0, 1}
struct tetrahedron
{
    const static Real3 p1;
    const static Real3 p2;
    const static Real3 p3;
    const static Real3 p4;

    static Polygon make()
    {
        const Triangle t1(p1, p2, p4);
        const Triangle t2(p1, p4, p3);
        const Triangle t3(p1, p3, p2);
        const Triangle t4(p2, p3, p4);

        std::vector<Triangle> triangles;
        triangles.push_back(t1);
        triangles.push_back(t2);
        triangles.push_back(t3);
        triangles.push_back(t4);

        return Polygon(Real3(10.0, 10.0, 10.0), triangles);
    }
};
const Real3 tetrahedron::p1 = Real3(0, 0, 0);
const Real3 tetrahedron::p2 = Real3(1, 0, 0);
const Real3 tetrahedron::p3 = Real3(0, 1, 0);
const Real3 tetrahedron::p4 = Real3(0, 0, 1);

BOOST_AUTO_TEST_CASE(Polygon_tetrahedron_construction_from_triangles)
{
    const Real pi = boost::math::constants::pi<Real>();
    const Polygon poly = tetrahedron::make();

    // check shape detection
    BOOST_CHECK_EQUAL(poly.face_size(), 4);
    BOOST_CHECK_EQUAL(poly.edge_size(), 6 * 2);
    BOOST_CHECK_EQUAL(poly.vertex_size(), 4);
    BOOST_CHECK_CLOSE(poly.total_area(), 1.5 + std::sqrt(3) / 2, 1e-8);

    // check vertex position
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(tetrahedron::p1)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(tetrahedron::p2)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(tetrahedron::p3)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(tetrahedron::p4)));

    typedef Polygon::vertex_id_type vertex_id_type;
    const vertex_id_type v1 = *poly.find_vertex(tetrahedron::p1);
    const vertex_id_type v2 = *poly.find_vertex(tetrahedron::p2);
    const vertex_id_type v3 = *poly.find_vertex(tetrahedron::p3);
    const vertex_id_type v4 = *poly.find_vertex(tetrahedron::p4);

    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v1), tetrahedron::p1) -
            tetrahedron::p1), 1e-8);
    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v2), tetrahedron::p2) -
            tetrahedron::p2), 1e-8);
    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v3), tetrahedron::p3) -
            tetrahedron::p3), 1e-8);
    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v4), tetrahedron::p4) -
            tetrahedron::p4), 1e-8);


    // check apex angle
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v1), 1.5     * pi, 1e-8);
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v2), 5.0/6.0 * pi, 1e-8);
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v3), 5.0/6.0 * pi, 1e-8);
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v4), 5.0/6.0 * pi, 1e-8);

    // check all the outgoing_edges are terminated at the correct vertex.
    {
        const std::vector<vertex_id_type> ans = boost::assign::list_of(v2)(v3)(v4);

        std::vector<vertex_id_type> result; result.reserve(3);
        for(typename std::vector<vertex_id_type>::const_iterator
            i(poly.outgoing_edges(v1).begin()), e(poly.outgoing_edges(v1).end());
            i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(boost::algorithm::is_permutation(
                    ans.begin(), ans.end(), result.begin()));
    }
    {
        const std::vector<vertex_id_type> ans = boost::assign::list_of(v1)(v3)(v4);

        std::vector<vertex_id_type> result; result.reserve(3);
        for(typename std::vector<vertex_id_type>::const_iterator
            i(poly.outgoing_edges(v2).begin()), e(poly.outgoing_edges(v2).end());
            i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(boost::algorithm::is_permutation(
                    ans.begin(), ans.end(), result.begin()));
    }
    {
        const std::vector<vertex_id_type> ans = boost::assign::list_of(v1)(v2)(v4);

        std::vector<vertex_id_type> result; result.reserve(3);
        for(typename std::vector<vertex_id_type>::const_iterator
            i(poly.outgoing_edges(v3).begin()), e(poly.outgoing_edges(v3).end());
            i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(boost::algorithm::is_permutation(
                    ans.begin(), ans.end(), result.begin()));
    }
    {
        const std::vector<vertex_id_type> ans = boost::assign::list_of(v1)(v2)(v3);

        std::vector<vertex_id_type> result; result.reserve(3);
        for(typename std::vector<vertex_id_type>::const_iterator
            i(poly.outgoing_edges(v4).begin()), e(poly.outgoing_edges(v4).end());
            i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(boost::algorithm::is_permutation(
                    ans.begin(), ans.end(), result.begin()));
    }

}

//! test data 2: octahedron
// below, the normal vector towords the depth of your display.
//
//       3
//       /\
// 3___1/__\4__3
//  \  /\  /\  /\
//   \/__\/__\/__\
//   2   5\  /6   2
//         \/
//          2
//
// p1 = {1, 1, 2}
// p2 = {2, 1, 1}
// p3 = {1, 2, 1}
// p4 = {0, 1, 1}
// p5 = {1, 0, 1}
// p6 = {1, 1, 0}
struct octahedron
{
    const static Real3 p1;
    const static Real3 p2;
    const static Real3 p3;
    const static Real3 p4;
    const static Real3 p5;
    const static Real3 p6;

    static Polygon make()
    {
        const Triangle t1(p1, p2, p3);
        const Triangle t3(p1, p4, p5);
        const Triangle t4(p1, p3, p4);
        const Triangle t2(p1, p5, p2);

        const Triangle t5(p6, p5, p4);
        const Triangle t6(p6, p4, p3);
        const Triangle t7(p6, p3, p2);
        const Triangle t8(p6, p2, p5);

        std::vector<Triangle> triangles;
        triangles.push_back(t1);
        triangles.push_back(t2);
        triangles.push_back(t3);
        triangles.push_back(t4);

        triangles.push_back(t5);
        triangles.push_back(t6);
        triangles.push_back(t7);
        triangles.push_back(t8);

        return Polygon(Real3(10.0, 10.0, 10.0), triangles);
    }
};
const Real3 octahedron::p1 = Real3(1, 1, 2);
const Real3 octahedron::p2 = Real3(2, 1, 1);
const Real3 octahedron::p3 = Real3(1, 2, 1);
const Real3 octahedron::p4 = Real3(0, 1, 1);
const Real3 octahedron::p5 = Real3(1, 0, 1);
const Real3 octahedron::p6 = Real3(1, 1, 0);




