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
typedef Polygon::face_id_type face_id_type;
typedef Polygon::edge_id_type edge_id_type;
typedef Polygon::vertex_id_type vertex_id_type;

bool check_equal(const Real3& lhs, const Real3& rhs, const Real tol)
{
    return std::abs(lhs[0] - rhs[0]) < tol &&
           std::abs(lhs[1] - rhs[1]) < tol &&
           std::abs(lhs[2] - rhs[2]) < tol;
}

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

    // check all the vertex are in contact with the correct set of faces.
    {
        const std::vector<face_id_type> ans = boost::assign::list_of(v1)(v2)(v3);

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

    // check all the edges exist and has correct next-edge
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v1, v2)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v1, v3)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v1, v4)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v2, v1)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v2, v3)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v2, v4)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v3, v1)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v3, v2)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v3, v4)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v4, v1)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v4, v2)));
    BOOST_CHECK(static_cast<bool>(poly.find_edge(v4, v3)));

    const edge_id_type e12 = *poly.find_edge(v1, v2);
    const edge_id_type e13 = *poly.find_edge(v1, v3);
    const edge_id_type e14 = *poly.find_edge(v1, v4);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e12))), v1);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e13))), v1);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e14))), v1);

    const edge_id_type e21 = *poly.find_edge(v2, v1);
    const edge_id_type e23 = *poly.find_edge(v2, v3);
    const edge_id_type e24 = *poly.find_edge(v2, v4);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e21))), v2);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e23))), v2);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e24))), v2);

    const edge_id_type e31 = *poly.find_edge(v3, v1);
    const edge_id_type e32 = *poly.find_edge(v3, v2);
    const edge_id_type e34 = *poly.find_edge(v3, v4);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e31))), v3);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e32))), v3);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e34))), v3);

    const edge_id_type e41 = *poly.find_edge(v4, v1);
    const edge_id_type e42 = *poly.find_edge(v4, v2);
    const edge_id_type e43 = *poly.find_edge(v4, v3);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e41))), v4);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e42))), v4);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e43))), v4);

    // check opposite edges
    BOOST_CHECK_EQUAL(poly.opposite_of(e12), e21);
    BOOST_CHECK_EQUAL(poly.opposite_of(e13), e31);
    BOOST_CHECK_EQUAL(poly.opposite_of(e14), e41);

    BOOST_CHECK_EQUAL(poly.opposite_of(e21), e12);
    BOOST_CHECK_EQUAL(poly.opposite_of(e23), e32);
    BOOST_CHECK_EQUAL(poly.opposite_of(e24), e42);

    BOOST_CHECK_EQUAL(poly.opposite_of(e31), e13);
    BOOST_CHECK_EQUAL(poly.opposite_of(e32), e23);
    BOOST_CHECK_EQUAL(poly.opposite_of(e34), e43);

    BOOST_CHECK_EQUAL(poly.opposite_of(e41), e14);
    BOOST_CHECK_EQUAL(poly.opposite_of(e42), e24);
    BOOST_CHECK_EQUAL(poly.opposite_of(e43), e34);

    // check face ids
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e12)),  poly.face_of(e12));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e12))), poly.face_of(e12));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e13)),  poly.face_of(e13));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e13))), poly.face_of(e13));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e14)),  poly.face_of(e14));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e14))), poly.face_of(e14));

    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e21)),  poly.face_of(e21));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e21))), poly.face_of(e21));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e23)),  poly.face_of(e23));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e23))), poly.face_of(e23));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e24)),  poly.face_of(e24));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e24))), poly.face_of(e24));

    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e31)),  poly.face_of(e31));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e31))), poly.face_of(e31));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e32)),  poly.face_of(e32));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e32))), poly.face_of(e32));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e34)),  poly.face_of(e34));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e34))), poly.face_of(e34));

    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e41)),  poly.face_of(e41));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e41))), poly.face_of(e41));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e42)),  poly.face_of(e42));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e42))), poly.face_of(e42));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(             e43)),  poly.face_of(e43));
    BOOST_CHECK_EQUAL(poly.face_of(poly.next_of(poly.next_of(e43))), poly.face_of(e43));

    BOOST_CHECK(check_equal(poly.direction_of(e12), poly.periodic_transpose(
        poly.position_at(v2), poly.position_at(v1)) - poly.position_at(v1), 1e-8));
    BOOST_CHECK(check_equal(poly.direction_of(e13), poly.periodic_transpose(
        poly.position_at(v3), poly.position_at(v1)) - poly.position_at(v1), 1e-8));
    BOOST_CHECK(check_equal(poly.direction_of(e14), poly.periodic_transpose(
        poly.position_at(v4), poly.position_at(v1)) - poly.position_at(v1), 1e-8));

    BOOST_CHECK(check_equal(poly.direction_of(e21), poly.periodic_transpose(
        poly.position_at(v1), poly.position_at(v2)) - poly.position_at(v2), 1e-8));
    BOOST_CHECK(check_equal(poly.direction_of(e23), poly.periodic_transpose(
        poly.position_at(v3), poly.position_at(v2)) - poly.position_at(v2), 1e-8));
    BOOST_CHECK(check_equal(poly.direction_of(e24), poly.periodic_transpose(
        poly.position_at(v4), poly.position_at(v2)) - poly.position_at(v2), 1e-8));

    BOOST_CHECK(check_equal(poly.direction_of(e31), poly.periodic_transpose(
        poly.position_at(v1), poly.position_at(v3)) - poly.position_at(v3), 1e-8));
    BOOST_CHECK(check_equal(poly.direction_of(e32), poly.periodic_transpose(
        poly.position_at(v2), poly.position_at(v3)) - poly.position_at(v3), 1e-8));
    BOOST_CHECK(check_equal(poly.direction_of(e34), poly.periodic_transpose(
        poly.position_at(v4), poly.position_at(v3)) - poly.position_at(v3), 1e-8));

    BOOST_CHECK(check_equal(poly.direction_of(e41), poly.periodic_transpose(
        poly.position_at(v1), poly.position_at(v4)) - poly.position_at(v4), 1e-8));
    BOOST_CHECK(check_equal(poly.direction_of(e42), poly.periodic_transpose(
        poly.position_at(v2), poly.position_at(v4)) - poly.position_at(v4), 1e-8));
    BOOST_CHECK(check_equal(poly.direction_of(e43), poly.periodic_transpose(
        poly.position_at(v3), poly.position_at(v4)) - poly.position_at(v4), 1e-8));
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

BOOST_AUTO_TEST_CASE(Polygon_octahedron_construction_from_triangles)
{
    const Real pi = boost::math::constants::pi<Real>();
    const Polygon poly = octahedron::make();

    // check shape detection
    BOOST_CHECK_EQUAL(poly.face_size(), 8);
    BOOST_CHECK_EQUAL(poly.edge_size(), 8 * 3);
    BOOST_CHECK_EQUAL(poly.vertex_size(), 6);
    BOOST_CHECK_CLOSE(poly.total_area(), 8 * std::sqrt(3.0) / 2, 1e-8);

    // check vertex position
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p1)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p2)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p3)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p4)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p5)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p6)));

    const vertex_id_type v1 = *poly.find_vertex(octahedron::p1);
    const vertex_id_type v2 = *poly.find_vertex(octahedron::p2);
    const vertex_id_type v3 = *poly.find_vertex(octahedron::p3);
    const vertex_id_type v4 = *poly.find_vertex(octahedron::p4);
    const vertex_id_type v5 = *poly.find_vertex(octahedron::p5);
    const vertex_id_type v6 = *poly.find_vertex(octahedron::p6);

    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v1), octahedron::p1) -
            octahedron::p1), 1e-8);
    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v2), octahedron::p2) -
            octahedron::p2), 1e-8);
    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v3), octahedron::p3) -
            octahedron::p3), 1e-8);
    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v4), octahedron::p4) -
            octahedron::p4), 1e-8);
    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v5), octahedron::p5) -
            octahedron::p5), 1e-8);
    BOOST_CHECK_SMALL(ecell4::length(
            poly.periodic_transpose(poly.position_at(v6), octahedron::p6) -
            octahedron::p6), 1e-8);

    // check apex angle
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v1), 4.0 / 3.0 * pi, 1e-8);
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v2), 4.0 / 3.0 * pi, 1e-8);
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v3), 4.0 / 3.0 * pi, 1e-8);
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v4), 4.0 / 3.0 * pi, 1e-8);
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v5), 4.0 / 3.0 * pi, 1e-8);
    BOOST_CHECK_CLOSE(poly.apex_angle_at(v6), 4.0 / 3.0 * pi, 1e-8);

    // check all the outgoing_edges are terminated at the correct vertex.
    /* v1 */{
        const std::vector<vertex_id_type> ans =
            boost::assign::list_of(v2)(v3)(v4)(v5);

        std::vector<vertex_id_type> result; result.reserve(4);
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
    /* v2 */{
        const std::vector<vertex_id_type> ans =
            boost::assign::list_of(v1)(v3)(v5)(v6);

        std::vector<vertex_id_type> result; result.reserve(4);
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
    /* v3 */{
        const std::vector<vertex_id_type> ans =
            boost::assign::list_of(v1)(v2)(v4)(v6);

        std::vector<vertex_id_type> result; result.reserve(4);
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
    /* v4 */{
        const std::vector<vertex_id_type> ans =
            boost::assign::list_of(v1)(v3)(v5)(v6);

        std::vector<vertex_id_type> result; result.reserve(4);
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
    /* v5 */{
        const std::vector<vertex_id_type> ans =
            boost::assign::list_of(v1)(v2)(v4)(v6);

        std::vector<vertex_id_type> result; result.reserve(4);
        for(typename std::vector<vertex_id_type>::const_iterator
            i(poly.outgoing_edges(v5).begin()), e(poly.outgoing_edges(v5).end());
            i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(boost::algorithm::is_permutation(
                    ans.begin(), ans.end(), result.begin()));
    }
    /* v6 */{
        const std::vector<vertex_id_type> ans =
            boost::assign::list_of(v2)(v3)(v4)(v5);

        std::vector<vertex_id_type> result; result.reserve(4);
        for(typename std::vector<vertex_id_type>::const_iterator
            i(poly.outgoing_edges(v6).begin()), e(poly.outgoing_edges(v6).end());
            i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(boost::algorithm::is_permutation(
                    ans.begin(), ans.end(), result.begin()));
    }
}

//! test data 3: plane
// below, the normal vector towords the depth of your display.
//
// each edge has length 2.
// +--> x
// | 0 __1__2__3__4__ 5
// |  |\ |\ |\ |\ |\ |
// v 6|_\|_\|_\|_\|_\|11
// y  |\ |\ |\ |\ |\ |
//   .|_\|_\|_\|_\|_\| .
//   .|\ |\ |\ |\ |\ | .
//   .|_\|_\|_\|_\|_\| .
//    |\ |\ |\ |\ |\ |
//  24|_\|_\|_\|_\|_\|29
//    |\ |\ |\ |\ |\ |
//  30|_\|_\|_\|_\|_\|35
struct plane
{
    const static Real3 edge_length;

    static Polygon make()
    {
        std::vector<Triangle> triangles;
        for(std::size_t j=0; j<5; ++j)
        {
            const Real y_up = 2.0 * (j+1);
            const Real y_lw = 2.0 *  j;
            for(std::size_t i=0; i<5; ++i)
            {
                const Real x_up = 2.0 * (i+1);
                const Real x_lw = 2.0 *  i;

                const Real3 uxuy = Real3(x_up, y_up, 5.0);
                const Real3 uxly = Real3(x_up, y_lw, 5.0);
                const Real3 lxuy = Real3(x_lw, y_up, 5.0);
                const Real3 lxly = Real3(x_lw, y_lw, 5.0);

                triangles.push_back(Triangle(lxly, uxly, uxuy));
                triangles.push_back(Triangle(uxuy, lxuy, lxly));
            }
        }
        return Polygon(edge_length, triangles);
    }
};
const Real3 plane::edge_length = Real3(10, 10, 10);

BOOST_AUTO_TEST_CASE(Polygon_plane_construction_from_triangles)
{
    const Real pi = boost::math::constants::pi<Real>();
    const Polygon poly = plane::make();

    // check shape detection
    BOOST_CHECK_EQUAL(poly.face_size(),   50);
    BOOST_CHECK_EQUAL(poly.edge_size(),   50 * 3);
    BOOST_CHECK_EQUAL(poly.vertex_size(), 25);
    BOOST_CHECK_CLOSE(poly.total_area(),  10 * 10, 1e-8);

    // check vertex positions
    for(std::size_t j=0; j<5; ++j)
    {
        const Real y = 2.0 * j;
        for(std::size_t i=0; i<5; ++i)
        {
            const Real x = 2.0 * i;
            BOOST_CHECK(static_cast<bool>(poly.find_vertex(Real3(x, y, 5.0))));
        }
    }

    // check all the vertices has 6 outgoing-edges under the PBC
    {
        const std::vector<vertex_id_type> vids = poly.list_vertex_ids();
        assert(vids.size() == poly.vertex_size());

        for(std::vector<vertex_id_type>::const_iterator
                i(vids.begin()), e(vids.end()); i!=e; ++i)
        {
            BOOST_CHECK_EQUAL(poly.outgoing_edges(*i).size(), 6);
        }
    }

    // check normal vector
    {
        const std::vector<face_id_type> fids = poly.list_face_ids();
        assert(fids.size() == poly.face_size());

        for(std::vector<face_id_type>::const_iterator
                i(fids.begin()), e(fids.end()); i!=e; ++i)
        {
            BOOST_CHECK_SMALL(poly.face_at(*i).triangle.normal()[0], 1e-8);
            BOOST_CHECK_SMALL(poly.face_at(*i).triangle.normal()[1], 1e-8);
            BOOST_CHECK_CLOSE(poly.face_at(*i).triangle.normal()[2], 1.0, 1e-8);
        }
    }
    // check areas
    {
        const std::vector<face_id_type> fids = poly.list_face_ids();
        assert(fids.size() == poly.face_size());

        for(std::vector<face_id_type>::const_iterator
                i(fids.begin()), e(fids.end()); i!=e; ++i)
        {
            BOOST_CHECK_CLOSE(poly.face_at(*i).triangle.area(), 2.0, 1e-8);
        }
    }


    // check opposite edges
    {
        const std::vector<vertex_id_type> vids = poly.list_vertex_ids();
        assert(vids.size() == poly.vertex_size());

        for(std::vector<vertex_id_type>::const_iterator
                i(vids.begin()), e(vids.end()); i!=e; ++i)
        {
            const vertex_id_type vid = *i;
            std::vector<edge_id_type> const& outs = poly.outgoing_edges(vid);
            for(std::vector<edge_id_type>::const_iterator
                    oi(outs.begin()), oe(outs.end()); oi != oe; ++oi)
            {
                BOOST_CHECK_EQUAL(poly.target_of(poly.opposite_of(*oi)), vid);
            }
        }
    }

    // check next edges
    {
        const std::vector<vertex_id_type> vids = poly.list_vertex_ids();
        assert(vids.size() == poly.vertex_size());

        for(std::vector<vertex_id_type>::const_iterator
                i(vids.begin()), e(vids.end()); i!=e; ++i)
        {
            const vertex_id_type vid = *i;
            std::vector<edge_id_type> const& outs = poly.outgoing_edges(vid);
            for(std::vector<edge_id_type>::const_iterator
                    oi(outs.begin()), oe(outs.end()); oi != oe; ++oi)
            {
                BOOST_CHECK_EQUAL(
                        poly.target_of(poly.next_of(poly.next_of(*oi))), vid);
            }
        }
    }

    // check apex angle; everywhere is flat
    {
        const std::vector<vertex_id_type> vids = poly.list_vertex_ids();
        assert(vids.size() == poly.vertex_size());

        for(std::vector<vertex_id_type>::const_iterator
                i(vids.begin()), e(vids.end()); i!=e; ++i)
        {
            BOOST_CHECK_CLOSE(poly.apex_angle_at(*i), 2.0 * pi, 1e-8);
        }
    }

    // check tilt angle
    {
        const std::vector<edge_id_type> eids = poly.list_edge_ids();
        assert(eids.size() == poly.edge_size());

        for(std::vector<edge_id_type>::const_iterator
                i(eids.begin()), e(eids.end()); i!=e; ++i)
        {
            BOOST_CHECK_SMALL(poly.tilt_angle_at(*i), 1e-8);
        }
    }
}
