#define BOOST_TEST_MODULE "Polygon_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/Polygon.hpp>
#include <boost/assign.hpp>
#include <utility>

using ecell4::Real;
using ecell4::Real3;
using ecell4::Triangle;
using ecell4::Polygon;
typedef Polygon::FaceID FaceID;
typedef Polygon::EdgeID EdgeID;
typedef Polygon::VertexID VertexID;

bool check_equal(const Real3& lhs, const Real3& rhs, const Real tol)
{
    return std::abs(lhs[0] - rhs[0]) < tol &&
           std::abs(lhs[1] - rhs[1]) < tol &&
           std::abs(lhs[2] - rhs[2]) < tol;
}

namespace ecell
{
// Q1. Why is this here? use boost::algorithm or c++11 is_permutation.
// A1. Boost 1.54 (Travis.CI default) uses boost/tr1/tuple to use std::tr1::tie
// in boost::algorithm::is_permutation. But ecell4 uses <tr1/tuple> through
// <tr1/functional>, and each library(GCC C++ standard library and Boost) has
// its original implementation for std::tr1::tuple. It cause multiple-definition
// problem! To avoid this, impelement is_permutation without std::tr1::tuple.
//
// Q2. Why is this in the `ecell` namespace?
// A2. If the environment has c++11 standard library, std::is_permutation become
// a candidate because of the ADL. The argument will be a std::vector::iterator,
// and std::is_permutation is found in the `std` namespace. To avoid this
// ambiguity, put it in the namespace ecell.
template<typename Iterator1, typename Iterator2>
bool is_permutation(const Iterator1 first1, const Iterator1 last1,
                    const Iterator2 first2, const Iterator2 last2)
{
    if(std::distance(first1, last1) != std::distance(first2, last2))
    {
        return false;
    }
    if(first1 == last1) {return true;}

    for(Iterator1 i(first1); i != last1; ++i)
    {
        const std::size_t num_in_1 = std::count(first1, last1, *i);
        const std::size_t num_in_2 = std::count(first2, last2, *i);
        if(num_in_1 != num_in_2) {return false;}
    }
    return true;
}
} // detail

//! test data 1: tetrahedron
// below, the normal vector towords the depth of your display.
//
//          _4
//    3__--- /
//   /|\  4 /
//  /3|1\  /
// /__|__\/
//4  1|  /2
//    |2/
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
    BOOST_CHECK_EQUAL(poly.face_size(), 4u);
    BOOST_CHECK_EQUAL(poly.edge_size(), 6u * 2u);
    BOOST_CHECK_EQUAL(poly.vertex_size(), 4u);
    BOOST_CHECK_CLOSE(poly.total_area(), 1.5 + std::sqrt(3) / 2, 1e-8);

    // check vertex position
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(tetrahedron::p1)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(tetrahedron::p2)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(tetrahedron::p3)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(tetrahedron::p4)));

    const VertexID v1 = *poly.find_vertex(tetrahedron::p1);
    const VertexID v2 = *poly.find_vertex(tetrahedron::p2);
    const VertexID v3 = *poly.find_vertex(tetrahedron::p3);
    const VertexID v4 = *poly.find_vertex(tetrahedron::p4);

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
        const std::vector<VertexID> ans = boost::assign::list_of(v2)(v3)(v4);

        std::vector<VertexID> result; result.reserve(3);
        const std::vector<EdgeID> es = poly.outgoing_edges(v1);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }
    {
        const std::vector<VertexID> ans = boost::assign::list_of(v1)(v3)(v4);

        std::vector<VertexID> result; result.reserve(3);
        const std::vector<EdgeID> es = poly.outgoing_edges(v2);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }
    {
        const std::vector<VertexID> ans = boost::assign::list_of(v1)(v2)(v4);

        std::vector<VertexID> result; result.reserve(3);
        const std::vector<EdgeID> es = poly.outgoing_edges(v3);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }
    {
        const std::vector<VertexID> ans = boost::assign::list_of(v1)(v2)(v3);

        std::vector<VertexID> result; result.reserve(3);
        const std::vector<EdgeID> es = poly.outgoing_edges(v4);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }

    // check all the vertex are in contact with the correct set of faces.
    {
        const std::vector<VertexID> ans = boost::assign::list_of(v1)(v2)(v3);

        std::vector<VertexID> result; result.reserve(3);
        const std::vector<EdgeID> es = poly.outgoing_edges(v4);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
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

    const EdgeID e12 = *poly.find_edge(v1, v2);
    const EdgeID e13 = *poly.find_edge(v1, v3);
    const EdgeID e14 = *poly.find_edge(v1, v4);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e12))), v1);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e13))), v1);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e14))), v1);

    const EdgeID e21 = *poly.find_edge(v2, v1);
    const EdgeID e23 = *poly.find_edge(v2, v3);
    const EdgeID e24 = *poly.find_edge(v2, v4);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e21))), v2);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e23))), v2);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e24))), v2);

    const EdgeID e31 = *poly.find_edge(v3, v1);
    const EdgeID e32 = *poly.find_edge(v3, v2);
    const EdgeID e34 = *poly.find_edge(v3, v4);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e31))), v3);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e32))), v3);
    BOOST_CHECK_EQUAL(poly.target_of(poly.next_of(poly.next_of(e34))), v3);

    const EdgeID e41 = *poly.find_edge(v4, v1);
    const EdgeID e42 = *poly.find_edge(v4, v2);
    const EdgeID e43 = *poly.find_edge(v4, v3);
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


    // test of distance
    const FaceID f1 = *(poly.find_face(v1, v2, v3));
    const FaceID f2 = *(poly.find_face(v1, v2, v4));
    const FaceID f3 = *(poly.find_face(v1, v3, v4));
    const FaceID f4 = *(poly.find_face(v2, v3, v4));
    {
        const Real3 p1(0, 0, 0);
        const Real3 p2(0, 0, 1);

        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p1, f1), std::make_pair(p2, f2)),
                1.0, 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p2, f2), std::make_pair(p1, f1)),
                1.0, 1e-8);

        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p1, f1), std::make_pair(p2, f3)),
                1.0, 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p2, f3), std::make_pair(p1, f1)),
                1.0, 1e-8);

        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p1, f1), std::make_pair(p2, f4)),
                1.0, 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p2, f4), std::make_pair(p1, f1)),
                1.0, 1e-8);
    }

    {
        const Real3 p1(1.0/3.0, 1.0/3.0, 0);
        const Real3 p2(1.0/3.0, 1.0/3.0, 1.0/3.0);

        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p1, f1), std::make_pair(p2, f4)),
                (std::sqrt(2) + std::sqrt(6)) / 6.0, 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p2, f4), std::make_pair(p1, f1)),
                (std::sqrt(2) + std::sqrt(6)) / 6.0, 1e-8);

    }
}

//! test data 2: octahedron
// below, the normal vector towords the depth of your display.
//            3
//            ^
//           / \
//          / 2 \
// 3______1/_____\4______3
//  \     /\     /\     /\
//   \ 1 /  \ 3 /  \ 7 /  \
//    \ / 4  \ /  6 \ /  8 \
//     v______v______v______\
//    2       5\     /6      2
//              \ 5 /
//               \ /
//                v
//                2
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
    BOOST_CHECK_EQUAL(poly.face_size(), 8u);
    BOOST_CHECK_EQUAL(poly.edge_size(), 8u * 3u);
    BOOST_CHECK_EQUAL(poly.vertex_size(), 6u);
    BOOST_CHECK_CLOSE(poly.total_area(), 8 * std::sqrt(3.0) / 2, 1e-8);

    // check vertex position
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p1)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p2)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p3)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p4)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p5)));
    BOOST_CHECK(static_cast<bool>(poly.find_vertex(octahedron::p6)));

    const VertexID v1 = *poly.find_vertex(octahedron::p1);
    const VertexID v2 = *poly.find_vertex(octahedron::p2);
    const VertexID v3 = *poly.find_vertex(octahedron::p3);
    const VertexID v4 = *poly.find_vertex(octahedron::p4);
    const VertexID v5 = *poly.find_vertex(octahedron::p5);
    const VertexID v6 = *poly.find_vertex(octahedron::p6);

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
        const std::vector<VertexID> ans =
            boost::assign::list_of(v2)(v3)(v4)(v5);

        std::vector<VertexID> result; result.reserve(4);
        const std::vector<EdgeID> es = poly.outgoing_edges(v1);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }
    /* v2 */{
        const std::vector<VertexID> ans =
            boost::assign::list_of(v1)(v3)(v5)(v6);

        std::vector<VertexID> result; result.reserve(4);
        const std::vector<EdgeID> es = poly.outgoing_edges(v2);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }
    /* v3 */{
        const std::vector<VertexID> ans =
            boost::assign::list_of(v1)(v2)(v4)(v6);

        std::vector<VertexID> result; result.reserve(4);
        const std::vector<EdgeID> es = poly.outgoing_edges(v3);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }
    /* v4 */{
        const std::vector<VertexID> ans =
            boost::assign::list_of(v1)(v3)(v5)(v6);

        std::vector<VertexID> result; result.reserve(4);
        const std::vector<EdgeID> es = poly.outgoing_edges(v4);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }
    /* v5 */{
        const std::vector<VertexID> ans =
            boost::assign::list_of(v1)(v2)(v4)(v6);

        std::vector<VertexID> result; result.reserve(4);
        const std::vector<EdgeID> es = poly.outgoing_edges(v5);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }
    /* v6 */{
        const std::vector<VertexID> ans =
            boost::assign::list_of(v2)(v3)(v4)(v5);

        std::vector<VertexID> result; result.reserve(4);
        const std::vector<EdgeID> es = poly.outgoing_edges(v6);
        for(typename std::vector<EdgeID>::const_iterator
                i(es.begin()), e(es.end()); i != e; ++i)
        {
            result.push_back(poly.target_of(*i));
        }
        BOOST_CHECK(ans.size() == result.size());
        BOOST_CHECK(ecell::is_permutation(
                    ans.begin(), ans.end(), result.begin(), result.end()));
    }

    // test of distance
    const FaceID f1 = *(poly.find_face(v1, v2, v3));
    const FaceID f2 = *(poly.find_face(v1, v3, v4));
    const FaceID f3 = *(poly.find_face(v1, v4, v5));
    const FaceID f4 = *(poly.find_face(v1, v5, v2));
    const FaceID f5 = *(poly.find_face(v6, v2, v5));
    const FaceID f6 = *(poly.find_face(v6, v5, v4));
    const FaceID f7 = *(poly.find_face(v6, v4, v3));
    const FaceID f8 = *(poly.find_face(v6, v3, v2));
    {
        const Real3 p1 = (octahedron::p1 + octahedron::p2 + octahedron::p3) / 3.0;
        const Real3 p2 = (octahedron::p1 + octahedron::p4 + octahedron::p5) / 3.0;

        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p1, f1), std::make_pair(p2, f3)),
                std::sqrt(2.0), 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p2, f3), std::make_pair(p1, f1)),
                std::sqrt(2.0), 1e-8);
    }
    {
        const Real3 p1 = (octahedron::p1 + octahedron::p3 + octahedron::p4) / 3.0;
        const Real3 p2 = (octahedron::p1 + octahedron::p2 + octahedron::p5) / 3.0;

        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p1, f2), std::make_pair(p2, f4)),
                std::sqrt(2.0), 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p2, f4), std::make_pair(p1, f2)),
                std::sqrt(2.0), 1e-8);
    }
    {
        const Real3 p1 = (octahedron::p2 + octahedron::p6 + octahedron::p5) / 3.0;
        const Real3 p2 = (octahedron::p3 + octahedron::p4 + octahedron::p6) / 3.0;

        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p1, f5), std::make_pair(p2, f7)),
                std::sqrt(2.0), 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p2, f7), std::make_pair(p1, f5)),
                std::sqrt(2.0), 1e-8);
    }
    {
        const Real3 p1 = (octahedron::p4 + octahedron::p5 + octahedron::p6) / 3.0;
        const Real3 p2 = (octahedron::p2 + octahedron::p3 + octahedron::p6) / 3.0;

        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p1, f6), std::make_pair(p2, f8)),
                std::sqrt(2.0), 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(
                poly.distance(std::make_pair(p2, f8), std::make_pair(p1, f6)),
                std::sqrt(2.0), 1e-8);
    }

    {
        const Real3 p1 = (octahedron::p1 + octahedron::p2 + octahedron::p3) / 3.0;
        const Real3 p2 = (octahedron::p4 + octahedron::p5 + octahedron::p6) / 3.0;

        BOOST_CHECK_EQUAL(
                poly.distance(std::make_pair(p1, f1), std::make_pair(p2, f6)),
                std::numeric_limits<Real>::infinity());
        BOOST_CHECK_EQUAL(
                poly.distance(std::make_pair(p2, f6), std::make_pair(p1, f1)),
                std::numeric_limits<Real>::infinity());
    }

    // travel
    {
        const Real3 p1 = (octahedron::p1 + octahedron::p2 + octahedron::p3) / 3.0;
        const Real3 p2 = (octahedron::p2 * 2 + octahedron::p3 * 2 - octahedron::p1) / 3.0;

        const std::pair<Real3, FaceID> p2_ =
            poly.travel(std::make_pair(p1, f1), p2 - p1);

        BOOST_CHECK_EQUAL(p2_.second, f8);

        const Real3 dst = (octahedron::p2 + octahedron::p3 + octahedron::p6) / 3.0;
        BOOST_CHECK_CLOSE(p2_.first[0], dst[0], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[1], dst[1], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[2], dst[2], 1e-6);

        const Real dist = poly.distance(std::make_pair(p1, f1), p2_);
        BOOST_CHECK_CLOSE(dist, length(p2 - p1), 1e-6);
    }
    {
        const Real3 p1 = (octahedron::p1 + octahedron::p2 + octahedron::p3) / 3.0;
        const Real3 p2 = (octahedron::p1 * 2 + octahedron::p2 * 2 - octahedron::p3) / 3.0;

        const std::pair<Real3, FaceID> p2_ =
            poly.travel(std::make_pair(p1, f1), p2 - p1);

        BOOST_CHECK_EQUAL(p2_.second, f4);

        const Real3 dst = (octahedron::p1 + octahedron::p2 + octahedron::p5) / 3.0;
        BOOST_CHECK_CLOSE(p2_.first[0], dst[0], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[1], dst[1], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[2], dst[2], 1e-6);

        const Real dist = poly.distance(std::make_pair(p1, f1), p2_);
        BOOST_CHECK_CLOSE(dist, length(p2 - p1), 1e-6);
    }
    {
        const Real3 p1 = (octahedron::p1 + octahedron::p2 + octahedron::p3) / 3.0;
        const Real3 p2 = (octahedron::p1 * 2 + octahedron::p3 * 2 - octahedron::p2) / 3.0;

        const std::pair<Real3, FaceID> p2_ =
            poly.travel(std::make_pair(p1, f1), p2 - p1);

        BOOST_CHECK_EQUAL(p2_.second, f2);

        const Real3 dst = (octahedron::p1 + octahedron::p3 + octahedron::p4) / 3.0;
        BOOST_CHECK_CLOSE(p2_.first[0], dst[0], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[1], dst[1], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[2], dst[2], 1e-6);

        const Real dist = poly.distance(std::make_pair(p1, f1), p2_);
        BOOST_CHECK_CLOSE(dist, length(p2 - p1), 1e-6);
    }

    // roll
    {
        const VertexID vid = v1;
        const Real       r = std::sqrt(2.0 / 3.0);
        const Real3 p1 = (octahedron::p1 + octahedron::p2 + octahedron::p3) / 3.0;
        const Real3 p2 = (octahedron::p1 + octahedron::p3 + octahedron::p4) / 3.0;
        const Real3 p3 = (octahedron::p1 + octahedron::p4 + octahedron::p5) / 3.0;
        const Real3 p4 = (octahedron::p1 + octahedron::p2 + octahedron::p5) / 3.0;

        {
            const std::pair<Real3, FaceID> no_move =
                ::ecell4::polygon::roll(poly, std::make_pair(p1, f1), v1, r, 0.0);
            BOOST_CHECK_EQUAL(no_move.second, f1);
            BOOST_CHECK_CLOSE(no_move.first[0], p1[0], 1e-6);
            BOOST_CHECK_CLOSE(no_move.first[1], p1[1], 1e-6);
            BOOST_CHECK_CLOSE(no_move.first[2], p1[2], 1e-6);
        }
        const std::pair<Real3, FaceID> p2_ = ::ecell4::polygon::roll(poly, std::make_pair(p1, f1), v1, r,  60.0 / 180.0 * pi);
        const std::pair<Real3, FaceID> p3_ = ::ecell4::polygon::roll(poly, std::make_pair(p1, f1), v1, r, 120.0 / 180.0 * pi);
        const std::pair<Real3, FaceID> p4_ = ::ecell4::polygon::roll(poly, std::make_pair(p1, f1), v1, r, 180.0 / 180.0 * pi);
        const std::pair<Real3, FaceID> p1_ = ::ecell4::polygon::roll(poly, std::make_pair(p1, f1), v1, r, 240.0 / 180.0 * pi);
        BOOST_CHECK_EQUAL(p2_.second, f2);
        BOOST_CHECK_EQUAL(p3_.second, f3);
        BOOST_CHECK_EQUAL(p4_.second, f4);
        BOOST_CHECK_EQUAL(p1_.second, f1);

        BOOST_CHECK_CLOSE(p1_.first[0], p1[0], 1e-6);
        BOOST_CHECK_CLOSE(p1_.first[1], p1[1], 1e-6);
        BOOST_CHECK_CLOSE(p1_.first[2], p1[2], 1e-6);

        BOOST_CHECK_CLOSE(p2_.first[0], p2[0], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[1], p2[1], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[2], p2[2], 1e-6);

        BOOST_CHECK_CLOSE(p3_.first[0], p3[0], 1e-6);
        BOOST_CHECK_CLOSE(p3_.first[1], p3[1], 1e-6);
        BOOST_CHECK_CLOSE(p3_.first[2], p3[2], 1e-6);

        BOOST_CHECK_CLOSE(p4_.first[0], p4[0], 1e-6);
        BOOST_CHECK_CLOSE(p4_.first[1], p4[1], 1e-6);
        BOOST_CHECK_CLOSE(p4_.first[2], p4[2], 1e-6);
    }

    {
        const VertexID vid = v1;
        const Real       r = std::sqrt(2.0 / 3.0);
        const Real3 p1 = (octahedron::p1 + octahedron::p2 + octahedron::p3) / 3.0;
        const Real3 p2 = (octahedron::p1 + octahedron::p3 + octahedron::p4) / 3.0;
        const Real3 p3 = (octahedron::p1 + octahedron::p4 + octahedron::p5) / 3.0;
        const Real3 p4 = (octahedron::p1 + octahedron::p2 + octahedron::p5) / 3.0;

        {
            const std::pair<Real3, FaceID> no_move =
                ::ecell4::polygon::roll(poly, std::make_pair(p3, f3), v1, r, 0.0);
            BOOST_CHECK_EQUAL(no_move.second, f3);
            BOOST_CHECK_CLOSE(no_move.first[0], p3[0], 1e-6);
            BOOST_CHECK_CLOSE(no_move.first[1], p3[1], 1e-6);
            BOOST_CHECK_CLOSE(no_move.first[2], p3[2], 1e-6);
        }
        const std::pair<Real3, FaceID> p4_ = ::ecell4::polygon::roll(poly, std::make_pair(p3, f3), v1, r,  60.0 / 180.0 * pi);
        const std::pair<Real3, FaceID> p1_ = ::ecell4::polygon::roll(poly, std::make_pair(p3, f3), v1, r, 120.0 / 180.0 * pi);
        const std::pair<Real3, FaceID> p2_ = ::ecell4::polygon::roll(poly, std::make_pair(p3, f3), v1, r, 180.0 / 180.0 * pi);
        const std::pair<Real3, FaceID> p3_ = ::ecell4::polygon::roll(poly, std::make_pair(p3, f3), v1, r, 240.0 / 180.0 * pi);
        BOOST_CHECK_EQUAL(p2_.second, f2);
        BOOST_CHECK_EQUAL(p3_.second, f3);
        BOOST_CHECK_EQUAL(p4_.second, f4);
        BOOST_CHECK_EQUAL(p1_.second, f1);

        BOOST_CHECK_CLOSE(p1_.first[0], p1[0], 1e-6);
        BOOST_CHECK_CLOSE(p1_.first[1], p1[1], 1e-6);
        BOOST_CHECK_CLOSE(p1_.first[2], p1[2], 1e-6);

        BOOST_CHECK_CLOSE(p2_.first[0], p2[0], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[1], p2[1], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[2], p2[2], 1e-6);

        BOOST_CHECK_CLOSE(p3_.first[0], p3[0], 1e-6);
        BOOST_CHECK_CLOSE(p3_.first[1], p3[1], 1e-6);
        BOOST_CHECK_CLOSE(p3_.first[2], p3[2], 1e-6);

        BOOST_CHECK_CLOSE(p4_.first[0], p4[0], 1e-6);
        BOOST_CHECK_CLOSE(p4_.first[1], p4[1], 1e-6);
        BOOST_CHECK_CLOSE(p4_.first[2], p4[2], 1e-6);
    }
}

//! test data 3: plane
// below, the normal vector towords the depth of your display.
//
// each edge has length 2.
// +--> x
// | 0 __1__2__3__4__ 0
// |  |\ |\ |\ |\ |\ |
// v 5|_\|_\|_\|_\9_\|5
// y  |\ |\ |\ |\ |\ |
//   .|_\|_\|_\|_\|_\| .
//   .|\ |\ |\ |\ |\ | .
//   .|_\|_\|_\|_\|_\| .
//    |\ |\ |\ |\ |\ |
//  20|_\|_\|_\|_\24\|20
//    |\ |\ |\ |\ |\ |
//   0|_\|_\|_\|_\|_\|0
//
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
    BOOST_CHECK_EQUAL(poly.face_size(),   50u);
    BOOST_CHECK_EQUAL(poly.edge_size(),   50u * 3u);
    BOOST_CHECK_EQUAL(poly.vertex_size(), 25u);
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
        const std::vector<VertexID> vids = poly.list_vertex_ids();
        assert(vids.size() == poly.vertex_size());

        for(std::vector<VertexID>::const_iterator
                i(vids.begin()), e(vids.end()); i!=e; ++i)
        {
            BOOST_CHECK_EQUAL(poly.outgoing_edges(*i).size(), 6u);
        }
    }

    // check normal vector
    {
        const std::vector<FaceID> fids = poly.list_face_ids();
        assert(fids.size() == poly.face_size());

        for(std::vector<FaceID>::const_iterator
                i(fids.begin()), e(fids.end()); i!=e; ++i)
        {
            BOOST_CHECK_SMALL(poly.triangle_at(*i).normal()[0], 1e-8);
            BOOST_CHECK_SMALL(poly.triangle_at(*i).normal()[1], 1e-8);
            BOOST_CHECK_CLOSE(poly.triangle_at(*i).normal()[2], 1.0, 1e-8);
        }
    }
    // check areas
    {
        const std::vector<FaceID> fids = poly.list_face_ids();
        assert(fids.size() == poly.face_size());

        for(std::vector<FaceID>::const_iterator
                i(fids.begin()), e(fids.end()); i!=e; ++i)
        {
            BOOST_CHECK_CLOSE(poly.triangle_at(*i).area(), 2.0, 1e-8);
        }
    }


    // check opposite edges
    {
        const std::vector<VertexID> vids = poly.list_vertex_ids();
        assert(vids.size() == poly.vertex_size());

        for(std::vector<VertexID>::const_iterator
                i(vids.begin()), e(vids.end()); i!=e; ++i)
        {
            const VertexID vid = *i;
            std::vector<EdgeID> const& outs = poly.outgoing_edges(vid);
            for(std::vector<EdgeID>::const_iterator
                    oi(outs.begin()), oe(outs.end()); oi != oe; ++oi)
            {
                BOOST_CHECK_EQUAL(poly.target_of(poly.opposite_of(*oi)), vid);
            }
        }
    }

    // check next edges
    {
        const std::vector<VertexID> vids = poly.list_vertex_ids();
        assert(vids.size() == poly.vertex_size());

        for(std::vector<VertexID>::const_iterator
                i(vids.begin()), e(vids.end()); i!=e; ++i)
        {
            const VertexID vid = *i;
            std::vector<EdgeID> const& outs = poly.outgoing_edges(vid);
            for(std::vector<EdgeID>::const_iterator
                    oi(outs.begin()), oe(outs.end()); oi != oe; ++oi)
            {
                BOOST_CHECK_EQUAL(
                        poly.target_of(poly.next_of(poly.next_of(*oi))), vid);
            }
        }
    }

    // check apex angle; everywhere is flat
    {
        const std::vector<VertexID> vids = poly.list_vertex_ids();
        assert(vids.size() == poly.vertex_size());

        for(std::vector<VertexID>::const_iterator
                i(vids.begin()), e(vids.end()); i!=e; ++i)
        {
            BOOST_CHECK_CLOSE(poly.apex_angle_at(*i), 2.0 * pi, 1e-8);
        }
    }

    // check tilt angle
    {
        const std::vector<EdgeID> eids = poly.list_edge_ids();
        assert(eids.size() == poly.edge_size());

        for(std::vector<EdgeID>::const_iterator
                i(eids.begin()), e(eids.end()); i!=e; ++i)
        {
            BOOST_CHECK_SMALL(poly.tilt_angle_at(*i), 1e-8);
        }
    }

    // distance stuff
    //
    // 24___20____21
    //  |\   |\   |
    //  | \12| \14|
    //  |11\ |13\ |
    // 4|___\0___\1____2
    //  |\   |\   |\   |
    //  | \6 | \2 | \4 |
    //  | 5\ | 1\ | 3\ |
    // 9|___\5___\6___\7
    //       |\   |\   |
    //       | \8 | \10|
    //       | 7\ | 9\ |
    //     10|___\11__\|12

    const VertexID v0  = *poly.find_vertex(Real3(0.0, 0.0, 5.0));
    const VertexID v1  = *poly.find_vertex(Real3(2.0, 0.0, 5.0));
    const VertexID v2  = *poly.find_vertex(Real3(4.0, 0.0, 5.0));
    const VertexID v4  = *poly.find_vertex(Real3(8.0, 0.0, 5.0));

    const VertexID v5  = *poly.find_vertex(Real3(0.0, 2.0, 5.0));
    const VertexID v6  = *poly.find_vertex(Real3(2.0, 2.0, 5.0));
    const VertexID v7  = *poly.find_vertex(Real3(4.0, 2.0, 5.0));
    const VertexID v9  = *poly.find_vertex(Real3(8.0, 2.0, 5.0));

    const VertexID v10 = *poly.find_vertex(Real3(0.0, 4.0, 5.0));
    const VertexID v11 = *poly.find_vertex(Real3(2.0, 4.0, 5.0));
    const VertexID v12 = *poly.find_vertex(Real3(4.0, 4.0, 5.0));

    const VertexID v20 = *poly.find_vertex(Real3(0.0, 8.0, 5.0));
    const VertexID v21 = *poly.find_vertex(Real3(2.0, 8.0, 5.0));
    const VertexID v24 = *poly.find_vertex(Real3(8.0, 8.0, 5.0));

    const FaceID f1 = *poly.find_face(v0, v5, v6);
    const FaceID f2 = *poly.find_face(v0, v1, v6);

    const FaceID f3 = *poly.find_face(v1, v6, v7);
    const FaceID f4 = *poly.find_face(v1, v2, v7);

    const FaceID f5 = *poly.find_face(v4, v5, v9);
    const FaceID f6 = *poly.find_face(v0, v4, v5);

    const FaceID f7 = *poly.find_face(v5, v10, v11);
    const FaceID f8 = *poly.find_face(v5, v6,  v11);

    const FaceID f9  = *poly.find_face(v6, v11, v12);
    const FaceID f10 = *poly.find_face(v6, v7,  v12);

    const FaceID f11 = *poly.find_face(v0, v4,  v24);
    const FaceID f12 = *poly.find_face(v0, v20, v24);

    const FaceID f13 = *poly.find_face(v0, v1,  v20);
    const FaceID f14 = *poly.find_face(v1, v20, v21);

    {
        const EdgeID e0_01 = *poly.find_edge(v0, v1);
        const EdgeID e0_20 = *poly.find_edge(v0, v20);
        const EdgeID e0_24 = *poly.find_edge(v0, v24);
        const EdgeID e0_04 = *poly.find_edge(v0, v4);
        const EdgeID e0_05 = *poly.find_edge(v0, v5);
        const EdgeID e0_06 = *poly.find_edge(v0, v6);
        BOOST_CHECK_EQUAL(poly.opposite_of(poly.next_of(poly.next_of(e0_01))), e0_06);
        BOOST_CHECK_EQUAL(poly.opposite_of(poly.next_of(poly.next_of(e0_06))), e0_05);
        BOOST_CHECK_EQUAL(poly.opposite_of(poly.next_of(poly.next_of(e0_05))), e0_04);
        BOOST_CHECK_EQUAL(poly.opposite_of(poly.next_of(poly.next_of(e0_04))), e0_24);
        BOOST_CHECK_EQUAL(poly.opposite_of(poly.next_of(poly.next_of(e0_24))), e0_20);
        BOOST_CHECK_EQUAL(poly.opposite_of(poly.next_of(poly.next_of(e0_20))), e0_01);
    }

    {
        const Real3 p1(0.5, 1.5, 5.0);
        const Real3 p2(1.5, 0.5, 5.0);
        BOOST_CHECK_CLOSE_FRACTION(poly.distance(
            std::make_pair(p1, f1), std::make_pair(p2, f2)),
            length(poly.periodic_transpose(p1, p2) - p2), 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(poly.distance(
            std::make_pair(p2, f2), std::make_pair(p1, f1)),
            length(poly.periodic_transpose(p1, p2) - p2), 1e-8);
    }

    {
        const Real3 p1(0.5, 1.5, 5.0);
        const Real3 p2(8.5, 9.5, 5.0);
        BOOST_CHECK_CLOSE_FRACTION(poly.distance(
            std::make_pair(p1, f1), std::make_pair(p2, f11)),
            length(poly.periodic_transpose(p1, p2) - p2), 1e-8);
        BOOST_CHECK_CLOSE_FRACTION(poly.distance(
            std::make_pair(p2, f11), std::make_pair(p1, f1)),
            length(poly.periodic_transpose(p1, p2) - p2), 1e-8);
    }

    {
        const Real3 p1(0.5, 1.5, 5.0);
        const Real3 p2(1.5, 0.5, 5.0);

        const Real3 v1to2 = poly.direction(
                std::make_pair(p1, f1), std::make_pair(p2, f2));
        const Real3 v2to1 = poly.direction(
                std::make_pair(p2, f2), std::make_pair(p1, f1));
        BOOST_CHECK_CLOSE(v1to2[0],  1.0, 1e-6);
        BOOST_CHECK_CLOSE(v1to2[1], -1.0, 1e-6);
        BOOST_CHECK_SMALL(v1to2[2], 1e-6);

        BOOST_CHECK_CLOSE(v2to1[0], -1.0, 1e-6);
        BOOST_CHECK_CLOSE(v2to1[1],  1.0, 1e-6);
        BOOST_CHECK_SMALL(v2to1[2], 1e-6);
    }

    {
        const Real3 p1(0.5, 1.5, 5.0);
        const Real3 p2(8.5, 9.5, 5.0);
        const Real3 v1to2 = poly.direction(
            std::make_pair(p1, f1), std::make_pair(p2, f11));
        const Real3 v2to1 = poly.direction(
            std::make_pair(p2, f11), std::make_pair(p1, f1));

        BOOST_CHECK_CLOSE(v1to2[0], -2.0, 1e-6);
        BOOST_CHECK_CLOSE(v1to2[1], -2.0, 1e-6);
        BOOST_CHECK_SMALL(v1to2[2], 1e-6);

        BOOST_CHECK_CLOSE(v2to1[0], 2.0, 1e-6);
        BOOST_CHECK_CLOSE(v2to1[1], 2.0, 1e-6);
        BOOST_CHECK_SMALL(v2to1[2], 1e-6);
    }

    // traveling ----------------------------------------------------------

    {
        const Real3 p1(0.5, 1.5, 5.0);
        const Real3 p2(1.5, 0.5, 5.0);

        const std::pair<Real3, FaceID> p2_ =
            poly.travel(std::make_pair(p1, f1), Real3(1.0, -1.0, 0.0));
        const std::pair<Real3, FaceID> p1_ =
            poly.travel(std::make_pair(p2, f2), Real3(-1.0, 1.0, 0.0));

        BOOST_CHECK_EQUAL(p1_.second, f1);
        BOOST_CHECK_EQUAL(p2_.second, f2);

        BOOST_CHECK_CLOSE(p1_.first[0], p1[0], 1e-6);
        BOOST_CHECK_CLOSE(p1_.first[1], p1[1], 1e-6);
        BOOST_CHECK_CLOSE(p1_.first[2], p1[2], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[0], p2[0], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[1], p2[1], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[2], p2[2], 1e-6);

    }

    {
        const Real3 p1(0.5, 1.5, 5.0);
        const Real3 p2(8.5, 9.5, 5.0);

        const std::pair<Real3, FaceID> p2_ =
            poly.travel(std::make_pair(p1, f1),  Real3(-2.0, -2.0, 0.0));
        const std::pair<Real3, FaceID> p1_ =
            poly.travel(std::make_pair(p2, f11), Real3( 2.0,  2.0, 0.0));

        BOOST_CHECK_EQUAL(p1_.second, f1);
        BOOST_CHECK_EQUAL(p2_.second, f11);

        BOOST_CHECK_CLOSE(p1_.first[0], p1[0], 1e-6);
        BOOST_CHECK_CLOSE(p1_.first[1], p1[1], 1e-6);
        BOOST_CHECK_CLOSE(p1_.first[2], p1[2], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[0], p2[0], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[1], p2[1], 1e-6);
        BOOST_CHECK_CLOSE(p2_.first[2], p2[2], 1e-6);
    }

    // roll ----------------------------------------------------------

    {
        // around v6, far from Boundary
        const Real3 p1(2.0 - std::sqrt(3.0), 1.0, 5.0);
        const std::pair<Real3, FaceID> p2 =
            ::ecell4::polygon::roll(poly, std::make_pair(p1, f1), v6, 2, 30.0 / 180.0 * pi);

        BOOST_CHECK_EQUAL(p2.second, f2);
        BOOST_CHECK_CLOSE(p2.first[0], 1.0,                  1e-6);
        BOOST_CHECK_CLOSE(p2.first[1], 2.0 - std::sqrt(3.0), 1e-6);
        BOOST_CHECK_CLOSE(p2.first[2], 5.0, 1e-6);

        const std::pair<Real3, FaceID> p10 =
            ::ecell4::polygon::roll(poly, std::make_pair(p1, f1), v6, 2, pi);
        BOOST_CHECK_EQUAL(p10.second, f10);
        BOOST_CHECK_CLOSE(p10.first[0], 2.0 + std::sqrt(3.0), 1e-6);
        BOOST_CHECK_CLOSE(p10.first[1], 3.0, 1e-6);
        BOOST_CHECK_CLOSE(p10.first[2], 5.0, 1e-6);
    }
    {
        const Real3 p1(1.0, std::sqrt(3.0), 5.0);
        const Real3 p2(std::sqrt(3.0), 1.0, 5.0);

        const std::pair<Real3, FaceID> p11 =
            ::ecell4::polygon::roll(poly, std::make_pair(p2, f2), v0, 2, pi);
        BOOST_CHECK_EQUAL(p11.second, f11);
        BOOST_CHECK_CLOSE(p11.first[0], 10.0 - std::sqrt(3.0), 1e-6);
        BOOST_CHECK_CLOSE(p11.first[1], 9.0, 1e-6);
        BOOST_CHECK_CLOSE(p11.first[2], 5.0, 1e-6);

        const std::pair<Real3, FaceID> p12 =
            ::ecell4::polygon::roll(poly, std::make_pair(p1, f1), v0, 2, pi);
        BOOST_CHECK_EQUAL(p12.second, f12);
        BOOST_CHECK_CLOSE(p12.first[0], 9.0, 1e-6);
        BOOST_CHECK_CLOSE(p12.first[1], 10.0 - std::sqrt(3.0), 1e-6);
        BOOST_CHECK_CLOSE(p12.first[2], 5.0, 1e-6);
    }
}
