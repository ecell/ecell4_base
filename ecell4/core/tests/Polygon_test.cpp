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
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/STLPolygonAdapter.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Triangle.hpp>
#include <utility>

struct dummy{};

template<typename T_fid, typename T_eid, typename T_vid>
struct index_generator
{
    index_generator(): face_(0), edge_(0), vertex_(0){}

    T_fid generate_face_id()  {return T_fid(face_++);}
    T_eid generate_edge_id()  {return T_eid(edge_++);}
    T_vid generate_vertex_id(){return T_vid(vertex_++);}

    T_fid face_;
    T_eid edge_;
    T_vid vertex_;
};

struct identity_converter
{
    template<typename T_id>
    std::size_t to_index(const T_id& id) const {return static_cast<std::size_t>(id);}

    template<typename T_id>
    T_id to_index(const std::size_t& i) const {return T_id(i);}

    template<typename T_id>
    void link(const T_id&, std::size_t){return;}
};

struct test_polygon_traits
{
    typedef std::size_t       index_type;
    typedef ecell4::Triangle  triangle_type;
    typedef dummy             vertex_descripter;
    typedef dummy             edge_descripter;
    typedef dummy             face_descripter;

    BOOST_STRONG_TYPEDEF(index_type, face_id_type)
    BOOST_STRONG_TYPEDEF(index_type, vertex_id_type)
    BOOST_STRONG_TYPEDEF(index_type, edge_id_type)

    typedef index_generator<face_id_type, edge_id_type, vertex_id_type>
        id_generator_type;
    typedef identity_converter converter_type;

    template<typename Tid>
    static inline Tid un_initialized()
    {
        return Tid(std::numeric_limits<std::size_t>::max());
    }
};

typedef ecell4::Real  Real;
typedef ecell4::Real3 Real3;
typedef ecell4::Triangle Triangle;
typedef ecell4::Polygon<test_polygon_traits> polygon_type;
typedef ecell4::STLPolygonAdapter<test_polygon_traits> adapter_type;
typedef polygon_type::index_type index_type;
typedef polygon_type::face_id_type face_id_type;
typedef polygon_type::edge_id_type edge_id_type;
typedef polygon_type::vertex_id_type vertex_id_type;

polygon_type make_cube()
{
    polygon_type retval;

    const Real3 v0(0, 0,  0);
    const Real3 v1(1, 0,  0);
    const Real3 v2(1, 1,  0);
    const Real3 v3(0, 1,  0);
    const Real3 v4(0, 0, -1);
    const Real3 v5(1, 0, -1);
    const Real3 v6(1, 1, -1);
    const Real3 v7(0, 1, -1);

    const Triangle  t1(v0, v1, v3);
    const Triangle  t2(v1, v2, v3);
    const Triangle  t3(v2, v6, v3);
    const Triangle  t4(v3, v6, v7);
    const Triangle  t5(v7, v6, v4);
    const Triangle  t6(v4, v6, v5);
    const Triangle  t7(v4, v5, v1);
    const Triangle  t8(v4, v1, v0);
    const Triangle  t9(v4, v0, v3);
    const Triangle t10(v4, v3, v7);
    const Triangle t11(v1, v5, v6);
    const Triangle t12(v1, v6, v2);

    std::vector<face_id_type> fids(12);
    fids[ 0] = retval.add_face( t1);
    fids[ 1] = retval.add_face( t2);
    fids[ 2] = retval.add_face( t3);
    fids[ 3] = retval.add_face( t4);
    fids[ 4] = retval.add_face( t5);
    fids[ 5] = retval.add_face( t6);
    fids[ 6] = retval.add_face( t7);
    fids[ 7] = retval.add_face( t8);
    fids[ 8] = retval.add_face( t9);
    fids[ 9] = retval.add_face(t10);
    fids[10] = retval.add_face(t11);
    fids[11] = retval.add_face(t12);

    adapter_type adapter;
    adapter.detect_edge_connections(retval, fids);
    adapter.detect_vertex_connections(retval, fids);

    return retval;
}

boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
find_neighbors(boost::array<face_id_type, 3> const& fs,
        const face_id_type& f1, const face_id_type& f2, const face_id_type& f3)
{
    boost::array<boost::array<face_id_type, 3>::const_iterator, 3> retval;
    retval[0] = std::find(fs.begin(), fs.end(), f1);
    retval[1] = std::find(fs.begin(), fs.end(), f2);
    retval[2] = std::find(fs.begin(), fs.end(), f3);
    return retval;
}

BOOST_AUTO_TEST_CASE(Polygon_num_stuff)
{
    const polygon_type poly = make_cube();

    BOOST_CHECK_EQUAL(poly.num_faces(),     12);
    BOOST_CHECK_EQUAL(poly.num_triangles(), 12);
    BOOST_CHECK_EQUAL(poly.num_edges(),     18);
    BOOST_CHECK_EQUAL(poly.num_vertices(),  8);
}


BOOST_AUTO_TEST_CASE(Polygon_neighbor_faces)
{
    using namespace boost::assign;
    typedef face_id_type f;
    const polygon_type poly = make_cube();
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(0));
        std::vector<face_id_type> expects;
        expects += f(1),f(2),f(3),f(6),f(7),f(8),f(9),f(10),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(1));
        std::vector<face_id_type> expects;
        expects += f(0),f(2),f(3),f(6),f(7),f(8),f(9),f(10),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(2));
        std::vector<face_id_type> expects;
        expects += f(0),f(1),f(3),f(4),f(5),f(8),f(9),f(10),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(3));
        std::vector<face_id_type> expects;
        expects += f(0),f(1),f(2),f(4),f(5),f(8),f(9),f(10),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(4));
        std::vector<face_id_type> expects;
        expects += f(2),f(3),f(5),f(6),f(7),f(8),f(9),f(10),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(5));
        std::vector<face_id_type> expects;
        expects += f(2),f(3),f(4),f(6),f(7),f(8),f(9),f(10),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(6));
        std::vector<face_id_type> expects;
        expects += f(0),f(1),f(4),f(5),f(7),f(8),f(9),f(10),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(7));
        std::vector<face_id_type> expects;
        expects += f(0),f(1),f(4),f(5),f(6),f(8),f(9),f(10),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(8));
        std::vector<face_id_type> expects;
        expects += f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(9);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(9));
        std::vector<face_id_type> expects;
        expects += f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(8);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(10));
        std::vector<face_id_type> expects;
        expects += f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(11);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
    {
        std::vector<face_id_type> const& neighbors = poly.neighbor_faces(f(11));
        std::vector<face_id_type> expects;
        expects += f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7),f(10);
        BOOST_CHECK_EQUAL(neighbors.size(), expects.size());
        const bool result = boost::algorithm::is_permutation(
                neighbors.begin(), neighbors.end(), expects.begin());
        BOOST_CHECK(result);
    }
}


BOOST_AUTO_TEST_CASE(Polygon_edge_connection)
{
    // TODO simplify
    const polygon_type poly = make_cube();
    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(0));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(1),
                                           face_id_type(7),
                                           face_id_type(8));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    //{{{
    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(1));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(0),
                                           face_id_type(2),
                                           face_id_type(11));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(2));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(1),
                                           face_id_type(3),
                                           face_id_type(11));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(3));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(2),
                                           face_id_type(4),
                                           face_id_type(9));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(4));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(3),
                                           face_id_type(5),
                                           face_id_type(9));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(5));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(4),
                                           face_id_type(6),
                                           face_id_type(10));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(6));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(7),
                                           face_id_type(10),
                                           face_id_type(5));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(7));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(0),
                                           face_id_type(8),
                                           face_id_type(6));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(8));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(0),
                                           face_id_type(7),
                                           face_id_type(9));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(9));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(3),
                                           face_id_type(4),
                                           face_id_type(8));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(10));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(5),
                                           face_id_type(6),
                                           face_id_type(11));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.adjacent_faces(face_id_type(11));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(1),
                                           face_id_type(2),
                                           face_id_type(10));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }
    //}}}
}

BOOST_AUTO_TEST_CASE(Polygon_distance_same_face)
{
    const polygon_type poly = make_cube();
    typedef polygon_type::face_id_type face_id_type;

    const Real dist_same_face1 = poly.distance(
            std::make_pair(Real3(0.5, 0.1, 0.), face_id_type(0)),
            std::make_pair(Real3(0.1, 0.5, 0.), face_id_type(0))
            );
    const Real ref_same_face1 = length(Real3(0.5, 0.1, 0.) - Real3(0.1, 0.5, 0.));
    BOOST_CHECK_CLOSE(dist_same_face1, ref_same_face1, 1e-12);

    const Real dist_same_face2 = poly.distance(
            std::make_pair(Real3(0.5, 0.0, 0.), face_id_type(0)),
            std::make_pair(Real3(0.0, 0.5, 0.), face_id_type(0))
            );
    const Real ref_same_face2 = length(Real3(0.5, 0.0, 0.) - Real3(0.0, 0.5, 0.));
    BOOST_CHECK_CLOSE(dist_same_face2, ref_same_face2, 1e-12);

    const Real dist_same_face3 = poly.distance(
            std::make_pair(Real3(1.0, 0.0, 0.), face_id_type(0)),
            std::make_pair(Real3(0.0, 1.0, 0.), face_id_type(0))
            );
    const Real ref_same_face3 = length(Real3(1.0, 0.0, 0.) - Real3(0.0, 1.0, 0.));
    BOOST_CHECK_CLOSE(dist_same_face3, ref_same_face3, 1e-12);
}

BOOST_AUTO_TEST_CASE(Polygon_distance_connected_by_edge)
{
    const polygon_type poly = make_cube();
    typedef polygon_type::face_id_type face_id_type;

    const Real dist1_1 = poly.distance(
            std::make_pair(Real3(0.3, 0.3, 0.), face_id_type(0)),
            std::make_pair(Real3(0.7, 0.7, 0.), face_id_type(1)));
    const Real dist1_2 = poly.distance(
            std::make_pair(Real3(0.7, 0.7, 0.), face_id_type(1)),
            std::make_pair(Real3(0.3, 0.3, 0.), face_id_type(0)));
    const Real ref1 = length(Real3(0.3, 0.3, 0.) - Real3(0.7, 0.7, 0.));
    BOOST_CHECK_CLOSE(dist1_1, ref1, 1e-12);
    BOOST_CHECK_CLOSE(dist1_2, ref1, 1e-12);

    const Real dist2_1 = poly.distance(
            std::make_pair(Real3(0.3, 0.3, 0.), face_id_type(0)),
            std::make_pair(Real3(0.3, 0.0, -0.3), face_id_type(7)));
    const Real dist2_2 = poly.distance(
            std::make_pair(Real3(0.3, 0.0, -0.3), face_id_type(7)),
            std::make_pair(Real3(0.3, 0.3, 0.), face_id_type(0)));
    const Real ref2 = 0.6;
    BOOST_CHECK_CLOSE(dist2_1, ref2, 1e-12);
    BOOST_CHECK_CLOSE(dist2_2, ref2, 1e-12);
}

BOOST_AUTO_TEST_CASE(Polygon_distance_connected_by_vertex)
{
    const polygon_type poly = make_cube();
    typedef polygon_type::face_id_type face_id_type;

    const Real dist1_1 = poly.distance(
            std::make_pair(Real3(0.3, 0.3,  0.0), face_id_type(0)),
            std::make_pair(Real3(1.0, 0.3, -0.7), face_id_type(10)));

    const Real dist1_2 = poly.distance(
            std::make_pair(Real3(1.0, 0.3, -0.7), face_id_type(10)),
            std::make_pair(Real3(0.3, 0.3,  0.0), face_id_type(0)));
    const Real ref1 = 1.4;
    BOOST_CHECK_CLOSE(dist1_1, ref1, 1e-12);
    BOOST_CHECK_CLOSE(dist1_2, ref1, 1e-12);
}

BOOST_AUTO_TEST_CASE(Polygon_developped_direction_same_face)
{
    const polygon_type poly = make_cube();
    typedef polygon_type::face_id_type face_id_type;

    const Real3 dir1 = poly.developed_direction(
            std::make_pair(Real3(0.5, 0.1, 0.), face_id_type(0)),
            std::make_pair(Real3(0.1, 0.5, 0.), face_id_type(0)));
    const Real dist1 = poly.distance(
            std::make_pair(Real3(0.5, 0.1, 0.), face_id_type(0)),
            std::make_pair(Real3(0.1, 0.5, 0.), face_id_type(0)));
    const std::pair<std::pair<Real3, face_id_type>, Real3>
        next1 = poly.move_next_face(
                std::make_pair(Real3(0.5, 0.1, 0.), face_id_type(0)), dir1);
    BOOST_CHECK_CLOSE(dist1, length(dir1), 1e-12);
    BOOST_CHECK_EQUAL(next1.first.second, face_id_type(0));
    BOOST_CHECK_SMALL(length(next1.second), 1e-12);
    BOOST_CHECK_SMALL(length(next1.first.first - Real3(0.1, 0.5, 0.)), 1e-12);

    const Real3 dir2 = poly.developed_direction(
            std::make_pair(Real3(0.5, 0.0, 0.), face_id_type(0)),
            std::make_pair(Real3(0.0, 0.5, 0.), face_id_type(0)));
    const Real dist2 = poly.distance(
            std::make_pair(Real3(0.5, 0.0, 0.), face_id_type(0)),
            std::make_pair(Real3(0.0, 0.5, 0.), face_id_type(0)));
    const std::pair<std::pair<Real3, face_id_type>, Real3>
        next2 = poly.move_next_face(
                std::make_pair(Real3(0.5, 0.0, 0.), face_id_type(0)), dir2);
    BOOST_CHECK_CLOSE(dist2, length(dir2), 1e-12);
    BOOST_CHECK_EQUAL(next2.first.second, face_id_type(0));
    BOOST_CHECK_SMALL(length(next2.second), 1e-12);
    BOOST_CHECK_SMALL(length(next2.first.first - Real3(0.0, 0.5, 0.)), 1e-12);


    const Real3 dir3 = poly.developed_direction(
            std::make_pair(Real3(1.0, 0.0, 0.), face_id_type(0)),
            std::make_pair(Real3(0.0, 1.0, 0.), face_id_type(0)));
    const Real dist3 = poly.distance(
            std::make_pair(Real3(1.0, 0.0, 0.), face_id_type(0)),
            std::make_pair(Real3(0.0, 1.0, 0.), face_id_type(0)));
    const std::pair<std::pair<Real3, face_id_type>, Real3>
        next3 = poly.move_next_face(
                std::make_pair(Real3(1.0, 0.0, 0.), face_id_type(0)), dir3);
    BOOST_CHECK_CLOSE(dist3, length(dir3), 1e-12);
    BOOST_CHECK_EQUAL(next3.first.second, face_id_type(0));
    BOOST_CHECK_SMALL(length(next3.second), 1e-12);
    BOOST_CHECK_SMALL(length(next3.first.first - Real3(0.0, 1.0, 0.)), 1e-12);
}


BOOST_AUTO_TEST_CASE(Polygon_developped_direction_connected_by_edge)
{
    const polygon_type poly = make_cube();
    typedef polygon_type::face_id_type face_id_type;

    const std::pair<Real3, face_id_type> start1 =
            std::make_pair(Real3(0.3, 0.3, 0.), face_id_type(0));
    const std::pair<Real3, face_id_type> term1  =
            std::make_pair(Real3(0.7, 0.7, 0.), face_id_type(1));

    const Real3 dir1_1 = poly.developed_direction(start1, term1);
    const Real dist1_1 = poly.distance(start1, term1);
    BOOST_CHECK_CLOSE(dist1_1, length(dir1_1), 1e-12);

    const std::pair<std::pair<Real3, face_id_type>, Real3>
        tmp1_1 = poly.move_next_face(start1, dir1_1);
    const std::pair<std::pair<Real3, face_id_type>, Real3>
        next1_1 = poly.move_next_face(tmp1_1.first, tmp1_1.second);
    BOOST_CHECK_SMALL(length(next1_1.second), 1e-12);
    BOOST_CHECK_EQUAL(next1_1.first.second, face_id_type(1));
    BOOST_CHECK_SMALL(length(next1_1.first.first - term1.first), 1e-12);

    const Real3 dir1_2 = poly.developed_direction(term1, start1);
    const Real dist1_2 = poly.distance(term1, start1);
    BOOST_CHECK_CLOSE(dist1_2, length(dir1_2), 1e-12);

    const std::pair<std::pair<Real3, face_id_type>, Real3>
        tmp1_2 = poly.move_next_face(term1, dir1_2);
    const std::pair<std::pair<Real3, face_id_type>, Real3>
        next1_2 = poly.move_next_face(tmp1_2.first, tmp1_2.second);
    BOOST_CHECK_SMALL(length(next1_2.second), 1e-12);
    BOOST_CHECK_EQUAL(next1_2.first.second, face_id_type(0));
    BOOST_CHECK_SMALL(length(next1_2.first.first - start1.first), 1e-12);

    // -------------------------------------------------------------------

    const std::pair<Real3, face_id_type> start2 =
            std::make_pair(Real3(0.3, 0.3, 0.), face_id_type(0));
    const std::pair<Real3, face_id_type> term2  =
            std::make_pair(Real3(0.3, 0.0, -0.3), face_id_type(7));

    const Real3 dir2_1 = poly.developed_direction(start2, term2);
    const Real dist2_1 = poly.distance(start2, term2);
    BOOST_CHECK_CLOSE(dist2_1, length(dir2_1), 1e-12);

    const std::pair<std::pair<Real3, face_id_type>, Real3>
        tmp2_1 = poly.move_next_face(start2, dir2_1);
    const std::pair<std::pair<Real3, face_id_type>, Real3>
        next2_1 = poly.move_next_face(tmp2_1.first, tmp2_1.second);
    BOOST_CHECK_SMALL(length(next2_1.second), 1e-12);
    BOOST_CHECK_EQUAL(next2_1.first.second, term2.second);
    BOOST_CHECK_SMALL(length(next2_1.first.first - term2.first), 1e-12);

    const Real3 dir2_2 = poly.developed_direction(term2, start2);
    const Real dist2_2 = poly.distance(term2, start2);
    BOOST_CHECK_CLOSE(dist2_2, length(dir2_2), 1e-12);

    const std::pair<std::pair<Real3, face_id_type>, Real3>
        tmp2_2 = poly.move_next_face(term2, dir2_2);
    const std::pair<std::pair<Real3, face_id_type>, Real3>
        next2_2 = poly.move_next_face(tmp2_2.first, tmp2_2.second);
    BOOST_CHECK_SMALL(length(next2_2.second), 1e-12);
    BOOST_CHECK_EQUAL(next2_2.first.second, start2.second);
    BOOST_CHECK_SMALL(length(next2_2.first.first - start2.first), 1e-12);
}

BOOST_AUTO_TEST_CASE(Polygon_developped_direction_connected_by_vertex)
{
    const polygon_type poly = make_cube();
    typedef polygon_type::face_id_type face_id_type;

    const std::pair<Real3, face_id_type> start1 =
            std::make_pair(Real3(0.3, 0.3,  0.0), face_id_type(0));
    const std::pair<Real3, face_id_type> term1 =
            std::make_pair(Real3(1.0, 0.3, -0.7), face_id_type(10));

    const Real3 dir1_1 = poly.developed_direction(start1, term1);
    const Real dist1_1 = poly.distance(start1, term1);
    BOOST_CHECK_CLOSE(dist1_1, length(dir1_1), 1e-12);

    std::pair<std::pair<Real3, face_id_type>, Real3>
        next1_1 = poly.move_next_face(start1, dir1_1);
    while(length(next1_1.second) < 1e-12)
    {
        BOOST_CHECK_EQUAL(next1_1.first.second, term1.second);
        BOOST_CHECK_SMALL(length(next1_1.first.first - term1.first), 1e-12);
    }

    // -------------------------------------------------------------------

    const Real3 dir1_2 = poly.developed_direction(
            std::make_pair(Real3(1.0, 0.3, -0.7), face_id_type(10)),
            std::make_pair(Real3(0.3, 0.3,  0.0), face_id_type(0)));
    const Real dist1_2 = poly.distance(
            std::make_pair(Real3(1.0, 0.3, -0.7), face_id_type(10)),
            std::make_pair(Real3(0.3, 0.3,  0.0), face_id_type(0)));
    BOOST_CHECK_CLOSE(dist1_2, length(dir1_2), 1e-12);

    std::pair<std::pair<Real3, face_id_type>, Real3>
        next1_2 = poly.move_next_face(term1, dir1_2);
    while(length(next1_2.second) < 1e-12)
    {
        BOOST_CHECK_EQUAL(next1_2.first.second, start1.second);
        BOOST_CHECK_SMALL(length(next1_2.first.first - start1.first), 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(Polygon_list_vertices_within_radius)
{
    const polygon_type polygon = make_cube();

    const Real3 pos1(0.2, 0.4, 0.0);
    const face_id_type fid1(0);
    const face_id_type fid2(1);
    const face_id_type fid8(7);
    const face_id_type fid9(8);
    std::pair<std::vector<std::pair<vertex_id_type, Real> >,
        std::pair<vertex_id_type, Real> > value =
            polygon.list_vertices_within_radius(
                    std::make_pair(pos1, fid1), std::sqrt(8.0/5.0) + 1e-8);

    const std::vector<std::pair<vertex_id_type, Real> >& list = value.first;
    const std::pair<vertex_id_type, Real>& nearest = value.second;
    BOOST_CHECK_EQUAL(nearest.first, polygon.get_vertex_id(std::make_pair(fid1, 0)));
    BOOST_CHECK_CLOSE_FRACTION(nearest.second, std::sqrt(1./5.), 1e-8);

    BOOST_CHECK_EQUAL(list.size(), 5);
    BOOST_CHECK_EQUAL(list.at(0).first, polygon.get_vertex_id(std::make_pair(fid1,0)));
    BOOST_CHECK_EQUAL(list.at(1).first, polygon.get_vertex_id(std::make_pair(fid1,2)));
    BOOST_CHECK_EQUAL(list.at(2).first, polygon.get_vertex_id(std::make_pair(fid1,1)));
    BOOST_CHECK_EQUAL(list.at(3).first, polygon.get_vertex_id(std::make_pair(fid2,1)));
    BOOST_CHECK_EQUAL(list.at(4).first, polygon.get_vertex_id(std::make_pair(fid9,0)));
    BOOST_CHECK_CLOSE_FRACTION(list.at(0).second, std::sqrt(1./5.), 1e-8);
    BOOST_CHECK_CLOSE_FRACTION(list.at(1).second, std::sqrt(2./5.), 1e-8);
    BOOST_CHECK_CLOSE_FRACTION(list.at(2).second, std::sqrt(4./5.), 1e-8);
    BOOST_CHECK_CLOSE_FRACTION(list.at(3).second, 1.0,              1e-8);
    BOOST_CHECK_CLOSE_FRACTION(list.at(4).second, std::sqrt(8./5.), 1e-8);
}
