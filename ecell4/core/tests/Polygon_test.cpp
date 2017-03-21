#define BOOST_TEST_MODULE "Polygon_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/serialization/strong_typedef.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/STLPolygonAdapter.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Triangle.hpp>
#include <utility>

struct dummy{};

struct test_polygon_traits
{
    typedef std::size_t       index_type;
    typedef ecell4::Triangle  triangle_type;
    typedef dummy             vertex_descripter;
    typedef dummy             edge_descripter;
    typedef dummy             face_descripter;

    BOOST_STRONG_TYPEDEF(index_type, local_idx_type)
    BOOST_STRONG_TYPEDEF(index_type, face_id_type)

    typedef std::pair<face_id_type, local_idx_type> fid_localidx_pair;

    BOOST_STRONG_TYPEDEF(fid_localidx_pair, vertex_id_type)
    BOOST_STRONG_TYPEDEF(fid_localidx_pair, edge_id_type)

    static inline vertex_id_type
    make_vid(const face_id_type& fid, const local_idx_type& lidx)
    {
        return vertex_id_type(std::make_pair(fid, lidx));
    }

    static inline edge_id_type
    make_eid(const face_id_type& fid, const local_idx_type& lidx)
    {
        return edge_id_type(std::make_pair(fid, lidx));
    }

    template<typename T>
    static inline face_id_type   get_face_id(const T& id)
    {
        return static_cast<fid_localidx_pair>(id).first;
    }
    template<typename T>
    static inline local_idx_type get_local_index(const T& id)
    {
        return static_cast<fid_localidx_pair>(id).second;
    }
};

typedef ecell4::Real3 Real3;
typedef ecell4::Triangle Triangle;
typedef ecell4::Polygon<test_polygon_traits> polygon_type;
typedef ecell4::STLPolygonAdapter<test_polygon_traits> adapter_type;
typedef typename polygon_type::index_type index_type;
typedef typename polygon_type::face_id_type face_id_type;
typedef typename polygon_type::edge_id_type edge_id_type;
typedef typename polygon_type::vertex_id_type vertex_id_type;

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

    retval.add_face( t1);
    retval.add_face( t2);
    retval.add_face( t3);
    retval.add_face( t4);
    retval.add_face( t5);
    retval.add_face( t6);
    retval.add_face( t7);
    retval.add_face( t8);
    retval.add_face( t9);
    retval.add_face(t10);
    retval.add_face(t11);
    retval.add_face(t12);

    adapter_type adapter;
    adapter.detect_edge_connections(retval);
    adapter.detect_vertex_connections(retval);

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


BOOST_AUTO_TEST_CASE(Polygon_edge_connection)
{
    // TODO simplify
    const polygon_type poly = make_cube();
    {
        boost::array<face_id_type, 3> const& faces =
            poly.connecting_faces(face_id_type(0));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(1),
                                           face_id_type(7),
                                           face_id_type(8));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.connecting_faces(face_id_type(1));

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
            poly.connecting_faces(face_id_type(2));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(1),
                                           face_id_type(3),
                                           face_id_type(11));
        BOOST_CHECK(result[0] != faces.end());// failed
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.connecting_faces(face_id_type(3));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(2),
                                           face_id_type(4),
                                           face_id_type(9));
        BOOST_CHECK(result[0] != faces.end()); // failed
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.connecting_faces(face_id_type(4));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(3),
                                           face_id_type(5),
                                           face_id_type(9));
        BOOST_CHECK(result[0] != faces.end());// failed
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.connecting_faces(face_id_type(5));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(4),
                                           face_id_type(6),
                                           face_id_type(10));
        BOOST_CHECK(result[0] != faces.end());// failed
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }

    {
        boost::array<face_id_type, 3> const& faces =
            poly.connecting_faces(face_id_type(6));

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
            poly.connecting_faces(face_id_type(7));

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
            poly.connecting_faces(face_id_type(8));

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
            poly.connecting_faces(face_id_type(9));

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
            poly.connecting_faces(face_id_type(10));

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
            poly.connecting_faces(face_id_type(11));

        boost::array<boost::array<face_id_type, 3>::const_iterator, 3>
            result = find_neighbors(faces, face_id_type(1),
                                           face_id_type(2),
                                           face_id_type(10));
        BOOST_CHECK(result[0] != faces.end());
        BOOST_CHECK(result[1] != faces.end());
        BOOST_CHECK(result[2] != faces.end());
    }
}
