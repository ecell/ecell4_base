#ifndef ECELL4_SGFRD_POLYGON_TRAITS
#define ECELL4_SGFRD_POLYGON_TRAITS
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/Polygon.hpp>
#include <boost/serialization/strong_typedef.hpp>
#include <utility>
#include <limits>
#include <vector>

namespace ecell4
{

namespace sgfrd
{

namespace polygon
{

BOOST_STRONG_TYPEDEF(std::size_t, face_id_type)
BOOST_STRONG_TYPEDEF(std::size_t, vertex_id_type)
BOOST_STRONG_TYPEDEF(std::size_t, edge_id_type)

struct vertex_descripter
{
    Real max_conical_shell_size;
    std::vector<face_id_type>   neighbor_faces;
    std::vector<vertex_id_type> neighbor_vertices; // include self
};

struct edge_descripter
{
    // empty.
};

struct face_descripter
{
    // developped edges of neighbor faces
    boost::array<std::pair<Real3, Real3>, 6> segments_must_not_collide;

    std::vector<face_id_type>   neighbor_faces;
    std::vector<vertex_id_type> neighbor_vertices;
};

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
    T_id to_id(const std::size_t& i) const {return T_id(i);}

    template<typename T_id>
    void link(const T_id&, std::size_t){return;}
};
}// namespace polygon

struct polygon_traits
{
    typedef std::size_t index_type;
    typedef ecell4::sgfrd::polygon::face_id_type   face_id_type;
    typedef ecell4::sgfrd::polygon::edge_id_type   edge_id_type;
    typedef ecell4::sgfrd::polygon::vertex_id_type vertex_id_type;

    typedef Triangle    triangle_type;
    typedef ecell4::sgfrd::polygon::vertex_descripter vertex_descripter;
    typedef ecell4::sgfrd::polygon::edge_descripter   edge_descripter;
    typedef ecell4::sgfrd::polygon::face_descripter   face_descripter;
    typedef ecell4::sgfrd::polygon::index_generator<
            face_id_type, edge_id_type, vertex_id_type> id_generator_type;
    typedef ecell4::sgfrd::polygon::identity_converter converter_type;

    template<typename Tid>
    static inline Tid un_initialized()
    {
        return Tid(std::numeric_limits<std::size_t>::max());
    }
};

//XXX impl of these function is dirty...
polygon::vertex_descripter
make_vertex_information(const ecell4::Polygon<polygon_traits>& poly,
                        const polygon::vertex_id_type vid);
polygon::face_descripter
make_face_information(const ecell4::Polygon<polygon_traits>& poly,
                      const polygon::face_id_type fid);

void setup_descriptors(ecell4::Polygon<polygon_traits>& poly);

} // sgfrd
} // ecell4

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<ecell4::sgfrd::polygon::face_id_type>
{
    typedef ecell4::sgfrd::polygon::face_id_type argument_type;

    std::size_t operator()(argument_type const& val) const
    {
        return hash<std::size_t>()(static_cast<std::size_t>(val));
    }
};

template<>
struct hash<ecell4::sgfrd::polygon_traits::edge_id_type>
{
    typedef ecell4::sgfrd::polygon::edge_id_type argument_type;

    std::size_t operator()(argument_type const& val) const
    {
        return hash<std::size_t>()(static_cast<std::size_t>(val));
    }
};

template<>
struct hash<ecell4::sgfrd::polygon::vertex_id_type>
{
    typedef ecell4::sgfrd::polygon::vertex_id_type argument_type;

    std::size_t operator()(argument_type const& val) const
    {
        return hash<std::size_t>()(static_cast<std::size_t>(val));
    }
};

ECELL4_DEFINE_HASH_END()

#endif // ECELL4_SGFRD_POLYGON_TRAITS
