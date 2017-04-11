#ifndef ECELL4_SGFRD_POLYGON_TRAITS
#define ECELL4_SGFRD_POLYGON_TRAITS
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Triangle.hpp>
#include <boost/serialization/strong_typedef.hpp>
#include <utility>
#include <vector>

namespace ecell4
{

namespace sgfrd
{

struct vertex_descripter
{
    Real distance_to_nearest_incenter;
};

struct edge_descripter
{
};

struct face_descripter
{
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

struct polygon_traits
{
    typedef std::size_t       index_type;
    typedef Triangle          triangle_type;
    typedef ::ecell4::sgfrd::vertex_descripter vertex_descripter;
    typedef ::ecell4::sgfrd::edge_descripter   edge_descripter;
    typedef ::ecell4::sgfrd::face_descripter   face_descripter;

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

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<ecell4::bd::polygon_traits::face_id_type>
{
    typedef ecell4::bd::polygon_traits::face_id_type argument_type;

    std::size_t operator()(argument_type const& val) const
    {
        return hash<std::size_t>()(static_cast<std::size_t>(val));
    }
};

template<>
struct hash<ecell4::bd::polygon_traits::edge_id_type>
{
    typedef ecell4::bd::polygon_traits::edge_id_type argument_type;

    std::size_t operator()(argument_type const& val) const
    {
        return hash<std::size_t>()(static_cast<std::size_t>(val));
    }
};

template<>
struct hash<ecell4::bd::polygon_traits::vertex_id_type>
{
    typedef ecell4::bd::polygon_traits::vertex_id_type argument_type;

    std::size_t operator()(argument_type const& val) const
    {
        return hash<std::size_t>()(static_cast<std::size_t>(val));
    }
};

ECELL4_DEFINE_HASH_END()

} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_POLYGON_TRAITS
