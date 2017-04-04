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

template<typename Tid>
struct index_generator
{
    index_generator(): current(0){}

    Tid operator()()
    {
        return Tid(current++);
    }

    std::size_t current;
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

    typedef index_generator<face_id_type>   face_id_generator_type;
    typedef index_generator<edge_id_type>   edge_id_generator_type;
    typedef index_generator<vertex_id_type> vertex_id_generator_type;
    typedef ecell4::utils::get_mapper_mf<face_id_type, std::size_t>::type
        face_id_idx_map_type;
    typedef ecell4::utils::get_mapper_mf<edge_id_type, std::size_t>::type
        edge_id_idx_map_type;
    typedef ecell4::utils::get_mapper_mf<vertex_id_type, std::size_t>::type
        vertex_id_idx_map_type;

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
