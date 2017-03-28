#ifndef ECELL4_BD_POLYGON
#define ECELL4_BD_POLYGON
#include <ecell4/core/Polygon.hpp>
#include <boost/serialization/strong_typedef.hpp>
#include <functional>
#include <utility>
#include <vector>
#include <stdint.h>

namespace ecell4
{
namespace bd
{

struct empty{};

struct polygon_traits
{
    typedef std::size_t      index_type;
    typedef ecell4::Triangle triangle_type;
    typedef empty            vertex_descripter;
    typedef empty            edge_descripter;
    typedef empty            face_descripter;

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

typedef Polygon<polygon_traits> BDPolygon;

}// bd
}// ecell4


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

ECELL4_DEFINE_HASH_END()

#endif //ECELL_BD_POLYGON
