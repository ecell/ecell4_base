#ifndef ECELL_BD_POLYGON
#define ECELL_BD_POLYGON
#include "Barycentric.hpp"
#include "map_const_at.hpp"
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/hash.hpp>
#include <functional>
#include <utility>
#include <vector>
#include <stdint.h>

namespace ecell4
{
namespace bd
{

// definition of face_id, edge_id, vertex_id and its utility function
struct polygon_traits
{
    typedef Triangle    face_type;
    typedef std::size_t face_index_type;
    typedef std::size_t face_id_type;
    typedef uint32_t    vertex_index_type;
    typedef uint32_t    edge_index_type;
    typedef std::pair<face_id_type, vertex_index_type> vertex_id_type;
    typedef std::pair<face_id_type, edge_index_type>   edge_id_type;

    static vertex_id_type
    make_vertex_id(face_id_type fid, vertex_index_type vidx)
    {
        return std::make_pair(fid, vidx);
    }

    static edge_id_type
    make_edge_id(face_id_type fid, edge_index_type eidx)
    {
        return std::make_pair(fid, eidx);
    }

    static edge_id_type vtoe(const vertex_id_type vid){return vid;}
    static vertex_id_type etov(const edge_id_type eid){return eid;}

    static face_id_type get_face_id(const vertex_id_type vid){return vid.first;}
//    static face_id_type get_face_id(const edge_id_type eid){return eid.first;}
    static vertex_index_type get_vertex_index(const vertex_id_type& vid)
    {
        return vid.second;
    }
    static edge_index_type get_edge_index(const edge_id_type& eid)
    {
        return eid.second;
    }
};

} // bd
} // ecell4

// define hash for vertex/edge_id_type
ECELL4_DEFINE_HASH_BEGIN()
template<>
struct hash<std::pair<std::size_t, uint32_t> >
{
    typedef std::pair<std::size_t, uint32_t> argument_type;

    std::size_t operator()(argument_type const& val) const
    {
        return hash<std::size_t>()(val.first) ^
               hash<uint32_t>()(val.second);
    }
};
ECELL4_DEFINE_HASH_END()

namespace ecell4
{
namespace bd
{

// template<typename T_traits>
class BDPolygon
{
  public:
    typedef polygon_traits traits; // hard-coded
    typedef traits traits_type;
    typedef /*typename*/ traits_type::face_type         face_type;
    typedef /*typename*/ traits_type::face_id_type      face_id_type;
    typedef /*typename*/ traits_type::face_index_type   face_index_type;
    typedef /*typename*/ traits_type::vertex_index_type vertex_index_type;
    typedef /*typename*/ traits_type::edge_index_type   edge_index_type;
    typedef /*typename*/ traits_type::vertex_id_type    vertex_id_type;
    typedef /*typename*/ traits_type::edge_id_type      edge_id_type;

    typedef std::vector<face_type>      face_container_type;
    typedef std::vector<vertex_id_type> vertex_id_list;
    typedef /*typename*/ utils::get_mapper_mf<edge_id_type, edge_id_type>::type
        edge_pair_type;
    // map from vertex_id_type --to-> pairof(vertex_id_list, whole_angle)
    typedef std::pair<vertex_id_list, Real> vertex_list_angle_pair;
    typedef /*typename*/ utils::get_mapper_mf<vertex_id_type,
            vertex_list_angle_pair>::type vertex_group_type;

    typedef std::pair<Real3, face_id_type> surface_position_type;
    typedef Barycentric<Real> barycentric_type;

  public:
    BDPolygon(){}
    ~BDPolygon(){}

    face_id_type add_face(const face_type& face)
    {
        const face_id_type id = faces_.size();
        faces_.push_back(face);
        return id;
    }

    void detect_connectivity(const Real tolerance = 1e-6);

    edge_id_type const&
    connecting_edge(const edge_id_type& eid) const
    {
        return const_at(edge_pairs_, eid);
    }

    vertex_id_list const&
    connecting_vertices(const vertex_id_type& vid) const
    {
        return const_at(vertex_groups_, vid).first;
    }

    std::pair<bool, edge_index_type>
    is_connected(const face_id_type& lhs, const face_id_type& rhs) const;

    std::pair<bool, vertex_index_type>
    is_share_vertex(const face_id_type& lhs, const face_id_type& rhs) const;

    Real distance_sq(const std::pair<Real3, face_id_type>& lhs,
                     const std::pair<Real3, face_id_type>& rhs) const;

    Real distance(const std::pair<Real3, face_id_type>& lhs,
                  const std::pair<Real3, face_id_type>& rhs) const
    {
        return std::sqrt(distance_sq(lhs, rhs));
    }

    Real3 inter_position_vector(const std::pair<Real3, face_id_type>& lhs,
                                const std::pair<Real3, face_id_type>& rhs) const;

    std::pair<std::pair<Real3, face_id_type>, Real3>
    move_next_face(const std::pair<Real3, face_id_type>& pos,
                   const Real3& disp) const;

    std::size_t num_faces() const {return faces_.size();}
    bool empty() const {return faces_.empty();}
    void clear() {return faces_.clear();}
    face_type&       operator[](const std::size_t i)       {return faces_[i];}
    face_type const& operator[](const std::size_t i) const {return faces_[i];}
    face_type&       at(const std::size_t i)       {return faces_.at(i);}
    face_type const& at(const std::size_t i) const {return faces_.at(i);}

  private:

    void detect_shared_vertices(const Real tolerance = 1e-6);
    void detect_shared_edges(const Real tolerance = 1e-6);

    std::pair<uint32_t, Real>
    crossed_edge(const barycentric_type& pos, const barycentric_type& disp) const;
    Real
    cross_section(const barycentric_type& pos, const barycentric_type& disp,
                  const edge_index_type eidx) const;

  private:

    face_container_type faces_;
    edge_pair_type      edge_pairs_;
    vertex_group_type   vertex_groups_;
};

}// bd

}// ecell
#endif //ECELL_BD_POLYGON
