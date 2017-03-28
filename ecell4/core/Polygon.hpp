#ifndef ECELL4_POLYGON
#define ECELL4_POLYGON
#include <utility>
#include <functional>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/comparators.hpp>
#include <ecell4/core/Shape.hpp>
#include <boost/array.hpp>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <limits>

namespace ecell4
{

template<typename T_traits>
class Polygon : public Shape
{
    template<typename vertexT> struct vertex_property;
    template<typename edgeT>   struct edge_property;
    template<typename faceT>   struct face_property;

  public:
    typedef T_traits traits_type;
    typedef typename traits_type::index_type        index_type;
    typedef typename traits_type::triangle_type     triangle_type;
    typedef typename traits_type::face_id_type      face_id_type;
    typedef typename traits_type::local_idx_type    local_idx_type;
    typedef typename traits_type::vertex_id_type    vertex_id_type;
    typedef typename traits_type::edge_id_type      edge_id_type;
    typedef typename traits_type::vertex_descripter vertex_descripter;
    typedef typename traits_type::edge_descripter   edge_descripter;
    typedef typename traits_type::face_descripter   face_descripter;
    typedef vertex_property<vertex_descripter> vertex_property_type;
    typedef edge_property<edge_descripter>     edge_property_type;
    typedef face_property<face_descripter>     face_property_type;
    typedef std::vector<vertex_property_type> vertex_container_type;
    typedef std::vector<edge_property_type>   edge_container_type;
    typedef std::vector<face_property_type>   face_container_type;
    typedef std::vector<triangle_type>   triangle_container_type;

    static inline vertex_id_type
    make_vid(const face_id_type& fid, const local_idx_type& lidx)
    {
        return traits_type::make_vid(fid, lidx);
    }

    static inline edge_id_type
    make_eid(const face_id_type& fid, const local_idx_type& lidx)
    {
        return traits_type::make_eid(fid, lidx);
    }

    template<typename T>
    static inline face_id_type
    get_face_id(const T& id)
    {
        return traits_type::template get_face_id<T>(id);
    }

    template<typename T>
    static inline local_idx_type
    get_local_index(const T& id)
    {
        return traits_type::template get_local_index<T>(id);
    }

  public:

    Polygon(){}
    ~Polygon(){}

    /* set connection --------------------------------------------------------*/
    index_type add_face(const triangle_type& tri);
    index_type add_face(const triangle_type& tri, const face_descripter& face);

    index_type connect_edges(const edge_id_type&, const edge_id_type&);
    index_type connect_edges(const edge_id_type&, const edge_id_type&,
                             const edge_descripter&);

    index_type connect_vertices(const std::vector<vertex_id_type>&);
    index_type connect_vertices(const std::vector<vertex_id_type>&,
                                const vertex_descripter&);

    /* get connecting information --------------------------------------------*/
    std::pair<bool, edge_id_type>
    is_connected_by_edge(const face_id_type& from, const face_id_type& to) const;

    std::pair<bool, vertex_id_type>
    is_connected_by_vertex(const face_id_type& from, const face_id_type& to) const;

    std::vector<face_id_type> const&
    neighbor_faces(const face_id_type& fid) const;

    boost::array<face_id_type, 3> const&
    connecting_faces(const face_id_type& fid) const;

    edge_id_type
    connecting_edge(const edge_id_type& eid) const;

    /* geometric functions ---------------------------------------------------*/
    Real distance_sq(const std::pair<Real3, face_id_type>& lhs,
                     const std::pair<Real3, face_id_type>& rhs) const;
    Real distance(const std::pair<Real3, face_id_type>& lhs,
                  const std::pair<Real3, face_id_type>& rhs) const;
    Real3
    developed_direction(const std::pair<Real3, face_id_type>& from,
                        const std::pair<Real3, face_id_type>& to) const;

    std::pair<std::vector<vertex_id_type>, std::pair<vertex_id_type, Real> >
    list_vertices_within_radius(const std::pair<Real3, face_id_type>& pos,
                                const Real radius);

    std::pair<std::vector<vertex_id_type>, std::pair<vertex_id_type, Real> >
    list_faces_within_radius(const Real3& pos, const Real radius);

    /* dynamics functions ----------------------------------------------------*/
    std::pair<std::pair<Real3, face_id_type>, Real3>
    move_next_face(const std::pair<Real3, face_id_type>& pos,
                   const Real3& disp) const;

    std::pair<Real3, face_id_type>
    rotate_around_vertex(const std::pair<Real3, face_id_type>& pos,
        const vertex_id_type& apex, const Real r, const Real angle) const;

    /* member access ---------------------------------------------------------*/
    std::size_t num_triangles() const {return triangles_.size();}
    std::size_t num_faces()     const {return faces_.size();}
    std::size_t num_edges()     const {return edges_.size();}
    std::size_t num_vertices()  const {return vertices_.size();}

    triangle_type&           triangle_at(const face_id_type i);
    triangle_type const&     triangle_at(const face_id_type i) const;
    face_descripter&         face_at(const face_id_type i);
    face_descripter const&   face_at(const face_id_type i) const;
    vertex_descripter&       vertex_at(const vertex_id_type&);
    vertex_descripter const& vertex_at(const vertex_id_type&) const;
    edge_descripter&         edge_at(const edge_id_type& i);
    edge_descripter const&   edge_at(const edge_id_type& i) const;

    /* required by shape -----------------------------------------------------*/
    dimension_kind dimension() const {return THREE;}

    Real is_inside(const Real3& coord) const;

    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const;

    bool test_AABB(const Real3& l, const Real3& u) const;

    void bounding_box(
            const Real3& edge_lengths, Real3& lower, Real3& upper) const;

  private:

    vertex_property_type&       vertex_prop_at(const vertex_id_type& vid);
    vertex_property_type const& vertex_prop_at(const vertex_id_type& vid) const;
    edge_property_type&         edge_prop_at(const edge_id_type& eid);
    edge_property_type const&   edge_prop_at(const edge_id_type& eid) const;
    face_property_type&         face_prop_at(const face_id_type& fid);
    face_property_type const&   face_prop_at(const face_id_type& fid) const;

  private:

    static const std::size_t un_initialized;

    template<typename vertexT>
    struct vertex_property
    {
        vertexT vertex;
        std::vector<std::size_t> edges; // idx of this->edges_
        // pairof (idx of this->faces_, local idx of this->faces_->vertices)
        std::vector<std::pair<std::size_t, std::size_t> > faces;
    };

    template<typename edgeT>
    struct edge_property
    {
        edgeT edge;
        std::pair<std::size_t, std::size_t> vertices;// idx of this->vertices_
        // pairof (idx of this->faces_, local idx of this->faces_->vertices)
        std::pair<std::pair<std::size_t, std::size_t>,
                  std::pair<std::size_t, std::size_t> > faces;
    };

    template<typename faceT>
    struct face_property
    {
        faceT      face;
        index_type triangle_index;               // idx of this->triangles_
        boost::array<std::size_t, 3>  vertices;  // idx of this->vertices_
        boost::array<std::size_t, 3>  edges;     // idx of this->edges_
        boost::array<face_id_type, 3> faces;     // idx of this->faces_
        std::vector<face_id_type>     neighbors; // idx of this->faces_
    };

    triangle_container_type triangles_;
    vertex_container_type   vertices_;
    edge_container_type     edges_;
    face_container_type     faces_;
};

template<typename T>
const std::size_t
Polygon<T>::un_initialized = std::numeric_limits<std::size_t>::max();

template<typename T>
inline Real
Polygon<T>::is_inside(const Real3& coord) const
{
    throw NotImplemented("ecell4::Polygon::is_inside");
}

template<typename T>
inline Real3
Polygon<T>::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    throw NotImplemented("ecell4::Polygon::draw_position");
}

template<typename T>
inline bool
Polygon<T>::test_AABB(const Real3& l, const Real3& u) const
{
    throw NotImplemented("ecell4::Polygon::test_AABB");
}

template<typename T>
inline void
Polygon<T>::bounding_box(const Real3& edge_lengths, Real3& l, Real3& u) const
{
    throw NotImplemented("ecell4::Polygon::bounding_box");
}

template<typename T>
typename Polygon<T>::index_type
Polygon<T>::add_face(const triangle_type& face)
{
    const index_type idx = this->triangles_.size();
    assert(idx == this->faces_.size());
    this->triangles_.push_back(face);

    face_property_type fp;
    fp.triangle_index = idx;
    fp.vertices.fill(un_initialized);
    fp.edges.fill(un_initialized);
    fp.faces.fill(face_id_type(un_initialized));
    this->faces_.push_back(fp);

    return idx;
}

template<typename T>
typename Polygon<T>::index_type
Polygon<T>::add_face(const triangle_type& tri, const face_descripter& face)
{
    const index_type idx = this->add_face(tri);
    this->faces_.at(idx).face = face;
    return idx;
}

template<typename T>
typename Polygon<T>::index_type
Polygon<T>::connect_edges(const edge_id_type& lhs, const edge_id_type& rhs)
{
    const index_type idx = this->edges_.size();
    const face_id_type fid1 = get_face_id(lhs);
    const face_id_type fid2 = get_face_id(rhs);

    edge_property_type ep;
    ep.faces.first  = std::make_pair(fid1, get_local_index(lhs));
    ep.faces.second = std::make_pair(fid2, get_local_index(rhs));
    ep.vertices.first  = un_initialized;
    ep.vertices.second = un_initialized;
    this->edges_.push_back(ep);

    //XXX update faces
    this->faces_.at(fid1).edges.at(get_local_index(lhs)) = idx;
    this->faces_.at(fid2).edges.at(get_local_index(rhs)) = idx;

    this->faces_.at(fid1).faces.at(get_local_index(lhs)) = fid2;
    this->faces_.at(fid2).faces.at(get_local_index(rhs)) = fid1;

    return idx;
}

template<typename T>
typename Polygon<T>::index_type
Polygon<T>::connect_edges(const edge_id_type& lhs, const edge_id_type& rhs,
                          const edge_descripter& edge)
{
    const index_type idx = this->connect_edges(lhs, rhs);
    this->edges_.at(idx).edge = edge;
    return idx;
}

template<typename T>
typename Polygon<T>::index_type
Polygon<T>::connect_vertices(const std::vector<vertex_id_type>& vtxs)
{
    const index_type idx = this->vertices_.size();
    vertex_property_type vp;
    vp.faces.reserve(vtxs.size());
    vp.edges.reserve(vtxs.size());

    for(typename std::vector<vertex_id_type>::const_iterator
        iter = vtxs.begin(); iter != vtxs.end(); ++iter)
    {
        const face_id_type   fid  = get_face_id(*iter);
        const local_idx_type lidx = get_local_index(*iter);
        const index_type     eidx = this->faces_.at(fid).edges.at(lidx);
        if(eidx == un_initialized)
            throw std::logic_error("edge connection is not specified");

        vp.faces.push_back(std::make_pair(fid, lidx));
        vp.edges.push_back(eidx);

        // XXX: update faces
        this->faces_.at(fid).vertices.at(lidx) = idx;

        // XXX: update edges
        if(this->edges_.at(eidx).vertices.first == un_initialized)
            this->edges_.at(eidx).vertices.first = idx;
        else if(this->edges_.at(eidx).vertices.second == un_initialized)
            this->edges_.at(eidx).vertices.second = idx;
        else
            throw std::invalid_argument(
                    "edges that connects to the vertex is already filled");
    }
    this->vertices_.push_back(vp);

    return idx;
}

template<typename T>
typename Polygon<T>::index_type
Polygon<T>::connect_vertices(const std::vector<vertex_id_type>& vtxs,
                             const vertex_descripter& vertex)
{
    const index_type idx = this->connect_vertices(vtxs);
    this->vertices_.at(idx).vertex = vertex;
    return idx;
}

template<typename T>
std::pair<bool, typename Polygon<T>::edge_id_type>
Polygon<T>::is_connected_by_edge(
        const face_id_type& lhs, const face_id_type& rhs) const
{
    const face_id_type fid1 = get_face_id(lhs);
    const face_id_type fid2 = get_face_id(rhs);

    const face_property_type& fp = this->faces_.at(fid1);
    for(std::size_t i=0; i<3; ++i)
    {
        const std::pair<std::size_t, std::size_t>& fs =
            this->edges_.at(fp.edges[i]).faces;

        if(get_face_id(fs.first) == fid2 || get_face_id(fs.second) == fid2)
            return std::make_pair(true, make_eid(fid1, local_idx_type(i)));
    }
    return std::make_pair(false,
            make_eid(face_id_type(un_initialized), local_idx_type(3)));
}

template<typename T>
std::pair<bool, typename Polygon<T>::vertex_id_type>
Polygon<T>::is_connected_by_vertex(
        const face_id_type& lhs, const face_id_type& rhs) const
{
    const face_id_type fid1 = get_face_id(lhs);
    const face_id_type fid2 = get_face_id(rhs);

    const face_property_type& fp = this->faces_.at(fid1);
    for(std::size_t i=0; i<3; ++i)
    {
        const std::vector<std::pair<std::size_t, std::size_t> >& fs =
            this->vertices_.at(fp.vertices[i]).faces;

        const typename std::vector<std::pair<std::size_t, std::size_t>
            >::const_iterator result = std::find_if(fs.begin(), fs.end(),
                utils::pair_first_element_unary_predicator<
                    std::size_t, std::size_t>(fid2));

        if(result != fs.end())
            return std::make_pair(true, make_vid(fid1, local_idx_type(i)));
        else
            continue;
    }
    return std::make_pair(false,
            make_vid(face_id_type(un_initialized), local_idx_type(3)));
}

template<typename T>
inline std::vector<typename Polygon<T>::face_id_type> const&
Polygon<T>::neighbor_faces(const face_id_type& fid) const
{
    return this->faces_.at(fid).neighbors;
}

template<typename T>
inline boost::array<typename Polygon<T>::face_id_type, 3> const&
Polygon<T>::connecting_faces(const face_id_type& fid) const
{
    return this->faces_.at(fid).faces;
}

template<typename T>
typename Polygon<T>::edge_id_type
Polygon<T>::connecting_edge(const edge_id_type& eid) const
{
    const face_id_type fid = get_face_id(eid);
    const edge_property_type& ep = this->edges_.at(
            this->faces_.at(get_face_id(eid)).edges.at(get_local_index(eid))
            );

    if(face_id_type(ep.faces.first.first) == get_face_id(eid))
        return make_eid(face_id_type(ep.faces.second.first),
                      local_idx_type(ep.faces.second.second));
    else if(face_id_type(ep.faces.second.first) == get_face_id(eid))
        return make_eid(face_id_type(ep.faces.first.first),
                      local_idx_type(ep.faces.first.second));
    else
        throw std::invalid_argument("no connected edge");
}


template<typename T>
Real Polygon<T>::distance_sq(const std::pair<Real3, face_id_type>& lhs,
                             const std::pair<Real3, face_id_type>& rhs) const
{
    throw NotImplemented("ecell4::Polygon::distance_sq");
}

template<typename T>
Real Polygon<T>::distance(const std::pair<Real3, face_id_type>& lhs,
                          const std::pair<Real3, face_id_type>& rhs) const
{
    throw NotImplemented("ecell4::Polygon::distance");
}

template<typename T>
Real3 Polygon<T>::developed_direction(
        const std::pair<Real3, face_id_type>& from,
        const std::pair<Real3, face_id_type>& to) const
{
    throw NotImplemented("ecell4::Polygon::developed_direction");
}

template<typename T>
std::pair<std::vector<typename Polygon<T>::vertex_id_type>,
                      std::pair<typename Polygon<T>::vertex_id_type, Real> >
Polygon<T>::list_vertices_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real radius)
{
    throw NotImplemented("ecell4::Polygon::list_vertices_within_radius");
}

template<typename T>
std::pair<std::vector<typename Polygon<T>::vertex_id_type>,
                      std::pair<typename Polygon<T>::vertex_id_type, Real> >
Polygon<T>::list_faces_within_radius(const Real3& pos, const Real radius)
{
    throw NotImplemented("ecell4::Polygon::list_face_within_radius");
}

template<typename T>
inline typename Polygon<T>::triangle_type&
Polygon<T>::triangle_at(const face_id_type i)
{
    return this->triangles_.at(i);
}

template<typename T>
inline typename Polygon<T>::triangle_type const&
Polygon<T>::triangle_at(const face_id_type i) const
{
    return this->triangles_.at(i);
}

template<typename T>
inline typename Polygon<T>::face_descripter&
Polygon<T>::face_at(const face_id_type i)
{
    return faces_.at(i).face;
}

template<typename T>
inline typename Polygon<T>::face_descripter const&
Polygon<T>::face_at(const face_id_type i) const
{
    return faces_.at(i).face;
}

template<typename T>
inline typename Polygon<T>::vertex_descripter&
Polygon<T>::vertex_at(const vertex_id_type& vid)
{
    return vertex_prop_at(vid).vertex;
}


template<typename T>
inline typename Polygon<T>::vertex_descripter const&
Polygon<T>::vertex_at(const vertex_id_type& vid) const
{
    return vertex_prop_at(vid).vertex;
}

template<typename T>
inline typename Polygon<T>::edge_descripter&
Polygon<T>::edge_at(const edge_id_type& eid)
{
    return edge_prop_at(eid).edge;
}

template<typename T>
inline typename Polygon<T>::edge_descripter const&
Polygon<T>::edge_at(const edge_id_type& eid) const
{
    return edge_prop_at(eid).edge;
}

template<typename T>
inline typename Polygon<T>::vertex_property_type&
Polygon<T>::vertex_prop_at(const vertex_id_type& vid)
{
    return vertices_.at(
        this->faces_.at(get_face_id(vid)).vertices.at(get_local_index(vid)));
}

template<typename T>
inline typename Polygon<T>::vertex_property_type const&
Polygon<T>::vertex_prop_at(const vertex_id_type& vid) const
{
    return vertices_.at(
        this->faces_.at(get_face_id(vid)).vertices.at(get_local_index(vid)));
}

template<typename T>
inline typename Polygon<T>::edge_property_type&
Polygon<T>::edge_prop_at(const edge_id_type& eid)
{
    return edges_.at(
        this->faces_.at(get_face_id(eid)).edges.at(get_local_index(eid)));
}

template<typename T>
inline typename Polygon<T>::edge_property_type const&
Polygon<T>::edge_prop_at(const edge_id_type& eid) const
{
    return edges_.at(
        this->faces_.at(get_face_id(eid)).edges.at(get_local_index(eid)));
}


} // ecell4
#endif// ECELL4_POLYGON
