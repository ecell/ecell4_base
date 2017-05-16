#ifndef ECELL4_POLYGON
#define ECELL4_POLYGON
#include <utility>
#include <functional>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/comparators.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/geometry.hpp>
#include <ecell4/core/Barycentric.hpp>

#include <boost/utility.hpp>
#include <boost/type_traits.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <limits>

namespace ecell4
{

/*! if the arg is already in container, do nothing. otherwise, insert it.
 * @return if inserted, return true.
 */
template<typename containerT, typename argumentT>
inline typename boost::enable_if<
    boost::is_convertible<typename containerT::value_type, argumentT>,
    bool>::type
uniquely_add(containerT& c, const argumentT& arg)
{
    if(std::find(c.begin(), c.end(), arg) != c.end()) return false;
    c.push_back(arg);
    return true;
}

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
    typedef typename traits_type::vertex_id_type    vertex_id_type;
    typedef typename traits_type::edge_id_type      edge_id_type;
    typedef typename traits_type::vertex_descripter vertex_descripter;
    typedef typename traits_type::edge_descripter   edge_descripter;
    typedef typename traits_type::face_descripter   face_descripter;
    typedef typename traits_type::converter_type    converter_type;
    typedef typename traits_type::id_generator_type id_generator_type;
    typedef std::pair<face_id_type, index_type>     local_index_type;
    typedef Barycentric<Real>                       barycentric_type;
    typedef vertex_property<vertex_descripter>      vertex_property_type;
    typedef edge_property<edge_descripter>          edge_property_type;
    typedef face_property<face_descripter>          face_property_type;
    typedef std::vector<vertex_property_type>       vertex_container_type;
    typedef std::vector<edge_property_type>         edge_container_type;
    typedef std::vector<face_property_type>         face_container_type;

  public:

    Polygon(){}
    ~Polygon(){}

    /* set connection --------------------------------------------------------*/
    face_id_type add_face(const triangle_type& tri);
    face_id_type add_face(const triangle_type& tri, const face_descripter& face);

    edge_id_type connect_edges(const local_index_type&, const local_index_type&);
    edge_id_type connect_edges(const local_index_type&, const local_index_type&,
                               const edge_descripter&);

    vertex_id_type connect_vertices(const std::vector<local_index_type>&);
    vertex_id_type connect_vertices(const std::vector<local_index_type>&,
                                    const vertex_descripter&);

    /* get connecting information --------------------------------------------*/
    std::pair<bool, edge_id_type>
    is_connected_by_edge(const face_id_type&, const face_id_type&) const;
    std::pair<bool, vertex_id_type>
    is_connected_by_vertex(const face_id_type&, const face_id_type&) const;

    //! return faceIDs that distance between (points on (it & arg)) is calculatable
    std::vector<face_id_type> const&
    neighbor_faces(const face_id_type& fid) const;
    //! return faceIDs that directory connects through the edges
    boost::array<face_id_type, 3> const&
    adjacent_faces(const face_id_type& fid) const;

    std::pair<face_id_type, face_id_type> const&
    connecting_faces(const edge_id_type& eid) const;
    std::vector<face_id_type> const&
    connecting_faces(const vertex_id_type& vid) const;

    boost::array<vertex_id_type, 3> const&
    connecting_vertices(const face_id_type& fid) const;
    std::pair<vertex_id_type, vertex_id_type> const&
    connecting_vertices(const edge_id_type& eid) const;

    boost::array<edge_id_type, 3> const&
    connecting_edges(const face_id_type& fid) const;
    std::vector<edge_id_type> const&
    connecting_edges(const vertex_id_type& vid) const;

    /* geometric functions ---------------------------------------------------*/
    Real apex_angle(const vertex_id_type& vid) const
    {return this->vertex_prop_at(vid).apex_angle;}

    Real distance_sq(const std::pair<Real3, face_id_type>& lhs,
                     const std::pair<Real3, face_id_type>& rhs) const;
    Real distance(const std::pair<Real3, face_id_type>& lhs,
                  const std::pair<Real3, face_id_type>& rhs) const;

    Real distance_sq(const std::pair<Real3, vertex_id_type>& lhs,
                     const std::pair<Real3, face_id_type>& rhs) const;
    Real distance(const std::pair<Real3, vertex_id_type>& lhs,
                  const std::pair<Real3, face_id_type>& rhs) const
    {return std::sqrt(distance_sq(lhs, rhs));}

    Real distance_sq(const std::pair<Real3, face_id_type>& lhs,
                     const std::pair<Real3, vertex_id_type>& rhs) const
    {return distance_sq(rhs, lhs);}
    Real distance(const std::pair<Real3, face_id_type>& lhs,
                  const std::pair<Real3, vertex_id_type>& rhs) const
    {return distance(rhs, lhs);}

    Real distance_sq(const std::pair<Real3, vertex_id_type>& lhs,
                     const std::pair<Real3, vertex_id_type>& rhs) const;
    Real distance(const std::pair<Real3, vertex_id_type>& lhs,
                  const std::pair<Real3, vertex_id_type>& rhs) const
    {return std::sqrt(distance_sq(lhs, rhs));}

    Real3
    developed_direction(const std::pair<Real3, face_id_type>& from,
                        const std::pair<Real3, face_id_type>& to) const;

    std::vector<std::pair<vertex_id_type, Real> >
    list_vertices_within_radius(const std::pair<Real3, face_id_type>& pos,
                                const Real radius) const;

    std::vector<std::pair<face_id_type, Real> >
    list_faces_within_radius(const Real3& pos, const Real radius) const;

    /* dynamics functions ----------------------------------------------------*/
    std::pair<std::pair<Real3, face_id_type>, Real3>
    move_next_face(const std::pair<Real3, face_id_type>& pos,
                   const Real3& disp) const;

    std::pair<Real3, face_id_type>
    rotate_around_vertex(const std::pair<Real3, face_id_type>& pos,
        const vertex_id_type& apex, const Real r, const Real theta) const;

    /* member access ---------------------------------------------------------*/
    std::size_t num_triangles() const {return faces_.size();}
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

    face_descripter&         at(const face_id_type& i)       {return face_at(i);}
    face_descripter const&   at(const face_id_type& i) const {return face_at(i);}
    vertex_descripter&       at(const vertex_id_type& i)       {return vertex_at(i);}
    vertex_descripter const& at(const vertex_id_type& i) const {return vertex_at(i);}
    edge_descripter&         at(const edge_id_type& i)       {return edge_at(i);}
    edge_descripter const&   at(const edge_id_type& i) const {return edge_at(i);}

    edge_id_type   get_edge_id(const local_index_type idx) const;
    vertex_id_type get_vertex_id(const local_index_type idx) const;

    std::vector<face_id_type>   list_face_id()   const;
    std::vector<vertex_id_type> list_vertex_id() const;
    std::vector<edge_id_type>   list_edge_id()   const;

    const converter_type& converter() const {return converter_;}

    /* required by shape -----------------------------------------------------*/
    dimension_kind dimension() const {return THREE;} // TWO?

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

//     Real distance_over_edge(const std::pair<Real3, face_id_type>& lhs,
//                             const std::pair<Real3, face_id_type>& rhs) const;
//     Real distance_over_vertex(const std::pair<Real3, face_id_type>& lhs,
//                               const std::pair<Real3, face_id_type>& rhs) const;
  private:

    face_id_type   generate_face_id()   {return idgen_.generate_face_id();}
    edge_id_type   generate_edge_id()   {return idgen_.generate_edge_id();}
    vertex_id_type generate_vertex_id() {return idgen_.generate_vertex_id();}

    template<typename T_id>
    index_type to_index(const T_id& id) const {return converter_.to_index(id);}
    template<typename T_id>
    T_id to_id(const index_type& i) const {return converter_.template to_id<T_id>(i);}
    template<typename T_id>
    void link(const T_id& id, index_type idx) {return converter_.link(id, idx);}

    template<typename Tid>
    static inline Tid un_initialized()
    {
        return traits_type::template un_initialized<Tid>();
    }

  private:

    template<typename vertexT>
    struct vertex_property
    {
        vertex_id_type id;
        vertexT vertex;
        Real apex_angle;
        std::vector<edge_id_type> edges;
        std::vector<face_id_type> faces;
        std::vector<index_type>   local_indices; //XXX: faces.vertex_at(this);
    };

    template<typename edgeT>
    struct edge_property
    {
        edge_id_type id;
        edgeT edge;
        Real tilt_angle;
        std::pair<vertex_id_type, vertex_id_type> vertices;
        std::pair<face_id_type,   face_id_type>   faces;
        std::pair<index_type,     index_type>     local_indices;
    };

    template<typename faceT>
    struct face_property
    {
        face_id_type  id;
        faceT         face;     // additional information for face
        triangle_type triangle; // geometric  information
        boost::array<vertex_id_type, 3> vertices;
        boost::array<edge_id_type, 3>   edges;
        boost::array<face_id_type, 3>   adjacents; // connected by edge
        std::vector<face_id_type>       neighbors; // connected by edge of vtx
    };

    converter_type        converter_;
    id_generator_type     idgen_;
    vertex_container_type vertices_;
    edge_container_type   edges_;
    face_container_type   faces_;
};

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

/* set connection ------------------------------------------------------------*/

template<typename T>
typename Polygon<T>::face_id_type
Polygon<T>::add_face(const triangle_type& triangle)
{
    const face_id_type fid = generate_face_id();

    const index_type idx = this->faces_.size();
    this->link(fid, idx);

    face_property_type fp;
    fp.id = fid;
    fp.triangle = triangle;
    fp.vertices.fill(un_initialized<vertex_id_type>());
    fp.edges.fill(un_initialized<edge_id_type>());
    fp.adjacents.fill(un_initialized<face_id_type>());
    this->faces_.push_back(fp);

    return fid;
}

template<typename T>
typename Polygon<T>::face_id_type
Polygon<T>::add_face(const triangle_type& tri, const face_descripter& face)
{
    const face_id_type fid = this->add_face(tri);
    this->face_at(fid) = face;
    return fid;
}

template<typename T>
typename Polygon<T>::edge_id_type
Polygon<T>::connect_edges(const local_index_type& lhs, const local_index_type& rhs)
{
    const edge_id_type eid  = generate_edge_id();

    const index_type idx = this->edges_.size();
    this->link(eid, idx);

    const face_id_type fid1 = lhs.first;
    const face_id_type fid2 = rhs.first;

    edge_property_type ep;
    ep.id            = eid;
    ep.faces         = std::make_pair(fid1, fid2);
    ep.local_indices = std::make_pair(lhs.second, rhs.second);
    ep.vertices      = std::make_pair(un_initialized<vertex_id_type>(),
                                      un_initialized<vertex_id_type>());
    ep.tilt_angle    = angle(this->triangle_at(fid1).normal(),
                             this->triangle_at(fid2).normal());
    this->edges_.push_back(ep);

    //XXX update faces
    this->face_prop_at(fid1).edges.at(lhs.second) = eid;
    this->face_prop_at(fid2).edges.at(rhs.second) = eid;
    this->face_prop_at(fid1).adjacents.at(lhs.second) = fid2;
    this->face_prop_at(fid2).adjacents.at(rhs.second) = fid1;
    uniquely_add(this->face_prop_at(fid1).neighbors, fid2);
    uniquely_add(this->face_prop_at(fid2).neighbors, fid1);
    return eid;
}

template<typename T>
typename Polygon<T>::edge_id_type
Polygon<T>::connect_edges(const local_index_type& lhs,
        const local_index_type& rhs, const edge_descripter& edge)
{
    const edge_id_type eid = this->connect_edges(lhs, rhs);
    this->edge_at(eid) = edge;
    return eid;
}

template<typename T>
typename Polygon<T>::vertex_id_type
Polygon<T>::connect_vertices(const std::vector<local_index_type>& vtxs)
{
    const vertex_id_type vid = generate_vertex_id();
    const index_type idx = this->vertices_.size();
    this->link(vid, idx);

    vertex_property_type vp;
    vp.id = vid;
    vp.faces.reserve(vtxs.size());
    vp.edges.reserve(vtxs.size());
    vp.local_indices.reserve(vtxs.size());
    vp.apex_angle = 0.;

    for(typename std::vector<local_index_type>::const_iterator
        iter = vtxs.begin(); iter != vtxs.end(); ++iter)
    {
        const face_id_type fid  = iter->first;
        const index_type   lidx = iter->second;
        const edge_id_type   eid = this->face_prop_at(fid).edges.at(lidx);
        if(eid == un_initialized<edge_id_type>())
            throw std::logic_error("edge connection is not specified");

        vp.faces.push_back(fid);
        vp.local_indices.push_back(lidx);
        vp.edges.push_back(eid);
        vp.apex_angle += this->triangle_at(fid).angle_at(lidx);

        // XXX: update faces
        this->face_prop_at(fid).vertices.at(lidx) = vid;

        // XXX: update edges
        if(this->edge_prop_at(eid).vertices.first ==
           un_initialized<vertex_id_type>())
        {
            this->edge_prop_at(eid).vertices.first = vid;
        }
        else if(this->edge_prop_at(eid).vertices.second ==
                un_initialized<vertex_id_type>())
        {
            this->edge_prop_at(eid).vertices.second = vid;
        }
        else
        {
            throw std::invalid_argument(
                    "edges that connects to the vertex is already filled");
        }
    }

    // XXX: update face_prop.neighbors
    for(typename std::vector<face_id_type>::const_iterator
        iter = vp.faces.begin(); iter != vp.faces.end(); ++iter)
    {
        face_property_type& fp = this->face_prop_at(*iter);
        for(typename std::vector<face_id_type>::const_iterator
            jter = vp.faces.begin(); jter != vp.faces.end(); ++jter)
        {
            if(*iter == *jter) continue;
            uniquely_add(fp.neighbors, *jter);
        }
    }

    this->vertices_.push_back(vp);

    return vid;
}

template<typename T>
typename Polygon<T>::vertex_id_type
Polygon<T>::connect_vertices(const std::vector<local_index_type>& vtxs,
                             const vertex_descripter& vertex)
{
    const vertex_id_type vid = this->connect_vertices(vtxs);
    this->vertex_at(vid) = vertex;
    return vid;
}

/* get connecting information ------------------------------------------------*/

template<typename T>
std::pair<bool, typename Polygon<T>::edge_id_type>
Polygon<T>::is_connected_by_edge(
        const face_id_type& fid1, const face_id_type& fid2) const
{
    const face_property_type& fp = this->face_prop_at(fid1);
    for(std::size_t i=0; i<3; ++i)
    {
        if(fp.adjacents[i] == fid2) return std::make_pair(true, fp.edges[i]);
    }
    return std::make_pair(false, un_initialized<edge_id_type>());
}

template<typename T>
std::pair<bool, typename Polygon<T>::vertex_id_type>
Polygon<T>::is_connected_by_vertex(
        const face_id_type& fid1, const face_id_type& fid2) const
{
    const face_property_type& fp = this->face_prop_at(fid1);
    for(std::size_t i=0; i<3; ++i)
    {
        const std::vector<face_id_type>& fs =
            this->vertices_.at(fp.vertices[i]).faces;

        if(std::find(fs.begin(), fs.end(), fid2) != fs.end())
            return std::make_pair(true, fp.vertices[i]);
    }
    return std::make_pair(false, un_initialized<vertex_id_type>());
}

template<typename T>
inline std::vector<typename Polygon<T>::face_id_type> const&
Polygon<T>::neighbor_faces(const face_id_type& fid) const
{
    return this->face_prop_at(fid).neighbors;
}

template<typename T>
inline boost::array<typename Polygon<T>::face_id_type, 3> const&
Polygon<T>::adjacent_faces(const face_id_type& fid) const
{
    return this->face_prop_at(fid).adjacents;
}

template<typename T>
inline std::pair<typename Polygon<T>::face_id_type,
                 typename Polygon<T>::face_id_type> const&
Polygon<T>::connecting_faces(const edge_id_type& eid) const
{
    return this->edge_prop_at(eid).faces;
}

template<typename T>
inline std::vector<typename Polygon<T>::face_id_type> const&
Polygon<T>::connecting_faces(const vertex_id_type& vid) const
{
    return this->vertex_prop_at(vid).faces;
}

template<typename T>
inline std::pair<typename Polygon<T>::vertex_id_type,
                 typename Polygon<T>::vertex_id_type> const&
Polygon<T>::connecting_vertices(const edge_id_type& eid) const
{
    return this->edge_prop_at(eid).vertices;
}

template<typename T>
inline boost::array<typename Polygon<T>::vertex_id_type, 3> const&
Polygon<T>::connecting_vertices(const face_id_type& fid) const
{
    return this->face_prop_at(fid).vertices;
}

template<typename T>
inline std::vector<typename Polygon<T>::edge_id_type> const&
Polygon<T>::connecting_edges(const vertex_id_type& vid) const
{
    return this->vertex_prop_at(vid).edges;
}

template<typename T>
inline boost::array<typename Polygon<T>::edge_id_type, 3> const&
Polygon<T>::connecting_edges(const face_id_type& fid) const
{
    return this->face_prop_at(fid).edges;
}


/* geometric functions -------------------------------------------------------*/

template<typename T>
Real Polygon<T>::distance_sq(const std::pair<Real3, face_id_type>& lhs,
                             const std::pair<Real3, face_id_type>& rhs) const
{
    if(lhs.second == rhs.second)
        return length_sq(lhs.first - rhs.first);

    const std::pair<bool, edge_id_type> edg =
        this->is_connected_by_edge(lhs.second, rhs.second);
    if(edg.first)
    {
        const edge_property_type& edge = this->edge_prop_at(edg.second);
        index_type lidx;
        if(edge.faces.first == lhs.second)
        {
            lidx = edge.local_indices.first;
        }
        else if(edge.faces.second == lhs.second)
        {
            lidx = edge.local_indices.second;
        }
        else
        {
            throw std::logic_error(
                    "Polygon::distance_sq: edge is not connected!");
        }

        const triangle_type&     lhs_t = this->triangle_at(lhs.second);
        const Real3 developped = lhs_t.vertex_at(lidx) +
            rotate(-1. * edge.tilt_angle,
                   lhs_t.edge_at(lidx) / lhs_t.length_of_edge_at(lidx),
                   rhs.first - lhs_t.vertex_at(lidx));

        return length_sq(lhs.first - developped);
    }

    const std::pair<bool, vertex_id_type> is_c_vtx =
        this->is_connected_by_vertex(lhs.second, rhs.second);

    if(is_c_vtx.first)
    {
        const vertex_property_type& vtx = this->vertex_prop_at(is_c_vtx.second);
        const std::vector<face_id_type>& fs = vtx.faces;
        const typename std::vector<face_id_type>::const_iterator result =
            std::find(fs.begin(), fs.end(), lhs.second);
        index_type idx = std::distance(fs.begin(), result);
        const index_type lidx = vtx.local_indices.at(idx);

        const Real3     vtx_position = triangle_at(lhs.second).vertex_at(lidx);
        const Real3       lhs_to_vtx = vtx_position - lhs.first;
        const Real3       vtx_to_rhs = rhs.first - vtx_position;
        const Real  lhs_to_vtx_lensq = length_sq(lhs_to_vtx);
        const Real  rhs_to_vtx_lensq = length_sq(vtx_to_rhs);
        const Real        apex_angle = vertex_prop_at(is_c_vtx.second).apex_angle;

        Real inter_angle = angle(lhs_to_vtx, triangle_at(lhs.second).edge_at(
                           (lidx == 0) ? 2 : lidx-1));

        // XXX: order of face idx
        bool round = false;
        ++idx;
        if(idx == fs.size())
        {
            idx = 0;
            round = true;
        }

        while(true)
        {
            const face_id_type fid(fs.at(idx));
            const index_type   lvidx(vtx.local_indices.at(idx));
            const triangle_type& f = this->triangle_at(fid);
            if(fid == rhs.second)
            {
                inter_angle += angle(vtx_to_rhs, f.edge_at(lvidx));
                break;
            }
            inter_angle += f.angle_at(lvidx);

            ++idx;
            if(idx == fs.size())
            {
                if(round)
                    throw std::logic_error("Polygon::distance: rhs not found");
                idx = 0;
                round = true;
            }
        }
        assert(inter_angle <= apex_angle);

        const Real min_angle = std::min(inter_angle, apex_angle - inter_angle);
        return lhs_to_vtx_lensq + rhs_to_vtx_lensq - 2. *
               std::sqrt(lhs_to_vtx_lensq * rhs_to_vtx_lensq) * std::cos(min_angle);
    }

    // lhs and rhs don't share edge nor vertex
    return std::numeric_limits<Real>::infinity();
}

template<typename T>
Real Polygon<T>::distance(const std::pair<Real3, face_id_type>& lhs,
                          const std::pair<Real3, face_id_type>& rhs) const
{
    return std::sqrt(this->distance_sq(lhs, rhs));
}

template<typename T>
Real Polygon<T>::distance_sq(const std::pair<Real3, vertex_id_type>& lhs,
                             const std::pair<Real3, face_id_type>& rhs) const
{
    const face_property_type& fp = face_prop_at(rhs.second);
    for(std::size_t i=0; i<3; ++i)
        if(fp.vertices[i] == lhs.second)
            return length_sq(lhs.first - rhs.first);

    //XXX it is enough for sgfrd,
    //    but basically one can measure the distance from more distant vertex
    return std::numeric_limits<Real>::max();
}

template<typename T>
Real Polygon<T>::distance_sq(const std::pair<Real3, vertex_id_type>& lhs,
                             const std::pair<Real3, vertex_id_type>& rhs) const
{
    const vertex_property_type& vp = vertex_prop_at(lhs.second);
    for(typename std::vector<edge_id_type>::const_iterator
        iter = vp.edges.begin(); iter != vp.edges.end(); ++iter)
    {
        const edge_property_type& ep = edge_prop_at(*iter);
        if(ep.vertices.first  == rhs.second || ep.vertices.second == rhs.second)
            return length_sq(lhs.first - rhs.first);
    }
    //XXX it is enough for sgfrd,
    //    but basically one can measure the distance from more distant vertex
    return std::numeric_limits<Real>::max();
}

template<typename T>
Real3 Polygon<T>::developed_direction(
        const std::pair<Real3, face_id_type>& lhs,
        const std::pair<Real3, face_id_type>& rhs) const
{
    if(lhs.second == rhs.second)
        return rhs.first - lhs.first;

    const std::pair<bool, edge_id_type> edg =
        this->is_connected_by_edge(lhs.second, rhs.second);
    if(edg.first)
    {
        const edge_property_type& edge = this->edge_prop_at(edg.second);
        index_type lidx = std::numeric_limits<index_type>::max();
        if(edge.faces.first == lhs.second)
        {
            lidx = edge.local_indices.first;
        }
        else if(edge.faces.second == lhs.second)
        {
            lidx = edge.local_indices.second;
        }
        else
        {
            throw std::logic_error(
                    "Polygon::distance_sq: edge is not connected!");
        }

        const triangle_type&     lhs_t = this->triangle_at(lhs.second);
        const Real3 developped = lhs_t.vertex_at(lidx) +
            rotate(-1. * edge.tilt_angle,
                   lhs_t.edge_at(lidx) / lhs_t.length_of_edge_at(lidx),
                   rhs.first - lhs_t.vertex_at(lidx));

        return developped - lhs.first;
    }

    const std::pair<bool, vertex_id_type> is_vtx =
        this->is_connected_by_vertex(lhs.second, rhs.second);

    if(is_vtx.first)
    {
        const vertex_property_type& vtx = this->vertex_prop_at(is_vtx.second);
        const std::vector<face_id_type>& fs = vtx.faces;
        const typename std::vector<face_id_type>::const_iterator result =
            std::find(fs.begin(), fs.end(), lhs.second);
        index_type idx = std::distance(fs.begin(), result);
        const index_type lidx = vtx.local_indices.at(idx);

        const Real3     vtx_position = triangle_at(lhs.second).vertex_at(lidx);
        const Real3       lhs_to_vtx = vtx_position - lhs.first;
        const Real3       vtx_to_rhs = rhs.first - vtx_position;
        const Real3           normal = this->triangle_at(lhs.second).normal();
        const Real  lhs_to_vtx_lensq = length_sq(lhs_to_vtx);
        const Real  rhs_to_vtx_lensq = length_sq(vtx_to_rhs);
        const Real        apex_angle = vertex_prop_at(is_vtx.second).apex_angle;

        Real inter_angle = angle(lhs_to_vtx, triangle_at(lhs.second).edge_at(
                           (lidx == 0) ? 2 : lidx-1));

        // XXX: order of face idx
        bool round = false;
        ++idx;
        if(idx == fs.size())
        {
            idx = 0;
            round = true;
        }

        while(true)
        {
            const face_id_type fid(fs.at(idx));
            const index_type   lvidx(vtx.local_indices.at(idx));
            const triangle_type& f = this->triangle_at(fid);
            if(fid == rhs.second)
            {
                inter_angle += angle(vtx_to_rhs, f.edge_at(lvidx));
                break;
            }
            inter_angle += f.angle_at(lvidx);

            ++idx;
            if(idx == fs.size())
            {
                if(round)
                    throw std::logic_error("Polygon::distance: rhs not found");
                idx = 0;
                round = true;
            }
        }
        assert(inter_angle <= apex_angle);

        const Real min_angle = std::min(inter_angle, apex_angle - inter_angle);

        return lhs_to_vtx + rotate(min_angle, normal, lhs_to_vtx * (-1.0)) *
                         std::sqrt(rhs_to_vtx_lensq / lhs_to_vtx_lensq);
    }
    throw NotSupported("Polygon::developped_direction couldn't determine the path");
}

template<typename T>
std::vector<std::pair<typename Polygon<T>::vertex_id_type, Real> >
Polygon<T>::list_vertices_within_radius(
        const std::pair<Real3, face_id_type>& pos, const Real radius) const
{
    // XXX: it searchs vertices only on the adjacent faces
    // because gfrd circular shell size is restricted so that the shell does not
    // contain the incenter of neighbor face, nomally this range is enough.

    std::vector<std::pair<vertex_id_type, Real> > list; list.reserve(6);

    const Real rad2 = radius * radius;
    triangle_type const& tri = triangle_at(pos.second);
    for(std::size_t i=0; i<3; ++i)
    {
        const Real dist2 = length_sq(pos.first - tri.vertex_at(i));
        if(dist2 <= rad2)
        {
            list.push_back(std::make_pair(
                face_prop_at(pos.second).vertices.at(i), std::sqrt(dist2)));
        }
    }

    boost::array<vertex_id_type, 3> const& known = face_prop_at(pos.second).vertices;
    boost::array<face_id_type, 3> const& adjs = face_prop_at(pos.second).adjacents;
    for(std::size_t i=0; i<3; ++i)
    {
        face_property_type const& f = face_prop_at(adjs[i]);
        const typename boost::array<vertex_id_type, 3>::const_iterator unknown =
            std::find_if(f.vertices.begin(), f.vertices.end(),
                         boost::lambda::_1 != known[0] &&
                         boost::lambda::_1 != known[1] &&
                         boost::lambda::_1 != known[2]);
        const std::size_t idx = std::distance(f.vertices.begin(), unknown);
        assert(idx != 3);
        const Real dist2 = this->distance_sq(
                pos, std::make_pair(f.triangle.vertex_at(idx), adjs[i]));

        if(dist2 <= rad2)
        {
            list.push_back(std::make_pair(*unknown, std::sqrt(dist2)));
        }
    }

    std::sort(list.begin(), list.end(), utils::pair_second_element_comparator<
              vertex_id_type, Real>());
    return list;
}

template<typename T>
std::vector<std::pair<typename Polygon<T>::face_id_type, Real> >
Polygon<T>::list_faces_within_radius(const Real3& pos, const Real radius) const
{
    // TODO: smarter (than brute force) spatial partition method is needed.
    std::vector<std::pair<face_id_type, Real> > list;

    const Real rad2 = radius * radius;
    for(typename face_container_type::const_iterator
            iter(faces_.begin()), end(faces_.end()); iter != end; ++iter)
    {
        const Real dist2 = distance_sq(pos, iter->triangle);
        if(dist2 < rad2)
        {
            list.push_back(std::make_pair(iter->id, std::sqrt(dist2)));
        }
    }
    std::sort(list.begin(), list.end(), utils::pair_second_element_comparator<
              face_id_type, Real>());
    return list;
}

template<typename T>
std::pair<std::pair<Real3, typename Polygon<T>::face_id_type>, Real3>
Polygon<T>::move_next_face(const std::pair<Real3, face_id_type>& pos,
                           const Real3& disp) const
{
    const triangle_type& t = this->triangle_at(pos.second);
    const barycentric_type newpos = to_barycentric(pos.first + disp, t);

    if(::ecell4::is_inside(newpos))
        return std::make_pair(
            std::make_pair(to_absolute(newpos, t), pos.second), Real3(0.,0.,0.));

    // calculate partial displacement to the edge
    const barycentric_type bpos = to_barycentric(pos.first, t);
    const barycentric_type bdis = newpos - bpos;

    const std::pair<std::size_t, Real> cross = first_cross_edge(bpos, bdis);
    const edge_id_type next_edge =
        this->face_prop_at(pos.second).edges.at(cross.first);

    const face_id_type next_face_id =
        this->face_prop_at(pos.second).adjacents.at(cross.first);
    const Real3 next_pos  = pos.first + disp * cross.second;

    // turn displacement
    const Real  ang  = edge_prop_at(next_edge).tilt_angle;
    const Real3 axis = t.edge_at(cross.first) / t.length_of_edge_at(cross.first);
    const Real3 next_disp = rotate(ang, axis, disp) * (1. - cross.second);

    return std::make_pair(std::make_pair(next_pos, next_face_id), next_disp);
}

template<typename T>
std::pair<Real3, typename Polygon<T>::face_id_type>
Polygon<T>::rotate_around_vertex(const std::pair<Real3, face_id_type>& pos,
    const vertex_id_type& apex, const Real r, const Real theta) const
{
    assert(theta <= this->apex_angle(apex));


    const vertex_property_type&     vtx = this->vertex_prop_at(apex);
    const std::vector<face_id_type>& fs = vtx.faces;
    const typename std::vector<face_id_type>::const_iterator initial =
        std::find(fs.begin(), fs.end(), pos.second);
    const std::size_t       idx = std::distance(fs.begin(), initial);
    const std::size_t      lidx = vtx.local_indices.at(idx);
    const triangle_type& initri = this->triangle_at(pos.second);
    const Real3    vtx_position = initri.vertex_at(lidx);
    const Real      inter_angle = angle(vtx_position - pos.first,
            initri.edge_at((lidx==0)?2:lidx-1));

    if(inter_angle >= theta)
    {
        const Real3 direction = rotate(theta, initri.normal(),
                                       pos.first - vtx_position);
        return std::make_pair(
            vtx_position + direction * (r / length(direction)), pos.second);
    }

    Real rest_angle = theta - inter_angle;
    typename std::vector<face_id_type>::const_iterator iter = initial + 1;
    while(iter != initial)
    {
        if(iter == fs.end()) iter = fs.begin();

        const std::size_t  i = std::distance(fs.begin(), iter);
        const std::size_t vi = vtx.local_indices.at(i);
        const triangle_type& tri = this->triangle_at(*iter);
        if(rest_angle <= tri.angle_at(vi))
        {
            const Real3 direction = rotate(rest_angle, tri.normal(),
                                           tri.edge_at(vi));
            return std::make_pair(
                vtx_position + direction * (r / length(direction)), *iter);
        }
        rest_angle -= tri.angle_at(vi);
        ++iter;
    }

    assert(rest_angle >= 0.);

    const Real3 direction = rotate(rest_angle, initri.normal(),
                                   initri.edge_at(lidx));
    return std::make_pair(
        vtx_position + direction * (r / length_sq(direction)), pos.second);

}

/* member access -------------------------------------------------------------*/

template<typename T>
inline typename Polygon<T>::triangle_type&
Polygon<T>::triangle_at(const face_id_type i)
{
    return face_prop_at(i).triangle;
}

template<typename T>
inline typename Polygon<T>::triangle_type const&
Polygon<T>::triangle_at(const face_id_type i) const
{
    return face_prop_at(i).triangle;
}

template<typename T>
inline typename Polygon<T>::face_descripter&
Polygon<T>::face_at(const face_id_type i)
{
    return face_prop_at(i).face;
}

template<typename T>
inline typename Polygon<T>::face_descripter const&
Polygon<T>::face_at(const face_id_type i) const
{
    return face_prop_at(i).face;
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
    return vertices_.at(to_index(vid));
}

template<typename T>
inline typename Polygon<T>::vertex_property_type const&
Polygon<T>::vertex_prop_at(const vertex_id_type& vid) const
{
    return vertices_.at(to_index(vid));
}

template<typename T>
inline typename Polygon<T>::edge_property_type&
Polygon<T>::edge_prop_at(const edge_id_type& eid)
{
    return edges_.at(to_index(eid));
}

template<typename T>
inline typename Polygon<T>::edge_property_type const&
Polygon<T>::edge_prop_at(const edge_id_type& eid) const
{
    return edges_.at(to_index(eid));
}

template<typename T>
inline typename Polygon<T>::face_property_type&
Polygon<T>::face_prop_at(const face_id_type& fid)
{
    return faces_.at(to_index(fid));
}

template<typename T>
inline typename Polygon<T>::face_property_type const&
Polygon<T>::face_prop_at(const face_id_type& fid) const
{
    return faces_.at(to_index(fid));
}

template<typename T>
inline typename Polygon<T>::edge_id_type
Polygon<T>::get_edge_id(const local_index_type idx) const
{
    return this->face_prop_at(idx.first).edges.at(idx.second);
}

template<typename T>
inline typename Polygon<T>::vertex_id_type
Polygon<T>::get_vertex_id(const local_index_type idx) const
{
    return this->face_prop_at(idx.first).vertices.at(idx.second);
}

template<typename T>
std::vector<typename Polygon<T>::face_id_type>
Polygon<T>::list_face_id()   const
{
    std::vector<face_id_type> retval; retval.reserve(faces_.size());
    for(typename face_container_type::const_iterator
        iter = faces_.begin(); iter != faces_.end(); ++iter)
    {
        retval.push_back(iter->id);
    }
    return retval;
}

template<typename T>
std::vector<typename Polygon<T>::vertex_id_type>
Polygon<T>::list_vertex_id() const
{
    std::vector<vertex_id_type> retval; retval.reserve(vertices_.size());
    for(typename vertex_container_type::const_iterator
        iter = vertices_.begin(); iter != vertices_.end(); ++iter)
    {
        retval.push_back(iter->id);
    }
    return retval;
}
template<typename T>
std::vector<typename Polygon<T>::edge_id_type>
Polygon<T>::list_edge_id() const
{
    std::vector<edge_id_type> retval; retval.reserve(vertices_.size());
    for(typename edge_container_type::const_iterator
        iter = vertices_.begin(); iter != vertices_.end(); ++iter)
    {
        retval.push_back(iter->id);
    }
    return retval;
}


} // ecell4
#endif// ECELL4_POLYGON
