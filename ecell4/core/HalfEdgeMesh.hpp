#ifndef ECELL4_POLYGON
#define ECELL4_POLYGON

#include <utility>
#include <functional>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/comparators.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/geometry.hpp>
#include <ecell4/core/triangle_geometry.hpp>
#include <ecell4/core/Barycentric.hpp>
#include <ecell4/core/IndexedContainer.hpp>

#include <boost/utility.hpp>
#include <boost/type_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/array.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/cstdint.hpp>

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <limits>

namespace ecell4
{

// Static Polygon. once made, the shape never change.
// (Face|Edge|Vertex)ID are just std::size_t.
class Polygon : public Shape
{
  public:

    static const Real absolute_tolerance;
    static const Real relative_tolerance;

    typedef std::size_t   face_id_type;
    typedef std::size_t   edge_id_type;
    typedef std::size_t vertex_id_type;

    struct vertex_data
    {
        bool  is_contiguous; // there are faces successively around the vertex
        Real  angle;         // total angle around this vertex
        Real3 position;      // position of this vertex
        std::vector<edge_id_type> outgoing_edges;
    };
    struct edge_data
    {
        Real           length;    // length of this edge
        Real           tilt;      // tilt from this->face to opposite face
        Real3          direction; // direction independent from boundary
        vertex_id_type target;    // destination vertex
        face_id_type   face;      // belonging face
        edge_id_type   next;      // edge on the same face, starting from target
        boost::optional<edge_id_type> opposite_edge;
    };
    struct face_data
    {
        Triangle triangle; // not considering Boundary, contains just the shape
        std::array<  edge_id_type, 3> edges;    // idx consistent with triangle
        std::array<vertex_id_type, 3> vertices; // idx consistent with triangle
    };

    typedef std::vector<vertex_data> vertex_container_type;
    typedef std::vector<  face_data>   face_container_type;
    typedef std::vector<  edge_data>   edge_container_type;

    // neighbor list (of faces) implementation
    typedef std::vector<std::size_t>                          nlist_idxs_type;
    typedef std::vector<std::pair<std::size_t, std::size_t> > nlist_ranges_type;

  public:

    Polygon(const Real3& edge_length)
        : edge_length_(edge_length), total_area(0.0),
    {}
    Polygon(const Real3& edge_length, const std::vector<Triangle>& ts)
        : edge_length_(edge_length), total_area(0.0),
    {
        this->assign(ts, tol_abs, tol_rel);
    }
    ~Polygon(){}

    // clear current polygon and assign different one.
    // tolerances are used to detect the same vertices in different triangles.
    void assign(const std::vector<Triangle>& ts);

    std::pair<Real3, face_id_type>
    move_on_surface(const std::pair<Real3, face_id_type>& pos,
        const Real3& disp, std::size_t& max_edges_traversed) const;

    std::pair<Real3, face_id_type>
    roll_around_vertex(const std::pair<Real3, face_id_type>& pos,
        const Real r, const Real theta, const vertex_id_type& vid) const;

    // pos1 -> pos2 <=> pos2 - pos1
    Real3 direction(const std::pair<Real3,   face_id_type>& pos1,
                    const std::pair<Real3,   face_id_type>& pos2) const;
    Real3 direction(const std::pair<Real3, vertex_id_type>& pos1,
                    const std::pair<Real3, vertex_id_type>& pos2) const;
    Real3 direction(const std::pair<Real3,   face_id_type>& pos1,
                    const std::pair<Real3, vertex_id_type>& pos2) const;
    Real3 direction(const std::pair<Real3, vertex_id_type>& pos1,
                    const std::pair<Real3,   face_id_type>& pos2) const;

    Real distance_sq(const std::pair<Real3,   face_id_type>& pos1,
                     const std::pair<Real3,   face_id_type>& pos2) const;
    Real distance_sq(const std::pair<Real3, vertex_id_type>& pos1,
                     const std::pair<Real3, vertex_id_type>& pos2) const;
    Real distance_sq(const std::pair<Real3,   face_id_type>& pos1,
                     const std::pair<Real3, vertex_id_type>& pos2) const;
    Real distance_sq(const std::pair<Real3, vertex_id_type>& pos1,
                     const std::pair<Real3,   face_id_type>& pos2) const;

    Real distance(const std::pair<Real3,   face_id_type>& pos1,
                  const std::pair<Real3,   face_id_type>& pos2) const
    {return this->std::sqrt(distance_sq(pos1, pos2));}
    Real distance(const std::pair<Real3, vertex_id_type>& pos1,
                  const std::pair<Real3, vertex_id_type>& pos2) const
    {return this->std::sqrt(distance_sq(pos1, pos2));}
    Real distance(const std::pair<Real3,   face_id_type>& pos1,
                  const std::pair<Real3, vertex_id_type>& pos2) const
    {return this->std::sqrt(distance_sq(pos1, pos2));}
    Real distance(const std::pair<Real3, vertex_id_type>& pos1,
                  const std::pair<Real3,   face_id_type>& pos2) const
    {return this->std::sqrt(distance_sq(pos1, pos2));}

    bool is_contiguous(const vertex_id_type& vid) const
    {return this->vertices_.at(vid).is_contiguous;}

    boost::optional<vertex_id_type>
    find_vertex(const Real& pos) const
    {
        const Real tol_rel2 = relative_tolerance_rel * relative_tolerance;
        const Real tol_abs2 = absolute_tolerance_abs * absolute_tolerance;

        for(std::size_t i=0; i<vertices_.size(); ++i)
        {
            const vertex_data& vd = vertices_[i];
            const Real dist_sq = length_sq(
                this->periodic_transpose(vd.position, pos) - pos);

            if(dist_sq < tol_abs2 || dist_sq < length_sq(pos) * tol_rel2)
            {
                return i;
            }
        }
        return boost::none;
    }

    // half-edge traverse
    // next edge: the edge belonging the same face,
    //            starting from the target of current edge
    edge_id_type
    next_of(const edge_id_type eid) const
    {return this->edge_at(eid).next;}

    boost::optional<edge_id_type>
    next_of(const boost::optional<edge_id_type>& eid) const
    {return (eid) ? (this->edge_at(*eid).next) : (boost::none);}

    // opposite edge: the edge that starts from the target of current edge,
    //                ends at the starting point of current edge.
    boost::optional<edge_id_type>
    opposite_of(const edge_id_type eid) const
    {return this->edge_at(eid).opposite_edge;}

    boost::optional<edge_id_type>
    opposite_of(const boost::optional<edge_id_type>& eid) const
    {return (eid) ? (this->edge_at(*eid).opposite_edge) : (boost::none);}

    // target vertex: the vertex that current edge stops at.
    vertex_id_type
    target_of(const edge_id_type eid) const
    {return this->edge_at(eid).target;}

    boost::optional<vertex_id_type>
    target_of(const boost::optional<edge_id_type>& eid) const
    {return (eid) ? (this->edge_at(*eid).target) : (boost::none);}

    // belonging face: the face that corresponds to the current edge.
    face_id_type
    face_of(const edge_id_type eid) const
    {return this->edge_at(eid).face;}

    boost::optional<face_id_type>
    face_of(const boost::optional<edge_id_type>& eid) const
    {return (eid) ? (this->edge_at(*eid).face) : (boost::none);}

    // edge ids that are corresponds to the face
    boost::array<edge_id_type, 3> const&
    edges_of(const face_id_type fid) const
    {return this->face_at(fid).edges;}

    // vertex ids that are corresponds to the face
    boost::array<vertex_id_type, 3> const&
    vertices_of(const face_id_type fid) const
    {return this->face_at(fid).vertices;}

    // edge ids that starts from the vertex
    std::vector<edge_id_type> const&
    outgoing_edges(const vertex_id_type vid) const
    {return this->vertices_.at(vid).outgoing_edges;}

    // accessor ---------------------------------------------------------------

    Real apex_angle_at(const vertex_id_type& vid) const
    {
        return this->vertices_.at(vid).angle;
    }
    Real position_at(const vertex_id_type& vid) const
    {
        return this->vertices_.at(vid).position;
    }
    std::vector<face_id_type>
    connecting_faces(const vertex_id_type& vid) const
    {
        std::vector<face_id_type> retval(outs.size());

        const std::vector<edge_id_type>& outs = this->outgoing_edges(vid);
        for(typename std::vector<edge_id_type>::const_iterator
            i(outs.begin()), e(outs.end()); i!=e; ++i)
        {
            retval.push_back(this->belonging_faces(*i));
        }
        return retval;
    }

    Real length_of(const edge_id_type& eid) const
    {
        return this->edge_at(eid).length;
    }
    Real tilt_angle_at(const edge_id_type& eid) const
    {
        return this->edge_at(eid).tilt;
    }
    boost::array<boost::optional<face_id_type>, 3>
    connecting_faces(const face_id_type& fid) const
    {
        boost::array<boost::optional<face_id_type>, 3> retval;
        const std::array<edge_id_type>& sides = this->edges_of(fid);
        for(std::size_t i=0; i<3; ++i)
        {
            retval[i] = this->face_of(this->opposite_of(sides[i]));
        }
        return retval;
    }
    boost::optional<face_id_type>
    connecting_faces(const face_id_type& fid, const std::size_t i) const
    {
        return this->face_of(this->opposite_of(this->edges_of(fid)[i]));
    }

    /* inherited from shape --------------------------------------------------*/
    dimension_kind dimension() const {return THREE;} // TWO?

    Real is_inside(const Real3& coord) const
    {
        throw NotImplemented("ecell4::Polygon::is_inside");
    }

    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("ecell4::Polygon::test_AABB");
    }

    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
    {
        face_id_type fid; // will be discarded
        const Real3 retval = this->draw_position(rng, fid);
        return retval;
    }

    void bounding_box(
            const Real3& edge_lengths, Real3& lower, Real3& upper) const
    {
        const Real maxr = std::numeric_limits<Real>::max();
        lower = Real3( maxr,  maxr,  maxr);
        upper = Real3(-maxr, -maxr, -maxr);

        for(typename vertex_container_type::const_iterator
            i(this->vertices_.begin()), e(this->vertices_.end()); i!=e; ++i)
        {
            lower[0] = std::min(lower[0], i->second.position[0]);
            lower[1] = std::min(lower[1], i->second.position[1]);
            lower[2] = std::min(lower[2], i->second.position[2]);

            upper[0] = std::max(upper[0], i->second.position[0]);
            upper[1] = std::max(upper[1], i->second.position[1]);
            upper[2] = std::max(upper[2], i->second.position[2]);
        }
        return;
    }

    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng,
                        face_id_type& fid) const
    {
        Real draw_triangle = rng->uniform(0.0, total_area_);

        for(std::size_t i=0; i<this->faces_.size(); ++i)
        {
            const face_data& fd = this->faces_[i];

            draw_triangle -= fd.triangle.area();
            if(draw_triangle <= 0.0)
            {
                fid = i;
                return fd.triangle.draw_position(rng);
            }
        }
        // if draw_triangle was positive throughout the loop,
        // maybe because of numerical error, put it on the last face.
        fid = faces_.back().first;
        return faces_.back().second.triangle.draw_position(rng);
    }

    // Boundary condition stuff ----------------------------------------------//

    // restrict position inside of the boundary.
    Real3 apply_boundary(const Real3& pos) const
    {
        return modulo(pos, edge_length_);
    }

    Real3 periodic_transpose(Real3 pos1, const Real3& pos2) const
    {
        const Real3 dpos = pos2 - pos1;
        const Real3 half = edge_length_ * 0.5;
        if     (half[0] <  dpos[0]) {pos1[0] += edge_length_[0];}
        else if(dpos[0] < -half[0]) {pos1[0] -= edge_length_[0];}
        if     (half[1] <  dpos[1]) {pos1[1] += edge_length_[1];}
        else if(dpos[1] < -half[1]) {pos1[1] -= edge_length_[1];}
        if     (half[2] <  dpos[2]) {pos1[2] += edge_length_[2];}
        else if(dpos[2] < -half[2]) {pos1[2] -= edge_length_[2];}
        return pos1;
    }

    std::size_t   face_size() const throw() {return    faces_.size();}
    std::size_t   edge_size() const throw() {return    edges_.size();}
    std::size_t vertex_size() const throw() {return vertices_.size();}

  private:

    vertex_data const& vertex_at(const vertex_id_type& vid) const
    {
        return this->vertices_[vid];
    }
    face_data const& face_at(const face_id_type& fid) const
    {
        return this->faces_[fid];
    }
    edge_data const& edge_at(const edge_id_type& eid) const
    {
        return this->edges_[eid];
    }

    Real distance_sq_connected_by_edges(
            const std::pair<Real3, face_id_type>& pos1,
            const std::pair<Real3, face_id_type>& pos2,
            const edge_id_type eid) const;
    Real distance_sq_connected_by_edges_impl(
            const std::pair<Real3, face_id_type>& pos1,
            const std::pair<Real3, face_id_type>& pos2,
            const edge_id_type eid,
            const vertex_data& vdata) const;

    Real distance_sq_connected_by_vertex(
            const std::pair<Real3, face_id_type>& pos1,
            const std::pair<Real3, face_id_type>& pos2,
            const vertex_id_type vid) const;

    boost::optional<edge_id_type>
    find_face_around_vtx(
            const edge_id_type start,
            const face_id_type goal,
            Real& theta) const;

    Real direction_connected_by_edges(
            const std::pair<Real3, face_id_type>& pos1,
            const std::pair<Real3, face_id_type>& pos2,
            const edge_id_type eid) const;
    Real direction_connected_by_edges_impl(
            const std::pair<Real3, face_id_type>& pos1,
            const std::pair<Real3, face_id_type>& pos2,
            const edge_id_type eid, const vertex_data& vdata) const;

    Real direction_connected_by_vertex(
            const std::pair<Real3, face_id_type>& pos1,
            const std::pair<Real3, face_id_type>& pos2,
            const vertex_id_type vid) const;

  private:

    Real  total_area_;
    Real3 edge_length_; // boundary({0,0,0}, {edge_length})

    vertex_container_type vertices_;
    face_container_type   faces_;
    edge_container_type   edges_;
    nlist_idxs_type       neighbor_list_;
    nlist_ranges_type     neighbor_list_range_;
};

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::distance_sq(
        const std::pair<Real3, face_id_type>& pos1,
        const std::pair<Real3, face_id_type>& pos2) const
{
    // case: positions are on the same face, 2D distance is same as 3D distance.
    if(pos1.second == pos2.second)
    {
        return length_sq(
            this->periodic_transpose(pos1.first, pos2.first) - pos2.first);
    }

    // case: 2 faces are connected by edges
    for(std::size_t i=0; i<3; ++i)
    {
        const edge_id_type eid = this->edges_of(pos1.second, i);
        if(this->face_of(this->opposite_of(eid)) == pos2.second)
        {
            return this->distance_sq_connected_by_edges(pos1, pos2, eid);
        }
    }

    // case: 2 positions are connected by vertex, not by edges
    const boost::array<vertex_id_type, 3>& vs = this->vertices_of(pos1.second);
    for(std::size_t i=0; i<3; ++i)
    {
        const std::vector<edge_id_type>& es = this->outgoing_edges(vs[i]);
        for(typename std::vector<edge_id_type>::const_iterator
                ei(es.begin()), ee(es.end()); ei!=ee; ++ei)
        {
            if(this->face_of(*ei) == pos2.second)
            {
                return this->distance_sq_connected_by_vertex(pos1, pos2, vs[i]);
            }
        }
    }
    // otherwise, distance cannot be calculated.
    return std::numeric_limits<Real>::infinity();
}

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::distance_sq_connected_by_edges(
        const std::pair<Real3, face_id_type>& pos1,
        const std::pair<Real3, face_id_type>& pos2,
        const edge_id_type eid) const
{
    const boost::optional<edge_id_type> opp = this->opposite_of(eid);
    if(!opp)
    {
        throw std::runtime_error("Polygon::distance_sq_connected_by_edges: "
            "invalid polygon structure");
    }

    const vertex_data& vdata1 = this->vertices_.at(this->target_of(eid));
    const vertex_data& vdata2 = this->vertices_.at(this->target_of(*opp));

    return std::min(
        this->distance_sq_connected_by_edges_impl(pos1, pos2,  eid, vdata1),
        this->distance_sq_connected_by_edges_impl(pos2, pos1, *opp, vdata2));
}

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::distance_sq_connected_by_edges_impl(
        const std::pair<Real3, face_id_type>& pos1,
        const std::pair<Real3, face_id_type>& pos2,
        const edge_id_type eid, const vertex_data& vdata) const
{
    const edge_data& edata = this->edge_at(eid);
    const Real tol_abs2 = tolerance_absolute * tolerance_absolute;
    const Real tol_rel2 = tolerance_relative * tolerance_relative;

    const Real3 v1     =
        this->periodic_transpose(vdata.position, pos1.first) - pos1.first;
    const Real3 v2     =
        this->periodic_transpose(pos2.first, vdata.position) - vdata.position;

    const Real  v1_sq  = length_sq(v1);
    const Real  v2_sq  = length_sq(v2);

    const Real  vposd2 = length_sq(vdata.position);
    if(v1_sq < tol_abs2 || v1_sq < tol_rel2 * vposd2)
    {
        return v2_len;
    }
    if(v2_sq < tol_abs2 || v2_sq < tol_rel2 * vposd2)
    {
        return v1_len;
    }
    const Real v1_len = std::sqrt(v1_sq);
    const Real v2_len = std::sqrt(v2_sq);

    Real theta =
        std::acos(dot_product(v1, edata.direction) / (v1_len * edata.length)) +
        std::acos(dot_product(v2, edata.direction) / (v2_len * edata.length));
    if(vdata.is_contiguous)
    {
        theta = std::min(theta, vdata.angle - theta);
    }
    return v1_sq + v2_sq - 2 * v1_len * v2_len * std::cos(theta);
}

template<typename T_fid, typename T_vid, typename T_eid>
boost::optional<typename Polygon<T_fid, T_vid, T_eid>::edge_id_type>
Polygon<T_fid, T_vid, T_eid>::find_face_around_vtx(
        const edge_id_type start, const face_id_type goal, Real& theta) const
{
    boost::optional<edge_id_type> eid = start;
    for(std::size_t i=0; i<10000; ++i)
    {
        eid = this->next_of(this->next_of(this->opposite_of(eid)));
        const face_id_type belonging_face_id = this->face_of(*eid);
        if(!eid || belonging_face_id == goal)
        {
            return eid;
        }

        bool found = false;
        const face_data& belonging_face_data = this->faces_.at(belonging_face_id);
        for(std::size_t i=0; i<3; ++i)
        {
            if(belonging_face_data.vertices[i] == vid)
            {
                theta += belonging_face_data.triangle.angle_at(i);
                found = true;
                break;
            }
        }
        if(!found)
        {
            throw std::logic_error("INTERNAL: ecell4::Polygon::"
                "find_face_around_vtx: invalid edge connection.");
        }

        if(*eid == start)
        {
            throw std::logic_error("INTERNAL: "
                "ecell4::Polygon::find_face_around_vtx: eid and start "
                "cannot be the same, because distance_sq already checked it.");
        }
    }
    throw std::runtime_error("INTERNAL: ecell4::Polygon::find_face_around_vtx: "
        "too many (over 10000) edges are connected to vertex.");
}

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::distance_sq_connected_by_vertex(
        const std::pair<Real3, face_id_type>& pos1,
        const std::pair<Real3, face_id_type>& pos2,
        const vertex_id_type vid) const
{
    const vertex_data& vdata = this->vertices_.at(vid);
    const Real3 v1    =
        this->periodic_transpose(vdata.position, pos1.first) - pos1.first;
    const Real3 v2    =
        this->periodic_transpose(pos2.first, vdata.position) - vdata.position;
    const Real  v1_sq = length_sq(v1);
    const Real  v2_sq = length_sq(v2);

    const Real  vposd2 = length_sq(vdata.position);
    if(v1_sq < tol_abs2 || v1_sq < tol_rel2 * vposd2)
    {
        return v2_len;
    }
    if(v2_sq < tol_abs2 || v2_sq < tol_rel2 * vposd2)
    {
        return v1_len;
    }
    const Real v1_len = std::sqrt(v1_sq);
    const Real v2_len = std::sqrt(v2_sq);

    // ------------------------------------------------------------------------
    const face_data& f1data = this->faces_.at(pos1.second);

    edge_id_type start;
    for(std::size_t i=0; i<3; ++i)
    {
        start = f1data.edges[i];
        if(this->target_of(start) == vid)
        {
            break;
        }
    }

    Real theta = std::acos(dot_product(v1, this->edge_at(start).direction) /
                           v1_len * this->edge_at(start).length);

    boost::optional<edge_id_type> eid = this->find_face_around_vtx(
            start, pos2.second, theta);

    if(eid) // found. *eid now belongs FaceID(pos2.second).
    {
        const edge_data& last_edge = this->edge_at(this->next_of(*eid));
        theta += std::acos(
            dot_product(v2, last_edge.direction) / v2_len * last_edge.length);

        if(vdata.is_contiguous)
        {
            theta = std::min(theta, vdata.angle - theta);
        }
        return v1_sq + v2_sq - 2 * v1_len * v2_len * std::cos(theta);
    }
    // else, there is a gap between p1 and p2, try reverse side...

    const face_data& f2data = this->face_at(pos2.second);
    for(std::size_t i=0; i<3; ++i)
    {
        start = f2data.edges[i];
        if(this->target_of(start) == vid)
        {
            break;
        }
    }
    theta = std::acos(dot_product(v1, this->edge_at(start).direction) /
                      v1_len * this->edge_at(start).length);

    eid = this->find_face_around_vtx(start, pos1.second, theta);
    if(eid) // found. *eid now belongs FaceID(pos1.second).
    {
        const edge_data& last_edge = this->edge_at(this->next_of(*eid));
        theta += std::acos(
            dot_product(v1, last_edge.direction) / v1_len * last_edge.length);

        assert(!vdata.is_contiguous); // gap has been found

        return v1_sq + v2_sq - 2 * v1_len * v2_len * std::cos(theta);
    }
    // couldn't find well-bahaved path between p1 and p2.
    return std::numeric_limits<Real>::infinity();
}

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::distance_sq(
        const std::pair<Real3, face_id_type>&   pos1,
        const std::pair<Real3, vertex_id_type>& pos2) const
{
    return this->distance_sq(pos2, pos1);
}

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::distance_sq(
        const std::pair<Real3, vertex_id_type>& pos1,
        const std::pair<Real3, face_id_type>&   pos2) const
{
    // it only calculates distance between vertex and position on the face
    // that connects to the vertex,
    // because the simulation algorithms does not require more than that.
    const std::vector<edge_id_type>& outs = this->outgoing_edges(pos1.second);
    for(typename std::vector<edge_id_type>::const_iterator
            i(outs.begin()), e(outs.end()); i!=e; ++i)
    {
        if(this->face_of(*i) == pos2.second)
        {
            return length_sq(this->periodic_transpose(pos1, pos2) - pos2);
        }
    }
    return std::numeric_limits<Real>::infinity();
}

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::distance_sq(
        const std::pair<Real3, vertex_id_type>& pos1,
        const std::pair<Real3, vertex_id_type>& pos2) const
{
    // it only calculates distance between vertices that are next to each other.
    // because the simulation algorithms does not require more than that.
    const std::vector<edge_id_type>& outs = this->outgoing_edges(pos1.second);
    for(typename std::vector<edge_id_type>::const_iterator
            i(outs.begin()), e(outs.end()); i!=e; ++i)
    {
        if(this->target_of(*i) == pos2.second)
        {
            const Real d = this->edge_at(*i).length;
            return d * d;
        }
    }
    return std::numeric_limits<Real>::infinity();
}

// direction ----------------------------------------------------------------

template<typename T_fid, typename T_vid, typename T_eid>
Real3 Polygon<T_fid, T_vid, T_eid>::direction(
        const std::pair<Real3, face_id_type>& pos1,
        const std::pair<Real3, face_id_type>& pos2) const
{
    // case: positions are on the same face, 2D distance is same as 3D distance.
    if(pos1.second == pos2.second)
    {
        return this->periodic_transpose(pos1.first, pos2.first) - pos2.first;
    }

    // case: 2 faces are connected by edges
    for(std::size_t i=0; i<3; ++i)
    {
        const edge_id_type eid = this->edges_of(pos1.second, i);
        if(this->face_of(this->opposite_of(eid)) == pos2.second)
        {
            return this->direction_connected_by_edges(pos1, pos2, eid);
        }
    }

    // case: 2 positions are connected by vertex, not by edges
    const boost::array<vertex_id_type, 3>& vs = this->vertices_of(pos1.second);
    for(std::size_t i=0; i<3; ++i)
    {
        const std::vector<edge_id_type>& es = this->outgoing_edges(vs[i]);
        for(typename std::vector<edge_id_type>::const_iterator
                ei(es.begin()), ee(es.end()); ei!=ee; ++ei)
        {
            if(this->face_of(*ei) == pos2.second)
            {
                return this->direction_connected_by_vertex(pos1, pos2, vs[i]);
            }
        }
    }
    throw std::out_of_range("ecell4::Polygon: direction cannot be calculated");
}

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::direction_connected_by_edges(
            const std::pair<Real3, face_id_type>& pos1,
            const std::pair<Real3, face_id_type>& pos2,
            const edge_id_type eid) const
{
    const boost::optional<edge_id_type> opp = this->opposite_of(eid);
    if(!opp)
    {
        throw std::runtime_error("Polygon::distance_sq_connected_by_edges: "
            "invalid polygon structure");
    }

    const vertex_data& vdata1 = this->vertices_.at(this->target_of(eid));
    const vertex_data& vdata2 = this->vertices_.at(this->target_of(*opp));

    const Real3 d1 =
        this->direction_connected_by_edges_impl(pos1, pos2,  eid, vdata1);
    const Real3 d2 =
        this->direction_connected_by_edges_impl(pos2, pos1, *opp, vdata2);

    return (length_sq(d1) < length_sq(d2)) ? d1 : d2;
}

template<typename T_fid, typename T_vid, typename T_eid>
Real3 Polygon<T_fid, T_vid, T_eid>::direction(
        const std::pair<Real3, face_id_type>&   pos1,
        const std::pair<Real3, vertex_id_type>& pos2) const
{
    return this->direction(pos2, pos1) * (-1.0);
}

template<typename T_fid, typename T_vid, typename T_eid>
Real3 Polygon<T_fid, T_vid, T_eid>::direction(
        const std::pair<Real3, vertex_id_type>& pos1,
        const std::pair<Real3, face_id_type>&   pos2) const
{
    // it only calculates distance between vertex and position on the face
    // that connects to the vertex,
    // because the simulation algorithms does not require more than that.
    const std::vector<edge_id_type>& outs = this->outgoing_edges(pos1.second);
    for(typename std::vector<edge_id_type>::const_iterator
            i(outs.begin()), e(outs.end()); i!=e; ++i)
    {
        if(this->face_of(*i) == pos2.second)
        {
            return this->periodic_transpose(pos2, pos1) - pos1;
        }
    }
    throw std::out_of_range("ecell4::Polygon: direction cannot be calculated");
}

template<typename T_fid, typename T_vid, typename T_eid>
Real Polygon<T_fid, T_vid, T_eid>::direction(
        const std::pair<Real3, vertex_id_type>& pos1,
        const std::pair<Real3, vertex_id_type>& pos2) const
{
    // it only calculates distance between vertices that are next to each other.
    // because the simulation algorithms does not require more than that.
    const std::vector<edge_id_type>& outs = this->outgoing_edges(pos1.second);
    for(typename std::vector<edge_id_type>::const_iterator
            i(outs.begin()), e(outs.end()); i!=e; ++i)
    {
        if(this->target_of(*i) == pos2.second)
        {
            return this->edge_at(*i).direction;
        }
    }
    throw std::out_of_range("ecell4::Polygon: direction cannot be calculated");
}

/********************** free functions to use a Polygon **********************/

namespace polygon
{

template<typename T1, typename T2, typename T3, typename s1idT, typename s2idT>
inline Real distance(const Polygon<T1, T2, T3>& p,
        const std::pair<Real3, s1idT>& p1, const std::pair<Real3, s2idT>& p2)
{
    return p.distance(p1, p2);
}
template<typename T1, typename T2, typename T3, typename s1idT, typename s2idT>
inline Real distance_sq(const Polygon<T1, T2, T3>& p,
        const std::pair<Real3, s1idT>& p1, const std::pair<Real3, s2idT>& p2)
{
    return p.distance_sq(p1, p2);
}
template<typename T1, typename T2, typename T3, typename s1idT, typename s2idT>
inline Real3 direction(const Polygon<T1, T2, T3>& p,
        const std::pair<Real3, sidT1>& p1, const std::pair<Real3, sidT2>& p2)
{
    return p.direction(p1, p2);
}

template<typename T1, typename T2, typename T3>
std::size_t travel(const Polygon<T1, T2, T3>& p,
        std::pair<Real3, typename Polygon<T1, T2, T3>::face_id_type>& pos,
        Real3& disp, std::size_t max_move_count = 10)
{
    pos = p.move_on_surface(pos, disp, max_move_count);
    return max_move_count;
}

template<typename T1, typename T2, typename T3>
inline std::pair<Real3, typename Polygon<T1, T2, T3>::face_id_type>
roll(const Polygon<T1, T2, T3>& p,
     const std::pair<Real3, typename Polygon<T1, T2, T3>::face_id_type>& pos,
     const typename Polygon<T1, T2, T3>::vertex_id_type vid,
     const Real r, const Real theta)
{
    return p.roll_around_vertex(pos, r, theta, vid);
}

} // polygon
} // ecell4
#endif// ECELL4_POLYGON
