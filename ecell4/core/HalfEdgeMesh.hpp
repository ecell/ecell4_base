#ifndef ECELL4_HALFEDGE_POLYGON
#define ECELL4_HALFEDGE_POLYGON

#include <utility>
#include <functional>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/comparators.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/geometry.hpp>
#include <ecell4/core/triangle_geometry.hpp>
#include <ecell4/core/Barycentric.hpp>

#include <boost/utility.hpp>
#include <boost/type_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/array.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/cstdint.hpp>

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <limits>

#ifndef ECELL4_STRONG_TYPEDEF // {{{
#  if __cplusplus < 201103L
#    define ECELL4_STRONG_TYPEDEF(UNDERLYING, NAME)                     \
      struct NAME {                                                     \
        UNDERLYING v_;                                                  \
        NAME(){}                                                        \
        ~NAME(){}                                                       \
        NAME(const NAME& rhs): v_(rhs.v_) {}                            \
        explicit NAME(UNDERLYING const& v): v_(v){}                     \
        NAME& operator=(const NAME& rhs){v_ = rhs.v_; return *this;}    \
        NAME& operator=(const UNDERLYING& rhs){v_ = rhs; return *this;} \
        operator UNDERLYING  () const {return v_;}                      \
        operator UNDERLYING& ()       {return v_;}                      \
        bool operator==(const NAME& rhs) const {return v_ == rhs.v_;}   \
        bool operator!=(const NAME& rhs) const {return v_ != rhs.v_;}   \
        bool operator< (const NAME& rhs) const {return v_ <  rhs.v_;}   \
        bool operator> (const NAME& rhs) const {return v_ >  rhs.v_;}   \
        bool operator<=(const NAME& rhs) const {return v_ <= rhs.v_;}   \
        bool operator>=(const NAME& rhs) const {return v_ >= rhs.v_;}   \
      };                                                                \
      /**/
#  else // c++11 enabled
#    define ECELL4_STRONG_TYPEDEF(UNDERLYING, NAME)\
       enum class NAME : UNDERLYING {};
#  endif
#endif// }}} ECELL4_STRONG_TYPEDEF

namespace ecell4
{

// Static Polygon. once made, the shape never change.
// (Face|Edge|Vertex)ID are just std::size_t.
// Any hole or duplicated vertex will be treated as invalid structure.
class HalfEdgePolygon : public Shape
{
  public:

    static const Real absolute_tolerance;
    static const Real relative_tolerance;

    // (strong|opaque) typedef
    ECELL4_STRONG_TYPEDEF(std::size_t,   face_id_type)
    ECELL4_STRONG_TYPEDEF(std::size_t,   edge_id_type)
    ECELL4_STRONG_TYPEDEF(std::size_t, vertex_id_type)

    struct vertex_data
    {
        Real  apex_angle; // total angle around this vertex
        Real3 position;   // position of this vertex
        std::vector<std::pair<edge_id_type, Real> > outgoing_edges;
        // these are garanteed to be sorted in the order that
        // (next of next of opposite of one edge on Halfedge polygon).
    };
    struct edge_data
    {
        Real           length;    // length of this edge
        Real           tilt;      // tilt from this->face to opposite face
        Real3          direction; // direction independent from boundary
        vertex_id_type target;    // destination vertex
        face_id_type   face;      // belonging face
        edge_id_type   next;      // edge on the same face, starting from target
        edge_id_type   opposite_edge;
    };
    struct face_data
    {
        Triangle triangle; // not considering Boundary, contains just the shape
        boost::array<  edge_id_type, 3> edges;    // idx consistent with triangle
        boost::array<vertex_id_type, 3> vertices; // idx consistent with triangle

        std::size_t index_of(const vertex_id_type& vid) const
        {
            if(vertices[0] == vid){return 0;}
            if(vertices[1] == vid){return 1;}
            if(vertices[2] == vid){return 2;}
            return 3;
        }
        std::size_t index_of(const edge_id_type& eid) const
        {
            if(edges[0] == eid){return 0;}
            if(edges[1] == eid){return 1;}
            if(edges[2] == eid){return 2;}
            return 3;
        }

        // neighbor list; that has pairs of {Fid, unfolded Triangle}.
        // each index corresponds to that of vertices.
        boost::array<std::vector<std::pair<face_id_type, Triangle> >, 3> neighbor_ccw;
        boost::array<std::vector<std::pair<face_id_type, Triangle> >, 3> neighbor_cw;
        // XXX: distance calculation between points on a polygon is complicated.
        // 1. there are several `local-minima` between any combination of points.
        //    so all the minima should be calculated and the shortest path
        //    should be found from them.
        // 2. each `local-minima` corresponds to a straight line along the
        //    unfolded triangles around corresponding vertices (excluding the
        //    case that is described in `3`).
        //    * if the two faces that two positions locates on are not connected
        //      by one vertices, it is not able to find the shortest path w/o
        //      employing some costly minimization technique.
        // 3. if the angle between p1-v-p2 exceeds pi, the shortest path become
        //    the path that goes through the vertex.
        //    * this means that before calculating distance, the angle should be
        //      calculated and checked whether it exceeds pi or not.
        // 4. an apex angle around a vertex can be any value (includeing inf).
        //    * by adding a pyramid instead of one face around a vertex, the
        //      angle will increase. We can put the pyramid that infinitely thin.
        //      it enables us to put any number of pyramids around one vertex,
        //      increasing the apex angle between the vertex.
    };

    typedef std::vector<vertex_data> vertex_container_type;
    typedef std::vector<  face_data>   face_container_type;
    typedef std::vector<  edge_data>   edge_container_type;

  public:

    HalfEdgePolygon(const Real3& edge_length)
        : total_area_(0.0), edge_length_(edge_length)
    {}
    HalfEdgePolygon(const Real3& edge_length, const std::vector<Triangle>& ts)
        : total_area_(0.0), edge_length_(edge_length)
    {
        this->assign(ts);
    }
    ~HalfEdgePolygon(){}

    // clear current polygon and assign different one.
    // tolerances are used to detect the same vertices in different triangles.
    void assign(const std::vector<Triangle>& ts);

//     // move `pos` to `pos + disp`.
//     std::pair<Real3, face_id_type>
//     travel(const std::pair<Real3, face_id_type>& pos, const Real3& disp) const;

    // pos1 -> pos2 <=> pos2 - pos1
//     Real3 direction (const std::pair<Real3, face_id_type>& pos1,
//                      const std::pair<Real3, face_id_type>& pos2) const;
    Real distance_sq(const std::pair<Real3, face_id_type>& pos1,
                     const std::pair<Real3, face_id_type>& pos2) const;
    Real distance   (const std::pair<Real3, face_id_type>& pos1,
                     const std::pair<Real3, face_id_type>& pos2) const
    {
        return std::sqrt(this->distance_sq(pos1, pos2));
    }

    // half-edge traverse
    // next edge: the edge belonging the same face,
    //            starting from the target of current edge
    edge_id_type next_of(const edge_id_type eid) const
    {return this->edge_at(eid).next;}

    // opposite edge: the edge that starts from the target of current edge,
    //                ends at the starting point of current edge.
    edge_id_type opposite_of(const edge_id_type eid) const
    {return this->edge_at(eid).opposite_edge;}

    // target vertex: the vertex that current edge stops at.
    vertex_id_type target_of(const edge_id_type eid) const
    {return this->edge_at(eid).target;}

    // belonging face: the face that corresponds to the current edge.
    face_id_type face_of(const edge_id_type eid) const
    {return this->edge_at(eid).face;}

    // edge ids that are corresponds to the face
    boost::array<edge_id_type, 3> const&
    edges_of(const face_id_type fid) const
    {return this->face_at(fid).edges;}

    // vertex ids that are corresponds to the face
    boost::array<vertex_id_type, 3> const&
    vertices_of(const face_id_type fid) const
    {return this->face_at(fid).vertices;}

    // edge ids that starts from the vertex
    // XXX: it returns rvalue, you cannot use it with the following form
    // for(std::vector<edge_id_type>::const_iterator
    //     i(poly.outgoing_edges().begin()), e(poly.outgoing_edges().end();
    //     i!=e; ++i)
    // {
    //     /* do some stuff ...*/;
    // }
    std::vector<edge_id_type>
    outgoing_edges(const vertex_id_type vid) const
    {
        const std::vector<std::pair<edge_id_type, Real> >&
            oedges = this->vertex_at(vid).outgoing_edges;

        std::vector<edge_id_type> retval(oedges.size());
        for(std::size_t i=0; i<oedges.size(); ++i)
        {
            retval[i] = oedges[i].first;
        }
        return retval;
    }

    // accessor ---------------------------------------------------------------

    Real apex_angle_at(const vertex_id_type& vid) const
    {
        return this->vertex_at(vid).apex_angle;
    }
    Real3 position_at(const vertex_id_type& vid) const
    {
        return this->vertex_at(vid).position;
    }
    std::vector<face_id_type>
    connecting_faces(const vertex_id_type& vid) const
    {
        std::vector<face_id_type> retval;
        const std::vector<edge_id_type>& outs = this->outgoing_edges(vid);
        for(typename std::vector<edge_id_type>::const_iterator
                i(outs.begin()), e(outs.end()); i!=e; ++i)
        {
            retval.push_back(this->face_of(*i));
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
    Real3 direction_of(const edge_id_type& eid) const
    {
        return this->edge_at(eid).direction;
    }

    Triangle const& triangle_at(const face_id_type& fid) const
    {
        return this->face_at(fid).triangle;
    }


    /* inherited from shape --------------------------------------------------*/
    dimension_kind dimension() const {return THREE;} // TWO?

    Real is_inside(const Real3&) const
    {
        throw NotImplemented("ecell4::Polygon::is_inside");
    }

    bool test_AABB(const Real3&, const Real3&) const
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

        for(vertex_container_type::const_iterator
                i(this->vertices_.begin()), e(this->vertices_.end()); i!=e; ++i)
        {
            lower[0] = std::min(lower[0], i->position[0]);
            lower[1] = std::min(lower[1], i->position[1]);
            lower[2] = std::min(lower[2], i->position[2]);

            upper[0] = std::max(upper[0], i->position[0]);
            upper[1] = std::max(upper[1], i->position[1]);
            upper[2] = std::max(upper[2], i->position[2]);
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
                fid = face_id_type(i);
                return fd.triangle.draw_position(rng);
            }
        }
        // if draw_triangle was positive throughout the loop,
        // maybe because of numerical error, put it on the last face.
        fid = face_id_type(faces_.size() - 1);
        return faces_.back().triangle.draw_position(rng);
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

    // size stuff

    Real        total_area()  const throw() {return total_area_;}
    std::size_t   face_size() const throw() {return faces_.size();}
    std::size_t   edge_size() const throw() {return edges_.size();}
    std::size_t vertex_size() const throw() {return vertices_.size();}

    std::vector<vertex_id_type> list_vertex_ids() const
    {
        std::vector<vertex_id_type> retval; retval.reserve(this->vertex_size());
        for(std::size_t i=0; i<this->vertex_size(); ++i)
        {
            retval.push_back(vertex_id_type(i));
        }
        return retval;
    }
    std::vector<face_id_type> list_face_ids() const
    {
        std::vector<face_id_type> retval; retval.reserve(this->face_size());
        for(std::size_t i=0; i<this->face_size(); ++i)
        {
            retval.push_back(face_id_type(i));
        }
        return retval;
    }
    std::vector<edge_id_type> list_edge_ids() const
    {
        std::vector<edge_id_type> retval; retval.reserve(this->edge_size());
        for(std::size_t i=0; i<this->edge_size(); ++i)
        {
            retval.push_back(edge_id_type(i));
        }
        return retval;
    }

    boost::optional<vertex_id_type>
    find_vertex(const Real3& pos) const
    {
        const Real tol_rel2 = relative_tolerance * relative_tolerance;
        const Real tol_abs2 = absolute_tolerance * absolute_tolerance;

        for(std::size_t i=0; i<vertices_.size(); ++i)
        {
            const vertex_data& vd = vertices_[i];
            const Real dist_sq = length_sq(
                this->periodic_transpose(vd.position, pos) - pos);

            if(dist_sq < tol_abs2 || dist_sq < length_sq(pos) * tol_rel2)
            {
                return vertex_id_type(i);
            }
        }
        return boost::none;
    }

    boost::optional<edge_id_type>
    find_edge(const vertex_id_type start, const vertex_id_type stop) const
    {
        const vertex_data& vd = this->vertices_.at(start);
        for(std::vector<std::pair<edge_id_type, Real> >::const_iterator
            i(vd.outgoing_edges.begin()), e(vd.outgoing_edges.end()); i!=e; ++i)
        {
            const edge_id_type eid = i->first;
            if(edges_.at(eid).target == stop)
            {
                return eid;
            }
        }
        return boost::none;
    }

  private:

    vertex_data const& vertex_at(const vertex_id_type& vid) const
    {
        return this->vertices_.at(static_cast<std::size_t>(vid));
    }
    face_data const& face_at(const face_id_type& fid) const
    {
        return this->faces_.at(static_cast<std::size_t>(fid));
    }
    edge_data const& edge_at(const edge_id_type& eid) const
    {
        return this->edges_.at(static_cast<std::size_t>(eid));
    }

    vertex_data& vertex_at(const vertex_id_type& vid)
    {
        return this->vertices_.at(static_cast<std::size_t>(vid));
    }
    face_data& face_at(const face_id_type& fid)
    {
        return this->faces_.at(static_cast<std::size_t>(fid));
    }
    edge_data& edge_at(const edge_id_type& eid)
    {
        return this->edges_.at(static_cast<std::size_t>(eid));
    }

  private:

    Real  total_area_;
    Real3 edge_length_; // boundary({0,0,0}, {edge_length})

    vertex_container_type vertices_;
    face_container_type   faces_;
    edge_container_type   edges_;
};


// /********************** free functions to use a Polygon **********************/
//
// namespace polygon
// {
//
// template<typename T1, typename T2, typename T3, typename s1idT, typename s2idT>
// inline Real distance(const Polygon<T1, T2, T3>& p,
//         const std::pair<Real3, s1idT>& p1, const std::pair<Real3, s2idT>& p2)
// {
//     return p.distance(p1, p2);
// }
// template<typename T1, typename T2, typename T3, typename s1idT, typename s2idT>
// inline Real distance_sq(const Polygon<T1, T2, T3>& p,
//         const std::pair<Real3, s1idT>& p1, const std::pair<Real3, s2idT>& p2)
// {
//     return p.distance_sq(p1, p2);
// }
// template<typename T1, typename T2, typename T3, typename s1idT, typename s2idT>
// inline Real3 direction(const Polygon<T1, T2, T3>& p,
//         const std::pair<Real3, sidT1>& p1, const std::pair<Real3, sidT2>& p2)
// {
//     return p.direction(p1, p2);
// }
//
// template<typename T1, typename T2, typename T3>
// std::size_t travel(const Polygon<T1, T2, T3>& p,
//         std::pair<Real3, typename Polygon<T1, T2, T3>::face_id_type>& pos,
//         Real3& disp, std::size_t max_move_count = 10)
// {
//     pos = p.move_on_surface(pos, disp, max_move_count);
//     return max_move_count;
// }
//
// template<typename T1, typename T2, typename T3>
// inline std::pair<Real3, typename Polygon<T1, T2, T3>::face_id_type>
// roll(const Polygon<T1, T2, T3>& p,
//      const std::pair<Real3, typename Polygon<T1, T2, T3>::face_id_type>& pos,
//      const typename Polygon<T1, T2, T3>::vertex_id_type vid,
//      const Real r, const Real theta)
// {
//     return p.roll_around_vertex(pos, r, theta, vid);
// }
//
// } // polygon
} // ecell4

#undef ECELL4_STRONG_TYPEDEF
#endif// ECELL4_POLYGON
