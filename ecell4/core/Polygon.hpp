#ifndef ECELL4_HALFEDGE_POLYGON_HPP
#define ECELL4_HALFEDGE_POLYGON_HPP

#include <utility>
#include <functional>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/comparators.hpp>
#include <ecell4/core/Integer3.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/geometry.hpp>
#include <ecell4/core/triangle_geometry.hpp>
#include <ecell4/core/Barycentric.hpp>

#ifdef WITH_HDF5
#include "PolygonHDF5Writer.hpp"
#endif

#include <boost/optional.hpp>

#include <algorithm>
#include <stdexcept>
#include <ostream>
#include <vector>
#include <limits>

namespace ecell4
{

// Static Polygon. once made, the shape never change.
// (Face|Edge|Vertex)ID are just std::size_t.
// Any hole or duplicated vertex will be treated as invalid structure.
class Polygon : public Shape
{
  public:

    static const Real absolute_tolerance;
    static const Real relative_tolerance;

    // (strong|opaque) typedef
    enum class FaceID   : std::size_t {};
    enum class EdgeID   : std::size_t {};
    enum class VertexID : std::size_t {};

    struct vertex_data
    {
        Real  apex_angle; // total angle around this vertex
        Real3 position;   // position of this vertex
        std::vector<std::pair<EdgeID, Real> > outgoing_edges;
        // these are garanteed to be sorted in the order that
        // (next of next of opposite of one edge on Halfedge polygon).
    };
    struct edge_data
    {
        Real     length;    // length of this edge
        Real     tilt;      // tilt from this->face to opposite face
        Real3    direction; // direction independent from boundary
        VertexID target;    // destination vertex
        FaceID   face;      // belonging face
        EdgeID   next;      // edge on the same face, starting from target
        EdgeID   opposite_edge;
    };
    struct face_data
    {
        Triangle triangle; // not considering Boundary, contains just the shape
        std::array<  EdgeID, 3> edges;    // idx consistent with triangle
        std::array<VertexID, 3> vertices; // idx consistent with triangle

        std::size_t index_of(const VertexID& vid) const
        {
            if(vertices[0] == vid){return 0;}
            if(vertices[1] == vid){return 1;}
            if(vertices[2] == vid){return 2;}
            return 3;
        }
        std::size_t index_of(const EdgeID& eid) const
        {
            if(edges[0] == eid){return 0;}
            if(edges[1] == eid){return 1;}
            if(edges[2] == eid){return 2;}
            return 3;
        }

        std::vector<FaceID> neighbors; // for searching objects on it
        // neighbor list; that has pairs of {Fid, unfolded Triangle}.
        // each index corresponds to that of vertices.
        std::array<std::vector<std::pair<FaceID, Triangle> >, 3> neighbor_ccw;
        std::array<std::vector<std::pair<FaceID, Triangle> >, 3> neighbor_cw;
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
        // 4. an apex angle around a vertex can be any value.
        //    * by adding a pyramid instead of one face around a vertex, the
        //      angle will increase. We can put the pyramid that infinitely thin.
        //      it enables us to put any number of pyramids around one vertex,
        //      increasing the apex angle between the vertex.
    };

    typedef std::vector<vertex_data> vertex_container_type;
    typedef std::vector<  face_data>   face_container_type;
    typedef std::vector<  edge_data>   edge_container_type;

  public:

    Polygon(const Real3&    edge_lengths =    Real3(1, 1, 1),
            const Integer3& matrix_sizes = Integer3(3, 3, 3))
        : total_area_(0.0), edge_length_(edge_lengths)
    {
        const std::size_t x_size = matrix_sizes[0];
        const std::size_t y_size = matrix_sizes[1];

        assert(x_size != 0);
        assert(y_size != 0);

        const Real dx = edge_lengths[0] / x_size;
        const Real dy = edge_lengths[1] / y_size;
        const Real z  = edge_lengths[2] / 2.0; // at the middle of Z-axis

        std::vector<Triangle> ts;
        ts.reserve(x_size * y_size * 2);

        for(std::size_t yi = 0; yi < y_size; ++yi)
        {
            for(std::size_t xi = 0; xi < x_size; ++xi)
            {
                //   4___ 3
                // y  | /|  upper left  = {1, 3, 4}
                // ^  |/_|  lower right = {1, 2, 3}
                // | 1    2
                // |
                // --> x

                const Real3 v1(dx *  xi   , dy *  yi   , z);
                const Real3 v2(dx * (xi+1), dy *  yi   , z);
                const Real3 v3(dx * (xi+1), dy * (yi+1), z);
                const Real3 v4(dx *  xi   , dy * (yi+1), z);

                ts.push_back(Triangle(v1, v3, v4));
                ts.push_back(Triangle(v1, v2, v3));
            }
        }
        this->assign(ts);
    }
    Polygon(const Real3& edge_length, const std::vector<Triangle>& ts)
        : total_area_(0.0), edge_length_(edge_length)
    {
        this->assign(ts);
    }
    ~Polygon(){}

    Polygon(const Polygon& rhs)
        : total_area_(rhs.total_area_), edge_length_(rhs.edge_length_),
          vertices_(rhs.vertices_), faces_(rhs.faces_), edges_(rhs.edges_)
    {}
    Polygon& operator=(const Polygon& rhs)
    {
        this->total_area_  = rhs.total_area_;
        this->edge_length_ = rhs.edge_length_;
        this->vertices_    = rhs.vertices_;
        this->faces_       = rhs.faces_;
        this->edges_       = rhs.edges_;
        return *this;
    }

    // clear current polygon and assign different one.
    // tolerances are used to detect the same vertices in different triangles.
    void assign(const std::vector<Triangle>& ts);

    // move `pos` to `pos + disp`.
    std::pair<Real3, FaceID>
    travel(const std::pair<Real3, FaceID>& pos, const Real3& disp) const;
    std::pair<Real3, FaceID>
    travel(const std::pair<Real3, FaceID>& pos, const Real3& disp,
           const std::size_t restraint) const;

    // pos1 -> pos2 <=> pos2 - pos1
    Real3 direction (const std::pair<Real3, FaceID>& pos1,
                     const std::pair<Real3, FaceID>& pos2) const;

    Real distance_sq(const std::pair<Real3, FaceID>& pos1,
                     const std::pair<Real3, FaceID>& pos2) const;
    Real distance   (const std::pair<Real3, FaceID>& pos1,
                     const std::pair<Real3, FaceID>& pos2) const
    {
        return std::sqrt(this->distance_sq(pos1, pos2));
    }

    Real distance_sq(const std::pair<Real3, VertexID>& pos1,
                     const std::pair<Real3, FaceID>&   pos2) const;
    Real distance   (const std::pair<Real3, VertexID>& pos1,
                     const std::pair<Real3, FaceID>&   pos2) const
    {
        return std::sqrt(this->distance_sq(pos1, pos2));
    }
    Real distance_sq(const std::pair<Real3, FaceID>& pos1,
                     const std::pair<Real3, VertexID>& pos2) const
    {
        return this->distance_sq(pos2, pos1);
    }
    Real distance   (const std::pair<Real3, FaceID>&   pos1,
                     const std::pair<Real3, VertexID>& pos2) const
    {
        return this->distance(pos2, pos1);
    }

    Real distance_sq(const std::pair<Real3, VertexID>& pos1,
                     const std::pair<Real3, VertexID>& pos2) const;
    Real distance(const std::pair<Real3, VertexID>& pos1,
                  const std::pair<Real3, VertexID>& pos2) const;

    // half-edge traverse
    // next edge: the edge belonging the same face,
    //            starting from the target of current edge
    EdgeID next_of(const EdgeID eid) const
    {return this->edge_at(eid).next;}

    // opposite edge: the edge that starts from the target of current edge,
    //                ends at the starting point of current edge.
    EdgeID opposite_of(const EdgeID eid) const
    {return this->edge_at(eid).opposite_edge;}

    // target vertex: the vertex that current edge stops at.
    VertexID target_of(const EdgeID eid) const
    {return this->edge_at(eid).target;}

    // belonging face: the face that corresponds to the current edge.
    FaceID face_of(const EdgeID eid) const
    {return this->edge_at(eid).face;}

    // edge ids that are corresponds to the face
    std::array<EdgeID, 3> const& edges_of(const FaceID fid) const
    {return this->face_at(fid).edges;}

    // vertex ids that are corresponds to the face
    std::array<VertexID, 3> const& vertices_of(const FaceID fid) const
    {return this->face_at(fid).vertices;}

    // edge ids that starts from the vertex
    // XXX: it returns rvalue, so the iterator cannot be used
    std::vector<EdgeID> outgoing_edges(const VertexID vid) const
    {
        const std::vector<std::pair<EdgeID, Real> >&
            oedges = this->vertex_at(vid).outgoing_edges;

        std::vector<EdgeID> retval(oedges.size());
        for(std::size_t i=0; i<oedges.size(); ++i)
        {
            retval[i] = oedges[i].first;
        }
        return retval;
    }

    std::vector<std::pair<EdgeID, Real> > const&
    outgoing_edge_and_angles(const VertexID vid) const
    {
        return this->vertex_at(vid).outgoing_edges;
    }

    std::vector<FaceID> neighbor_faces_of(const FaceID&   fid) const
    {
        return this->face_at(fid).neighbors;
    }
    std::vector<FaceID> neighbor_faces_of(const VertexID& vid) const
    {
        const std::vector<std::pair<EdgeID, Real> >&
            outs = this->vertex_at(vid).outgoing_edges;

        std::vector<FaceID> fids;
        fids.reserve(outs.size() * 2);
        for(std::vector<std::pair<EdgeID, Real> >::const_iterator
                i(outs.begin()), e(outs.end()); i!=e; ++i)
        {
            fids.push_back(this->face_of(i->first));
            fids.push_back(this->face_of(
                this->opposite_of(this->next_of(i->first))));
        }
        return fids;
    }

    std::vector<VertexID> neighbor_vertices_of(const FaceID&   fid) const
    {
        std::vector<VertexID> vids;
        vids.reserve(6);

        std::array<EdgeID, 3> const& eids = this->edges_of(fid);
        for(std::size_t i=0; i<3; ++i)
        {
            vids.push_back(this->target_of(eids[i]));
            vids.push_back(this->target_of(
                this->next_of(this->opposite_of(eids[i]))));
        }
        return vids;
    }
    std::vector<VertexID> neighbor_vertices_of(const VertexID& vid) const
    {
        const std::vector<std::pair<EdgeID, Real> >&
            outs = this->vertex_at(vid).outgoing_edges;
        std::vector<VertexID> vids;
        vids.reserve(outs.size() * 2);
        for(std::vector<std::pair<EdgeID, Real> >::const_iterator
                i(outs.begin()), e(outs.end()); i!=e; ++i)
        {
            vids.push_back(this->target_of(i->first));
            vids.push_back(this->target_of(
                this->next_of(this->opposite_of(this->next_of(i->first)))));
        }
        return vids;
    }

    // accessor ---------------------------------------------------------------

    Real apex_angle_at(const VertexID& vid) const
    {
        return this->vertex_at(vid).apex_angle;
    }
    Real3 position_at(const VertexID& vid) const
    {
        return this->vertex_at(vid).position;
    }
    std::vector<FaceID> connecting_faces(const VertexID& vid) const
    {
        std::vector<FaceID> retval;
        const std::vector<std::pair<EdgeID, Real> >&
            outs = this->outgoing_edge_and_angles(vid);
        for(typename std::vector<std::pair<EdgeID, Real> >::const_iterator
                i(outs.begin()), e(outs.end()); i!=e; ++i)
        {
            retval.push_back(this->face_of(i->first));
        }
        return retval;
    }

    Real length_of(const EdgeID& eid) const
    {
        return this->edge_at(eid).length;
    }
    Real tilt_angle_at(const EdgeID& eid) const
    {
        return this->edge_at(eid).tilt;
    }
    Real3 direction_of(const EdgeID& eid) const
    {
        return this->edge_at(eid).direction;
    }

    Triangle const& triangle_at(const FaceID& fid) const
    {
        return this->face_at(fid).triangle;
    }
    std::vector<FaceID> const& neighbors(const FaceID& fid) const
    {
        return this->face_at(fid).neighbors;
    }

    std::vector<std::pair<std::pair<FaceID, Triangle>, Real>>
    list_faces_within_radius(const Real3& pos, const Real& radius) const
    {
        return this->list_faces_within_radius_impl(pos, radius,
                [](const FaceID&) noexcept -> bool {
                    return false;
                });
    }
    std::vector<std::pair<std::pair<FaceID, Triangle>, Real>>
    list_faces_within_radius(const Real3& pos, const Real& radius,
                             const FaceID& ignore1) const
    {
        return this->list_faces_within_radius_impl(pos, radius,
                [=](const FaceID& fid) noexcept -> bool {
                    return fid == ignore1;
                });
    }
    std::vector<std::pair<std::pair<FaceID, Triangle>, Real>>
    list_faces_within_radius(const Real3& pos, const Real& radius,
                             const FaceID& ignore1, const FaceID& ignore2) const
    {
        return this->list_faces_within_radius_impl(pos, radius,
                [=](const FaceID& fid) noexcept -> bool {
                    return fid == ignore1 || fid == ignore2;
                });
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

    Real3 draw_position(std::shared_ptr<RandomNumberGenerator>& rng) const
    {
        FaceID fid; // will be discarded
        const Real3 retval = this->draw_position(rng, fid);
        return retval;
    }

    //XXX This function randomly choose a face and does not consider the
    //    initial value of `fid`. An ID of randomly-chosen face will be
    //    written in `fid` when it returns.
    Real3 draw_position(std::shared_ptr<RandomNumberGenerator>& rng,
                        FaceID& fid) const
    {
        Real draw_triangle = rng->uniform(0.0, total_area_);

        for(std::size_t i=0; i<this->faces_.size(); ++i)
        {
            const face_data& fd = this->faces_[i];

            draw_triangle -= fd.triangle.area();
            if(draw_triangle <= 0.0)
            {
                fid = FaceID(i);
                return fd.triangle.draw_position(rng);
            }
        }
        // if draw_triangle was positive throughout the loop,
        // maybe because of numerical error, put it on the last face.
        fid = FaceID(faces_.size() - 1);
        return faces_.back().triangle.draw_position(rng);
    }

    // XXX This function considers `fid`.
    Real3 draw_position_on_face(std::shared_ptr<RandomNumberGenerator>& rng,
                                const FaceID& fid) const
    {
        return this->triangle_at(fid).draw_position(rng);
    }

    // to save/load the shape ------------------------------------------------//

    void reset(const Real3& edge_lengths)
    {
        this->edge_length_ = edge_lengths;
        return;
    }

    std::vector<Triangle> triangles() const
    {
        std::vector<Triangle> retval;
        retval.reserve(this->faces_.size());
        for(const auto& f : this->faces_)
        {
            retval.push_back(f.triangle);
        }
        return retval;
    }

#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const
    {
        save_triangles_polygon(*this, root);
    }

    void load_hdf5(const H5::Group& root)
    {
        load_triangles_polygon(root, this);
    }
#endif

    // Boundary condition stuff ----------------------------------------------//

    Real3 const& edge_lengths() const throw() {return this->edge_length_;}

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

    bool is_inside_of_boundary(const Real3& pos) const
    {
        const Real3& edge = this->edge_length_;
        if     (edge[0] < pos[0]) {return false;}
        else if( pos[0] <    0.0) {return false;}
        if     (edge[1] < pos[1]) {return false;}
        else if( pos[1] <    0.0) {return false;}
        if     (edge[2] < pos[2]) {return false;}
        else if( pos[2] <    0.0) {return false;}
        return true;
    }

    // size stuff

    Real        total_area()  const throw() {return total_area_;}
    std::size_t   face_size() const throw() {return faces_.size();}
    std::size_t   edge_size() const throw() {return edges_.size();}
    std::size_t vertex_size() const throw() {return vertices_.size();}

    std::vector<VertexID> list_vertex_ids() const
    {
        std::vector<VertexID> retval; retval.reserve(this->vertex_size());
        for(std::size_t i=0; i<this->vertex_size(); ++i)
        {
            retval.push_back(VertexID(i));
        }
        return retval;
    }
    std::vector<FaceID> list_face_ids() const
    {
        std::vector<FaceID> retval; retval.reserve(this->face_size());
        for(std::size_t i=0; i<this->face_size(); ++i)
        {
            retval.push_back(FaceID(i));
        }
        return retval;
    }
    std::vector<EdgeID> list_edge_ids() const
    {
        std::vector<EdgeID> retval; retval.reserve(this->edge_size());
        for(std::size_t i=0; i<this->edge_size(); ++i)
        {
            retval.push_back(EdgeID(i));
        }
        return retval;
    }

    boost::optional<VertexID>
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
                return VertexID(i);
            }
        }
        return boost::none;
    }

    boost::optional<EdgeID>
    find_edge(const VertexID start, const VertexID stop) const
    {
        const vertex_data& vd = this->vertex_at(start);
        for(std::vector<std::pair<EdgeID, Real> >::const_iterator
            i(vd.outgoing_edges.begin()), e(vd.outgoing_edges.end()); i!=e; ++i)
        {
            const EdgeID eid = i->first;
            if(this->target_of(eid) == stop)
            {
                return eid;
            }
        }
        return boost::none;
    }

    boost::optional<FaceID>
    find_face(const VertexID v1, const VertexID v2,
              const VertexID v3) const
    {
        const std::vector<EdgeID>& outs = this->outgoing_edges(v1);
        for(typename std::vector<EdgeID>::const_iterator
                i(outs.begin()), e(outs.end()); i!=e; ++i)
        {
            if(target_of(*i) == v2 && target_of(next_of(*i)) == v3)
            {
                return face_of(*i);
            }
            if(target_of(*i) == v3 && target_of(next_of(*i)) == v2)
            {
                return face_of(*i);
            }
        }
        return boost::none;
    }


  private:

    vertex_data const& vertex_at(const VertexID& vid) const
    {
        return this->vertices_.at(static_cast<std::size_t>(vid));
    }
    face_data const& face_at(const FaceID& fid) const
    {
        return this->faces_.at(static_cast<std::size_t>(fid));
    }
    edge_data const& edge_at(const EdgeID& eid) const
    {
        return this->edges_.at(static_cast<std::size_t>(eid));
    }

    vertex_data& vertex_at(const VertexID& vid)
    {
        return this->vertices_.at(static_cast<std::size_t>(vid));
    }
    face_data& face_at(const FaceID& fid)
    {
        return this->faces_.at(static_cast<std::size_t>(fid));
    }
    edge_data& edge_at(const EdgeID& eid)
    {
        return this->edges_.at(static_cast<std::size_t>(eid));
    }

    template<typename Filter>
    std::vector<std::pair<std::pair<FaceID, Triangle>, Real>>
    list_faces_within_radius_impl(
            const Real3& pos, const Real& radius, Filter filter) const
    {
        std::vector<std::pair<std::pair<FaceID, Triangle>, Real>> retval;

        const Real radius_sq = radius * radius;
        for(std::size_t i=0; i<faces_.size(); ++i)
        {
            const FaceID fid = static_cast<FaceID>(i);
            if(filter(fid))
            {
                continue;
            }

            const Triangle& tri = faces_[i].triangle;
            const Real dist_sq  = distance_sq_point_Triangle(pos, tri);
            if(dist_sq <= radius_sq)
            {
                retval.emplace_back(std::make_pair(fid, tri), std::sqrt(dist_sq));
            }
        }
        std::sort(retval.begin(), retval.end(),
            utils::pair_second_element_comparator<std::pair<FaceID, Triangle>, Real>{});

        return retval;
    }

  private:

    Real  total_area_;
    Real3 edge_length_; // boundary({0,0,0}, {edge_length})

    vertex_container_type vertices_;
    face_container_type   faces_;
    edge_container_type   edges_;
};

template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Polygon::FaceID& fid)
{
    os << "FID(" << static_cast<std::size_t>(fid) << ")";
    return os;
}

template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Polygon::EdgeID& eid)
{
    os << "EID(" << static_cast<std::size_t>(eid) << ")";
    return os;
}

template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Polygon::VertexID& vid)
{
    os << "VID(" << static_cast<std::size_t>(vid) << ")";
    return os;
}

/********************** free functions to use a Polygon **********************/

namespace polygon
{

template<typename ID1, typename ID2>
inline Real distance(const Polygon& p,
        const std::pair<Real3, ID1>& p1, const std::pair<Real3, ID2>& p2)
{
    return p.distance(p1, p2);
//     const Real d12 = p.distance(p1, p2);
//     const Real d21 = p.distance(p2, p1);
//     if(std::abs(d12 / d21) - 1.0 > 1e-12)
//     {
//         std::cerr << "distance btw 1-2 = " << std::setw(20) << d12 << ", 2-1 = " << std::setw(20) << d21
//                   << ", the difference = " << std::setw(20) << std::abs(d12 / d21 - 1.0);
//         std::cerr << ",  p1 = (" << p1.first << ", " << p1.second
//                   << "), p2 = (" << p2.first << ", " << p2.second
//                   << ')' << std::endl;
//     }
//     return std::min(d12, d21);
}

template<typename ID1, typename ID2>
inline Real distance_sq(const Polygon& p,
        const std::pair<Real3, ID1>& p1, const std::pair<Real3, ID2>& p2)
{
    return p.distance_sq(p1, p2);
//     const Real d12 = p.distance_sq(p1, p2);
//     const Real d21 = p.distance_sq(p2, p1);
//     if(std::abs(d12 / d21) - 1.0 > 1e-12)
//     {
//         std::cerr << "distance_sq btw 1-2 = " << std::setw(20) << d12 << ", 2-1 = " << std::setw(20) << d21
//                   << ", the difference = " << std::setw(20) << std::abs(d12 / d21 - 1.0);
//         std::cerr << ",  p1 = (" << p1.first << ", " << p1.second
//                   << "), p2 = (" << p2.first << ", " << p2.second
//                   << ')' << std::endl;
//     }
//     return std::min(d12, d21);
}

template<typename ID1, typename ID2>
inline Real3 direction(const Polygon& p,
        const std::pair<Real3, ID1>& p1, const std::pair<Real3, ID2>& p2)
{
    return p.direction(p1, p2);
}

inline std::pair<Real3, Polygon::FaceID> travel(const Polygon& p,
        const std::pair<Real3, Polygon::FaceID>& pos, const Real3& disp)
{
    return p.travel(pos, disp);
}

inline std::pair<Real3, Polygon::FaceID> travel(const Polygon& p,
        const std::pair<Real3, Polygon::FaceID>& pos, const Real3& disp,
        const std::size_t edge_restraint)
{
    return p.travel(pos, disp, edge_restraint);
}

// for sGFRD
inline std::pair<Real3, Polygon::FaceID>
roll(const Polygon& poly,
     const std::pair<Real3, Polygon::FaceID>& pos,
     const Polygon::VertexID vid, const Real r, const Real theta_)
{
    const Real3 vpos = poly.periodic_transpose(poly.position_at(vid), pos.first);
    if(theta_ == 0.0)
    {
        const Real3  disp = pos.first - vpos;
        const Real3 nwpos = poly.apply_boundary(vpos + disp * (r / length(disp)));
        assert(::ecell4::is_inside(
               ::ecell4::to_barycentric(nwpos, poly.triangle_at(pos.second))));
        return std::make_pair(nwpos, pos.second);
    }

    const Real apex_angle = poly.apex_angle_at(vid);
    assert(std::abs(theta_) < apex_angle);

    typedef Polygon::FaceID FaceID;
    typedef Polygon::EdgeID EdgeID;

    std::vector<std::pair<EdgeID, Real> > const& outedges =
        poly.outgoing_edge_and_angles(vid);

    std::vector<std::pair<EdgeID, Real> >::const_iterator current_edge;
    Real current_angle = 0.0;
    for(std::vector<std::pair<EdgeID, Real> >::const_iterator
            i(outedges.begin()), e(outedges.end()); i!=e; ++i)
    {
        if(poly.face_of(i->first) == pos.second)
        {
            current_edge = i;
            break;
        }
        current_angle += i->second;
    }
    current_angle += calc_angle(
            poly.direction_of(current_edge->first), pos.first - vpos);

    const Real theta = ::ecell4::modulo(current_angle + theta_, apex_angle);
//     std::cerr << "initial position = " << pos.first << std::endl;
//     std::cerr << "initial first    = " << pos.second << ", " << poly.triangle_at(pos.second) << std::endl;
//     std::cerr << "current_angle    = " << current_angle << ", theta_ = " << theta_
//               << ", apex_angle = "  <<    apex_angle << std::endl;
//     std::cerr << "theta = " << theta << std::endl;

    current_angle = 0.0;
    for(std::vector<std::pair<EdgeID, Real> >::const_iterator
            i(outedges.begin()), e(outedges.end()); i!=e; ++i)
    {
//         {
//             std::cerr << "face " << poly.face_of(i->first) << " = "
//                       << poly.triangle_at(poly.face_of(i->first)) << std::endl;
//             std::cerr << "angle around vtx " << vid << " at " << vpos << " = "
//                       << i->second << std::endl;
//         }

        const Real next_angle = current_angle + i->second;
        if(theta < next_angle)
        {
            const FaceID next_fid = poly.face_of(i->first);
            const Real3& normal   = poly.triangle_at(next_fid).normal();
            const Real diff_angle = theta - current_angle;
            const Real3 direction = rotate(
                    diff_angle, normal, poly.direction_of(i->first));
            const Real3 next_pos  = poly.apply_boundary(
                    vpos + direction * (r / length(direction)));

            if(!::ecell4::is_inside(
                        to_barycentric(next_pos, poly.triangle_at(next_fid))))
            {
                std::cerr << "position = " << next_pos << ", face " << next_fid
                          << " = " <<  poly.triangle_at(next_fid) << std::endl;
                assert(false);
            }
            return std::make_pair(next_pos, next_fid);
        }
        current_angle = next_angle;
    }
    throw std::logic_error("ecell4::polygon::roll never reach here");
}

} // polygon
} // ecell4

namespace std {

template<>
struct hash<ecell4::Polygon::FaceID>
{
    typedef std::size_t result_type;
    typedef ecell4::Polygon::FaceID argument_type;
    result_type operator()(const argument_type& val) const
    {
        return static_cast<std::size_t>(val);
    }
};
template<>
struct hash<ecell4::Polygon::EdgeID>
{
    typedef std::size_t result_type;
    typedef ecell4::Polygon::EdgeID argument_type;
    result_type operator()(const argument_type& val) const
    {
        return static_cast<std::size_t>(val);
    }
};
template<>
struct hash<ecell4::Polygon::VertexID>
{
    typedef std::size_t result_type;
    typedef ecell4::Polygon::VertexID argument_type;
    result_type operator()(const argument_type& val) const
    {
        return static_cast<std::size_t>(val);
    }
};
} // std

#endif// ECELL4_POLYGON
