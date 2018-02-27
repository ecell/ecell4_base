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
#include <ecell4/core/ObjectIDContainer.hpp>

#include <boost/utility.hpp>
#include <boost/type_traits.hpp>
#include <boost/array.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/cstdint.hpp>

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <limits>

namespace ecell4
{

struct VertexID : public Identifier<VertexID, boost::uint64_t, boost::uint8_t>
{
    typedef Identifier<VertexID, boost::uint64_t, boost::uint8_t> base_type;
    VertexID(const value_type& v = value_type(0, 0)) : base_type(v){}
};
template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>&
operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const VertexID& v)
{
    strm << "VtxID(" << static_cast<boost::uint32_t>(v.lot())
         << ':' << v.serial() << ')';
    return strm;
}

struct FaceID : public Identifier<FaceID, boost::uint64_t, boost::uint8_t>
{
    typedef Identifier<FaceID, boost::uint64_t, boost::uint8_t> base_type;
    FaceID(const value_type& v = value_type(0, 0)) : base_type(v){}
};
template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>&
operator<<(std::basic_ostream<Tstrm_, Ttraits_>& strm, const FaceID& v)
{
    strm << "FaceID(" << static_cast<boost::uint32_t>(v.lot())
         << ':' << v.serial() << ')';
    return strm;
}
} // ecell4

ECELL4_DEFINE_HASH_BEGIN()
template<>
struct hash<ecell4::FaceID>
{
    std::size_t operator()(const ecell4::FaceID& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};

template<>
struct hash<ecell4::VertexID>
{
    std::size_t operator()(const ecell4::VertexID& val) const
    {
        return static_cast<std::size_t>(val().first ^ val().second);
    }
};
ECELL4_DEFINE_HASH_END()

namespace ecell4
{

// Static Polygon. once made, the shape never change.
template<typename T_face_id,
         typename T_fid_gen   = SerialIDGenerator<T_face_id>,
         typename T_vertex_id = VertexID,
         typename T_vid_gen   = SerialIDGenerator<T_vertex_id> >
class Polygon : public Shape
{
  public:
    typedef T_vertex_id vertex_id_type;
    typedef T_face_id   face_id_type;
    typedef T_vid_gen   vertex_id_generator_type;
    typedef T_fid_gen   face_id_generator_type;

    struct vertex_data
    {
        Real  angle; // total angle around this vertices
        Real3 position;
        // faces connecting to this vertex.
        // if there are blank (non-connecting faces), boost::none will be
        // inserted between them.
        std::vector<boost::optional<std::pair<face_id_type, std::size_t> >
            > faces;
    };
    struct face_data
    {
        Triangle     triangle;
        boost::array<vertex_id_type, 3> vertices;
        boost::array<boost::optional<std::pair<face_id_type, Real> >, 3> faces;
        // {id, tilt}. if there are no face connected, it become boost::none.
    };

    typedef ObjectIDContainer<vertex_id_type, vertex_data> vertex_container_type;
    typedef ObjectIDContainer<face_id_type,   face_data>   face_container_type;

  public:

    // tolerances are used to detect the same vertices in different triangles.
    Polygon(const Real3& edge_length, const std::vector<Triangle>& ts,
            const Real tol_abs = 1e-12, const Real tol_rel = 1e-8);
    ~Polygon(){}

    std::pair<Real3, face_id_type>
    travel(const std::pair<Real3, face_id_type>& pos, const Real3& disp,
           std::size_t& max_edges_traversed);

    std::pair<Real3, face_id_type>
    roll_around(const std::pair<Real3, face_id_type>& pos,
                const Real r, const Real theta, const vertex_id_type& vid);

    Real3 direction(const std::pair<Real3,   face_id_type>& pos1,
                    const std::pair<Real3,   face_id_type>& pos2);
    Real3 direction(const std::pair<Real3, vertex_id_type>& pos1,
                    const std::pair<Real3, vertex_id_type>& pos2);
    Real3 direction(const std::pair<Real3,   face_id_type>& pos1,
                    const std::pair<Real3, vertex_id_type>& pos2);
    Real3 direction(const std::pair<Real3, vertex_id_type>& pos1,
                    const std::pair<Real3,   face_id_type>& pos2);

    Real distance_sq(const std::pair<Real3,   face_id_type>& pos1,
                     const std::pair<Real3,   face_id_type>& pos2);
    Real distance_sq(const std::pair<Real3, vertex_id_type>& pos1,
                     const std::pair<Real3, vertex_id_type>& pos2);
    Real distance_sq(const std::pair<Real3,   face_id_type>& pos1,
                     const std::pair<Real3, vertex_id_type>& pos2);
    Real distance_sq(const std::pair<Real3, vertex_id_type>& pos1,
                     const std::pair<Real3,   face_id_type>& pos2);

    Real distance(const std::pair<Real3,   face_id_type>& pos1,
                  const std::pair<Real3,   face_id_type>& pos2)
    {return std::sqrt(distance_sq(pos1, pos2));}
    Real distance(const std::pair<Real3, vertex_id_type>& pos1,
                  const std::pair<Real3, vertex_id_type>& pos2)
    {return std::sqrt(distance_sq(pos1, pos2));}
    Real distance(const std::pair<Real3,   face_id_type>& pos1,
                  const std::pair<Real3, vertex_id_type>& pos2)
    {return std::sqrt(distance_sq(pos1, pos2));}
    Real distance(const std::pair<Real3, vertex_id_type>& pos1,
                  const std::pair<Real3,   face_id_type>& pos2)
    {return std::sqrt(distance_sq(pos1, pos2));}

    /* accessor --------------------------------------------------------------*/
    // pass std::nothrow explicitly to call noexcept version.

    Real apex_angle_at(const vertex_id_type& vid) const
    {
        return this->vertex_at(vid).angle;
    }
    Real apex_angle_at(const vertex_id_type& vid, std::nothrow_t) const throw()
    {
        return this->vertices_[vid].angle;
    }
    Real position_at(const vertex_id_type& vid) const
    {
        return this->vertex_at(vid).position;
    }
    Real position_at(const vertex_id_type& vid, std::nothrow_t) const throw()
    {
        return this->vertices_[vid].position;
    }
    std::vector<face_id_type> const&
    connecting_faces(const vertex_id_type& vid) const
    {
        return this->vertex_at(vid).faces;
    }
    std::vector<face_id_type> const&
    connecting_faces(const vertex_id_type& vid, std::nothrow_t) const throw()
    {
        return this->vertices_[vid].faces;
    }

    Triangle const& triangle_at(const face_id_type& fid) const
    {
        return this->faces_.at(fid).triangle;
    }
    Triangle const&
    triangle_at(const face_id_type& fid, std::nothrow_t) const throw()
    {
        return this->faces_[fid].triangle;
    }
    boost::array<vertex_id_type, 3> const&
    connecting_vertices(const face_id_type& fid) const
    {
        return this->faces_.at(fid).vertices;
    }
    boost::array<vertex_id_type, 3> const&
    connecting_vertices(const face_id_type& fid, std::nothrow_t) const throw()
    {
        return this->faces_[fid].vertices;
    }
    boost::array<boost::optional<std::pair<face_id_type, Real> >, 3> const&
    connecting_faces(const face_id_type& fid) const
    {
        return this->faces_.at(fid).faces;
    }
    boost::array<boost::optional<std::pair<face_id_type, Real> >, 3> const&
    connecting_faces(const face_id_type& fid, std::nothrow_t) const throw()
    {
        return this->faces_[fid].faces;
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
                i(vertices_.begin()), e(vertices_.end()); i != e; ++i)
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
        for(typename face_container_type::const_iterator
                i(faces_.begin()), e(faces_.end()); i != e; ++i)
        {
            draw_triangle -= i->second.triangle.area();
            if(draw_triangle <= 0.0)
            {
                fid = i->first;
                return i->second.triangle.draw_position(rng);
            }
        }
        // if draw_triangle was positive throughout the loop,
        // maybe because of numerical error, put it on the last face.
        fid = (faces_.end() - 1)->first;
        return (faces_.end() - 1)->second.triangle.draw_position(rng);
    }

    /* accessor: mainly for debug or internal use ----------------------------*/

    vertex_data const& vertex_at(const vertex_id_type& vid) const
    {
        return this->vertices_.at(vid);
    }
    face_data const& face_at(const face_id_type& fid) const
    {
        return this->faces_.at(fid);
    }

    vertex_data const&
    vertex_at(const vertex_id_type& vid, std::nothrow_t) const throw()
    {
        return this->vertices_[this->vertex_id_map_[vid]];
    }
    face_data const&
    face_at(const face_id_type& vid, std::nothrow_t) const throw()
    {
        return this->faces_[this->face_id_map_[vid]];
    }

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

    std::vector<face_id_type> list_face_id() const
    {
        std::vector<face_id_type> list; list.reserve(faces_.size());
        for(typename face_container_type::const_iterator
                i(faces_.begin()), e(faces_.end()); i!=e; ++i)
        {
            list.push_back(i->first);
        }
        return list;
    }

    std::vector<vertex_id_type> list_vertex_id() const
    {
        std::vector<vertex_id_type> list; list.reserve(vertices_.size());
        for(typename vertex_container_type::const_iterator
                i(vertices_.begin()), e(vertices_.end()); i!=e; ++i)
        {
            list.push_back(i->first);
        }
        return list;
    }

    Real area() const throw() {return this->total_area_;}

  private:

    Real  total_area_;
    Real3 edge_length_; // boundary({0,0,0}, {edge_length})

    vertex_id_generator_type vidgen_;
    vertex_container_type    vertices_;
    face_id_generator_type   fidgen_;
    face_container_type      faces_;
};

template<typename T_fid, typename T_fidgen, typename T_vid, typename T_vidgen>
Polygon<T_fid, T_fidgen, T_vid, T_vidgen>::Polygon(const Real3& el,
    const std::vector<Triangle>& ts, const Real tol_abs, const Real tol_rel)
    : edge_length_(el)
{
    typedef std::pair<face_id_type, std::size_t> fid_index_pair;
    std::vector<std::vector<fid_index_pair> >    identical_vtx_list;

    total_area_ = 0.0;
    for(typename std::vector<Triangle>::const_iterator
            ti(ts.begin()), te(ts.end()); ti != te; ++ti)
    {
        const Triangle& t = *ti;
        const face_id_type fid = fidgen_();

        face_data fd; fd.triangle = t;
        faces_.update(fid, fd);

        total_area_ += t.area();

        for(std::size_t idx=0; idx<3; ++idx)
        {
            const Real3& v = t.vertices()[idx];
            bool found = false;
            for(typename std::vector<std::vector<fid_index_pair> >::iterator
                    vi(identical_vtx_list.begin()), ve(identical_vtx_list.end());
                    vi != ve; ++vi)
            {
                const Real3 v_ = this->periodic_transpose(faces_.at(
                    vi->front().first).second.triangle.vertex_at(
                    vi->front().second), v);
                const Real l = length(v - v_);

                if(l < tol_abs || l / length(v) < tol_rel)
                {
                    vi->push_back(std::make_pair(fid, idx));
                    found = true;
                    break;
                }
            }
            if(!found)
            {
                identical_vtx_list.push_back(std::vector<fid_index_pair>(1,
                            std::make_pair(fid, idx)));
            }
        }
    }

    for(typename std::vector<std::vector<fid_index_pair> >::iterator
            vli(identical_vtx_list.begin()), vle(identical_vtx_list.end());
            vli != vle; ++vli)
    {
        const vertex_id_type vid = vidgen_();
        vertex_data vd;

        vd.position = Real3(0.0, 0.0, 0.0);
        for(typename std::vector<fid_index_pair>::const_iterator
                vi(vli->begin()), ve(vli->end()); vi != ve; ++vi)
        {
            vd.faces.push_back(*vi);
            vd.position += faces_[vi->first].second.triangle.vertices()[vi->second];
        }
        vd.position /= static_cast<Real>(vli->size());

        // update position of each triangle...
        for(typename std::vector<fid_index_pair>::const_iterator
                vi(vli->begin()), ve(vli->end()); vi != ve; ++vi)
        {
            boost::array<Real3, 3> vs = faces_[vi->first].second.triangle.vertices();
            vs[vi->second] = vd.position;

            faces_[vi->first].second.triangle = Triangle(vs);
            faces_[vi->first].second.vertices[vi->second] = vid;
        }

        vertices_.update(vid, vd);
    }

    // detect face-face connection
    for(typename face_container_type::iterator
            fi1(faces_.begin()), fe1(faces_.end()); fi1 != fe1; ++fi1)
    {
        const face_id_type self = fi1->first;
        for(std::size_t vidx1=0; vidx1 < 3; ++vidx1)
        {
            // if already found
            if(fi1->second.faces[vidx1]) {continue;}

            const vertex_id_type vid1 = fi1->second.vertices[vidx1];
            const vertex_id_type vid2 = fi1->second.vertices[
                (vidx1==2) ? 0 : vidx1+1];

            // so far, none of the elements of vertex_data.faces are boost::none.
            const std::vector<boost::optional<
                std::pair<face_id_type, std::size_t>
                > >& faces_v1 = vertices_[vid1].second.faces;
            const std::vector<boost::optional<
                std::pair<face_id_type, std::size_t>
                > >& faces_v2 = vertices_[vid2].second.faces;

            bool found = false;
            for(std::size_t i=0; i<faces_v1.size(); ++i)
            {
                if(faces_v1[i]->first == self) {continue;}
                const face_id_type cand = faces_v1[i]->first;
                for(std::size_t j=0; j<faces_v2.size(); ++j)
                {
                    if(faces_v2[j]->first == self) {continue;}
                    if(cand == faces_v2[j]->first)
                    {
                        const std::size_t vidx2 = faces_v2[j]->second;
                        const Triangle& partner = faces_.at(cand).second.triangle;
                        found = true;
                        //  ^n
                        // _I_....       n
                        //    \) theta   ^  / -theta
                        //     \        _I_/...
                        //
                        // to consider local convex / concave, use copysign.
                        Real agl =
                            angle(fi1->second.triangle.normal(), partner.normal());
                        agl = boost::math::copysign(agl, -1 * dot_product(
                            fi1->second.triangle.normal(),
                            partner.edges()[(vidx2==2) ? 0 : vidx2+1]));

                        fi1->second.faces[vidx1]         = std::make_pair(cand, agl);
                        faces_[cand].second.faces[vidx2] = std::make_pair(self, agl);
                        break;
                    }
                }
                if(found) {break;}
            }
            if(!found) {fi1->second.faces[vidx1] = boost::none;}
        }
    }

    // sort vertex_data::faces
    for(typename vertex_container_type::iterator
            vi(vertices_.begin()), ve(vertices_.end()); vi!=ve; ++vi)
    {
        std::vector<boost::optional<std::pair<face_id_type, std::size_t> >
            > original = vi->second.faces;
        std::vector<boost::optional<std::pair<face_id_type, std::size_t> >
            > sorted; sorted.reserve(vi->second.faces.size());

        while(true)
        {
            const boost::optional<std::pair<face_id_type, std::size_t> > start =
                original.back(); original.pop_back();
            sorted.push_back(start);

            boost::optional<std::pair<face_id_type, std::size_t> > current;
            {
                const boost::optional<std::pair<face_id_type, Real> >& next =
                    faces_[start->first].second
                    .faces[start->second == 0 ? 2 : start->second - 1];
                if(next)
                {
                    const boost::array<vertex_id_type, 3>& vs =
                        faces_[next->first].second.vertices;

                    for(std::size_t i=0; i<3; ++i)
                    {
                        if(vs[i] == vi->first)
                        {
                            current = std::make_pair(next->first, i);
                            break;
                        }
                    }
                    assert(current);
                }
                else
                {
                    current = boost::none;
                }
            }

            while(current && *current != *start)
            {
                typename std::vector<boost::optional<
                    std::pair<face_id_type, std::size_t> > >::iterator found =
                    std::find(original.begin(), original.end(), *current);

                assert(found != original.end());

                original.erase(found);
                sorted.push_back(*current);

                {
                    const boost::optional<std::pair<face_id_type, Real> >& next =
                        faces_[current->first].second
                        .faces[current->second == 0 ? 2 : current->second - 1];

                    current = boost::none;
                    if(next)
                    {
                        const boost::array<vertex_id_type, 3>& vs =
                            faces_[next->first].second.vertices;

                        for(std::size_t i=0; i<3; ++i)
                        {
                            if(vs[i] == vi->first)
                            {
                                current = std::make_pair(next->first, i);
                                break;
                            }
                        }
                        assert(current);
                    }
                }
            }

            if(original.empty())
            {
                break;
            }
            else
            {
                sorted.push_back(boost::none);
            }
        }
        vi->second.faces = sorted;
    }

    // calc total angle
    for(typename vertex_container_type::iterator
            vci(vertices_.begin()), vce(vertices_.end()); vci != vce; ++vci)
    {
        vci->second.angle = 0.0;
        for(typename std::vector<boost::optional<
                std::pair<face_id_type, std::size_t> > >::const_iterator
                vi(vci->second.faces.begin()), ve(vci->second.faces.end());
                vi != ve; ++vi)
        {
            boost::optional<std::pair<face_id_type, std::size_t> >const&
                f_idx = *vi;
            if(f_idx)
            {
                vci->second.angle +=
                    faces_[f_idx->first].second.triangle.angle_at(f_idx->second);
            }
        }
    }
}

/********************** free functions to use a Polygon **********************/

namespace polygon
{

template<typename T, typename s1idT, typename s2idT>
inline Real distance(const Polygon<T>& p,
        const std::pair<Real3, s1idT>& p1, const std::pair<Real3, s2idT>& p2)
{
    return p.distance(p1, p2);
}
template<typename T, typename s1idT, typename s2idT>
inline Real distance_sq(const Polygon<T>& p,
        const std::pair<Real3, s1idT>& p1, const std::pair<Real3, s2idT>& p2)
{
    return p.distance_sq(p1, p2);
}

template<typename T, typename sidT1, typename sidT2>
inline Real3 direction(const Polygon<T>& p,
        const std::pair<Real3, sidT1>& p1, const std::pair<Real3, sidT2>& p2)
{
    return p.direction(p1, p2);
}

template<typename T>
std::size_t
travel(const Polygon<T>& p, std::pair<Real3, typename T::face_id_type>& pos,
       Real3& disp, std::size_t max_move_count = 10)
{
    pos = p.travel(pos, disp, max_move_count);
    return max_move_count;
}

template<typename T>
inline std::pair<Real3, typename T::face_id_type>
roll(const Polygon<T>& p, const std::pair<Real3, typename T::face_id_type>& pos,
     typename T::vertex_id_type vid, const Real r, const Real theta)
{
    return p.roll_around(pos, r, theta, vid);
}

} // polygon
} // ecell4
#endif// ECELL4_POLYGON
