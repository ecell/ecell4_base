#ifndef ECELL4_EGFRD_POLYGON_HPP
#define ECELL4_EGFRD_POLYGON_HPP

#include <ecell4/core/Shape.hpp>
#include <ecell4/core/AABB.hpp>
#include <ecell4/core/Real3.hpp>
#include "exceptions.hpp"
#include "FaceTriangle.hpp"
#include <limits>
#include <algorithm>

namespace ecell4
{
namespace egfrd
{
template<typename coordT>
struct Polygon : public ecell4::Shape
{
    typedef coordT coordinate_type;
    typedef FaceTriangle<coordinate_type> face_type;
    typedef std::size_t face_id_type;

    static face_id_type make_nonsence_id()
    {
        return std::numeric_limits<std::size_t>::max();
    }

    Polygon(){}
    ~Polygon(){}

    // nearest
    // (intruders, (idx, (distance, max_radius)))
    std::pair<std::vector<std::size_t>, std::pair<std::size_t, std::pair<Real, Real> > >
    get_faces_within_radius(const coordinate_type& pos, const Real range) const;

    void emplace(const std::array<coordinate_type, 3>& vertices)
    {
        this->faces.push_back(face_type(vertices));
    }

// data member
    std::vector<face_type> faces;

// for shapes (not implemented yet)
    dimension_kind dimension() const {return THREE;}
    Real  is_inside(const coordinate_type& coord) const
    {
        throw ecell4::NotImplemented("polygon::is_inside");
    }
    coordinate_type draw_position(boost::shared_ptr<ecell4::RandomNumberGenerator>& rng) const
    {
        throw ecell4::NotImplemented("polygon::draw_position");
    }
    bool  test_AABB(const coordinate_type& l, const coordinate_type& u) const
    {
        throw ecell4::NotImplemented("polygon::test_AABB");
    }
};

template<typename coordT>
std::pair<std::vector<std::size_t>, std::pair<std::size_t, std::pair<Real, Real> > >
Polygon<coordT>::get_faces_within_radius(const coordinate_type& pos, const Real range) const
{
    std::vector<std::size_t> intruders;
    std::size_t nearest_idx = std::numeric_limits<std::size_t>::max();
    std::pair<Real, Real> min_dist(std::numeric_limits<Real>::max(), 0.);
    std::size_t idx = 0;
    for(typename std::vector<FaceTriangle<coordinate_type> >::const_iterator
            iter = faces.begin(); iter != faces.end(); ++iter)
    {
        std::pair<Real, Real> dist = distance(pos, *iter);
        if(dist.first <= range) // is intruder face
            intruders.push_back(idx);

        if(dist.first < min_dist.first) // is nearest one
        {
            min_dist = dist;
            nearest_idx = idx;
        }
        ++idx;
    }
    return std::make_pair(intruders, std::make_pair(nearest_idx, min_dist));
}

template<typename coordT>
std::pair<std::pair<coordT, coordT>, typename Polygon<coordT>::face_id_type>
apply_reflection(const Polygon<coordT>& poly,
                 const coordT& pos, const coordT& disp,
                 const typename Polygon<coordT>::face_id_type  intruder_face)
{
    using FaceID = typename Polygon<coordT>::face_id_type;

    const coordT stop = pos + disp;

    const std::pair<bool, coordT> test_result =
        test_intersect_segment_triangle(pos, stop, poly.faces.at(intruder_face));

    const coordT next_stop =
        reflect_plane(pos, stop, poly.faces.at(intruder_face));

    return std::make_pair(std::make_pair(test_result.second, next_stop),
                          intruder_face);
}

template<typename coordT>
std::pair<bool, std::pair<Real, typename Polygon<coordT>::face_id_type> >
intersect_ray(const Polygon<coordT>& poly,
              const coordT& pos, const coordT& disp,
              const typename Polygon<coordT>::face_id_type ignore_face)
{
    using FaceID = typename Polygon<coordT>::face_id_type;

    const coordT stop = pos + disp;
    const Real   len  = length(disp);

    bool   collide_face           = false;
    FaceID first_collide_face_idx = std::numeric_limits<std::size_t>::max();
    Real   first_collide_dist_sq  = len * len;

    for(const FaceID& intruder : poly.get_faces_within_radius(pos, len).first)
    {
        if(intruder == ignore_face)
        {
            continue;
        }
        const std::pair<bool, coordT> test_result =
            test_intersect_segment_triangle(pos, stop, poly.faces.at(intruder));

        if(test_result.first)
        {
            const Real distsq_to_face = length_sq(test_result.second - pos);
            if(distsq_to_face < first_collide_dist_sq)
            {
                collide_face           = true;
                first_collide_face_idx = intruder;
                first_collide_dist_sq  = distsq_to_face;
            }
        }
    }
    return std::make_pair(collide_face, std::make_pair(
                std::sqrt(first_collide_dist_sq), first_collide_face_idx));
}

} // egfrd
} // ecell4
#endif // ECELL4_EGFRD_POLYGON_HPP
