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
struct Polygon : public ecell4::Shape
{
    typedef FaceTriangle<Real3> face_type;
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
    get_faces_within_radius(const Real3& pos, const Real range) const
    {
        std::vector<std::size_t> intruders;
        std::size_t nearest_idx = std::numeric_limits<std::size_t>::max();
        std::pair<Real, Real> min_dist(std::numeric_limits<Real>::max(), 0.0);

        std::size_t idx = 0;
        for(const auto& face : this->faces)
        {
            std::pair<Real, Real> dist = distance(pos, face);
            if(dist.first <= range) // is intruder face
            {
                intruders.push_back(idx);
            }

            if(dist.first < min_dist.first) // is nearest one
            {
                min_dist = dist;
                nearest_idx = idx;
            }
            ++idx;
        }
        return std::make_pair(intruders, std::make_pair(nearest_idx, min_dist));
    }

    void emplace(const std::array<Real3, 3>& vertices)
    {
        this->faces.push_back(face_type(vertices));
    }

// data member
    std::vector<face_type> faces;

// for shapes (not implemented yet)
    dimension_kind dimension() const {return THREE;}
    Real  is_inside(const Real3& coord) const
    {
        throw ecell4::NotImplemented("polygon::is_inside");
    }
    Real3 draw_position(std::shared_ptr<ecell4::RandomNumberGenerator>& rng) const
    {
        throw ecell4::NotImplemented("polygon::draw_position");
    }
    bool  test_AABB(const Real3& l, const Real3& u) const
    {
        throw ecell4::NotImplemented("polygon::test_AABB");
    }
};

inline std::pair<std::pair<Real3, Real3>, Polygon::face_id_type>
apply_reflection(const Polygon& poly, const Real3& pos, const Real3& disp,
                 const Polygon::face_id_type intruder_face)
{
    using FaceID = Polygon::face_id_type;

    const Real3 stop = pos + disp;

    const std::pair<bool, Real3> test_result =
        test_intersect_segment_triangle(pos, stop, poly.faces.at(intruder_face));

    const Real3 next_stop =
        reflect_plane(pos, stop, poly.faces.at(intruder_face));

    return std::make_pair(std::make_pair(test_result.second, next_stop),
                          intruder_face);
}

inline std::pair<bool, std::pair<Real, Polygon::face_id_type> >
intersect_ray(const Polygon& poly, const Real3& pos, const Real3& disp,
              const Polygon::face_id_type ignore_face)
{
    using FaceID = Polygon::face_id_type;

    const Real3 stop = pos + disp;
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
        const std::pair<bool, Real3> test_result =
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
