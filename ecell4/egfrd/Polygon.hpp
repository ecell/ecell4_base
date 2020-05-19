#ifndef ECELL4_EGFRD_POLYGON_HPP
#define ECELL4_EGFRD_POLYGON_HPP

#include <ecell4/core/Shape.hpp>
#include <ecell4/core/AABB.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/comparators.hpp>
#include "exceptions.hpp"
#include "FaceTriangle.hpp"
#include <limits>
#include <algorithm>

namespace ecell4
{
namespace egfrd
{

inline std::pair<std::pair<Real3, Real3>, Polygon::FaceID>
apply_reflection(const Polygon& poly, const Real3& pos, const Real3& disp,
                 const Polygon::FaceID intruder_face)
{
    using FaceID = Polygon::FaceID;

    const Real3 stop = pos + disp;

    const std::pair<bool, Real3> test_result =
        test_intersect_segment_triangle(pos, stop, poly.triangle_at(intruder_face));

    const Real3 next_stop =
        reflect_plane(pos, stop, poly.triangle_at(intruder_face));

    return std::make_pair(std::make_pair(test_result.second, next_stop),
                          intruder_face);
}

inline std::pair<bool, std::pair<Real, boost::optional<Polygon::FaceID>>>
intersect_ray(const Polygon& poly, const Real3& pos, const Real3& disp,
              const boost::optional<Polygon::FaceID> ignore_face)
{
    using FaceID = Polygon::FaceID;

    const Real3 stop = pos + disp;
    const Real   len  = length(disp);

    bool collide_face          = false;
    Real first_collide_dist_sq = len * len;
    boost::optional<FaceID> first_collide_face_idx = boost::none;

    const auto intruders =
        ignore_face ? poly.list_faces_within_radius(pos, len, *ignore_face) :
                      poly.list_faces_within_radius(pos, len);

    for(const auto& intruder : intruders)
    {
        const FaceID&   fid = intruder.first.first;
        const Triangle& tri = intruder.first.second;

        const std::pair<bool, Real3> test_result =
            test_intersect_segment_triangle(pos, stop, tri);

        if(test_result.first)
        {
            const Real distsq_to_face = length_sq(test_result.second - pos);
            if(distsq_to_face < first_collide_dist_sq)
            {
                collide_face           = true;
                first_collide_face_idx = fid;
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
