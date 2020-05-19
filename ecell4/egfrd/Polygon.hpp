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
struct Polygon : public ecell4::Shape
{
    typedef Triangle face_type;
    typedef std::size_t FaceID;

    Polygon(){}
    ~Polygon(){}

    std::vector<std::pair<std::pair<FaceID, Triangle>, Real>>
    list_faces_within_radius(const Real3& pos, const Real range) const
    {
        std::vector<std::pair<std::pair<FaceID, Triangle>, Real>> intruders;

        for(std::size_t i=0; i<this->faces.size(); ++i)
        {
            const auto& face = this->faces.at(i);
            const Real dist = distance_point_to_triangle(pos, face);
            if(dist <= range) // is intruder face
            {
                intruders.emplace_back(std::make_pair(i, face), dist);
            }
        }
        std::sort(intruders.begin(), intruders.end(),
            utils::pair_second_element_comparator<
                std::pair<FaceID, Triangle>, Real>{});

        return intruders;
    }
    std::vector<std::pair<std::pair<FaceID, Triangle>, Real>>
    list_faces_within_radius(const Real3& pos, const Real range, const FaceID& ignore) const
    {
        std::vector<std::pair<std::pair<FaceID, Triangle>, Real>>
            intruders;

        for(std::size_t i=0; i<this->faces.size(); ++i)
        {
            if(i == ignore)
            {
                continue;
            }
            const auto& face = this->faces.at(i);
            const Real dist = distance_point_to_triangle(pos, face);
            if(dist <= range) // is intruder face
            {
                intruders.emplace_back(std::make_pair(i, face), dist);
            }
        }
        std::sort(intruders.begin(), intruders.end(),
            utils::pair_second_element_comparator<
                std::pair<FaceID, Triangle>, Real>{});

        return intruders;
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

inline std::pair<std::pair<Real3, Real3>, Polygon::FaceID>
apply_reflection(const Polygon& poly, const Real3& pos, const Real3& disp,
                 const Polygon::FaceID intruder_face)
{
    using FaceID = Polygon::FaceID;

    const Real3 stop = pos + disp;

    const std::pair<bool, Real3> test_result =
        test_intersect_segment_triangle(pos, stop, poly.faces.at(intruder_face));

    const Real3 next_stop =
        reflect_plane(pos, stop, poly.faces.at(intruder_face));

    return std::make_pair(std::make_pair(test_result.second, next_stop),
                          intruder_face);
}

inline std::pair<bool, std::pair<Real, Polygon::FaceID> >
intersect_ray(const Polygon& poly, const Real3& pos, const Real3& disp,
              const boost::optional<Polygon::FaceID> ignore_face)
{
    using FaceID = Polygon::FaceID;

    const Real3 stop = pos + disp;
    const Real   len  = length(disp);

    bool   collide_face           = false;
    FaceID first_collide_face_idx = std::numeric_limits<std::size_t>::max();
    Real   first_collide_dist_sq  = len * len;

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
