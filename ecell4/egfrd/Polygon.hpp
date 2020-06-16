#ifndef ECELL4_EGFRD_POLYGON_HPP
#define ECELL4_EGFRD_POLYGON_HPP

#include <ecell4/core/geometry.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/AABB.hpp>
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/comparators.hpp>
#include "exceptions.hpp"
#include <limits>
#include <algorithm>

namespace ecell4
{
namespace egfrd
{

namespace detail
{
inline Real3 closest_point(const Real3& pos, const std::array<Real3, 3>& vertices)
{
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.141-142

    const Real3 a = vertices[0];
    const Real3 b = vertices[1];
    const Real3 c = vertices[2];

    const Real3 ab = b - a;
    const Real3 ac = c - a;
    const Real3 ap = pos - a;
    const Real d1 = dot_product(ab, ap);
    const Real d2 = dot_product(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        return a;
    }

    const Real3 bp = pos - b;
    const Real d3 = dot_product(ab, bp);
    const Real d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
        return b;
    }

    const Real vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        Real v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const Real3 cp = pos - c;
    const Real d5 = dot_product(ab, cp);
    const Real d6 = dot_product(ac, cp);
    if (d6 >= 0.0 && d5 <= d6)
    {
        return c;
    }

    const Real vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const Real w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const Real va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    const Real denom = 1.0 / (va + vb + vc);
    const Real v = vb * denom;
    const Real w = vc * denom;
    return a + ab * v + ac * w;
}

inline std::pair<bool, Real3>
test_intersect_segment_triangle(const Real3& begin, const Real3& end,
                                const std::array<Real3, 3>& vertices)
{
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.190-194

    const Real3 line = begin - end;
    const Real3 ab = vertices[1] - vertices[0];
    const Real3 ac = vertices[2] - vertices[0];
    const Real3 normal = cross_product(ab, ac);

    const Real d = dot_product(line, normal);
    if(d < 0.0)
    {
        return std::make_pair(false, Real3(0.,0.,0.));
    }

    const Real3 ap = begin - vertices[0];
    const Real t = dot_product(ap, normal);
    if(t < 0.0 || d < t)
    {
        return std::make_pair(false, Real3(0.,0.,0.));
    }

    const Real3 e = cross_product(line, ap);
    Real v = dot_product(ac, e);
    if(v < 0. || d < v)
    {
        return std::make_pair(false, Real3(0.,0.,0.));
    }
    Real w = -1.0 * dot_product(ab, e);
    if(w < 0. || d < v + w)
    {
        return std::make_pair(false, Real3(0.,0.,0.));
    }

    const Real ood = 1. / d;
    v *= ood;
    w *= ood;
    const Real u = 1. - v - w;
    const Real3 intersect = vertices[0] * u + vertices[1] * v + vertices[2] * w;

    return std::make_pair(true, intersect);
}
} // detail

inline Real distance_point_to_triangle(const Real3& pos, const Triangle& face)
{
    std::array<Real3, 3> triangle = face.vertices();
    if(dot_product(pos - face.vertex_at(0), face.normal()) < 0)
    {
        triangle[0] = face.vertex_at(2);
        triangle[1] = face.vertex_at(1);
        triangle[2] = face.vertex_at(0);
    }
    return length(detail::closest_point(pos, triangle) - pos);
}

inline std::pair<bool, Real3>
test_intersect_segment_triangle(const Real3& begin, const Real3& end,
                                const Triangle& face)
{
    const Real3 line = end - begin;
    if(dot_product(line, face.normal()) < 0.0)
    {
        return detail::test_intersect_segment_triangle(begin, end, face.vertices());
    }
    else
    {
        std::array<Real3, 3> rev;
        rev[0] = face.vertex_at(2);
        rev[1] = face.vertex_at(1);
        rev[2] = face.vertex_at(0);
        return detail::test_intersect_segment_triangle(begin, end, rev);
    }
}

inline std::pair<std::pair<Real3, Real3>, FaceID>
apply_reflection(const Polygon& poly, const Real3& pos, const Real3& disp,
                 const FaceID intruder_face)
{
    const Real3 stop = pos + disp;
    const auto& tri  = poly.triangle_at(intruder_face);

    const std::pair<bool, Real3> test_result =
        test_intersect_segment_triangle(pos, stop, tri);

    const Real3 next_stop =
        reflect_plane(pos, stop, tri.normal(), tri.vertex_at(0));

    return std::make_pair(std::make_pair(test_result.second, next_stop),
                          intruder_face);
}

inline std::pair<bool, std::pair<Real, boost::optional<FaceID>>>
intersect_ray(const Polygon& poly, const Real3& pos, const Real3& disp,
              const boost::optional<FaceID> ignore_face)
{
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
