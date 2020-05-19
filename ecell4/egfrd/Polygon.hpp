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
template<typename coordT>
coordT closest_point(const coordT& pos, const std::array<coordT, 3>& vertices)
{
    typedef typename element_type_of<coordT>::type valueT;
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.141-142

    const coordT a = vertices[0];
    const coordT b = vertices[1];
    const coordT c = vertices[2];

    const coordT ab = b - a;
    const coordT ac = c - a;
    const coordT ap = pos - a;
    const valueT d1 = dot_product(ab, ap);
    const valueT d2 = dot_product(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        return a;
    }

    const coordT bp = pos - b;
    const valueT d3 = dot_product(ab, bp);
    const valueT d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
        return b;
    }

    const valueT vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        valueT v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const coordT cp = pos - c;
    const valueT d5 = dot_product(ab, cp);
    const valueT d6 = dot_product(ac, cp);
    if (d6 >= 0.0 && d5 <= d6)
    {
        return c;
    }

    const valueT vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const valueT w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const valueT va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const valueT w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    const valueT denom = 1.0 / (va + vb + vc);
    const valueT v = vb * denom;
    const valueT w = vc * denom;
    return a + ab * v + ac * w;
}

template<typename coordT>
std::pair<bool, coordT>
test_intersect_segment_triangle(const coordT& begin, const coordT& end,
                                const std::array<coordT, 3>& vertices)
{
    typedef typename element_type_of<coordT>::type valueT;
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.190-194

    const coordT line = begin - end;
    const coordT ab = vertices[1] - vertices[0];
    const coordT ac = vertices[2] - vertices[0];
    const coordT normal = cross_product(ab, ac);

    const valueT d = dot_product(line, normal);
    if(d < 0.0)
        return std::make_pair(false, coordT(0.,0.,0.));

    const coordT ap = begin - vertices[0];
    const valueT t = dot_product(ap, normal);
    if(t < 0.0 || d < t)
        return std::make_pair(false, coordT(0.,0.,0.));

    const coordT e = cross_product(line, ap);
    valueT v = dot_product(ac, e);
    if(v < 0. || d < v)
        return std::make_pair(false, coordT(0.,0.,0.));
    valueT w = -1.0 * dot_product(ab, e);
    if(w < 0. || d < v + w)
        return std::make_pair(false, coordT(0.,0.,0.));

    const valueT ood = 1. / d;
    v *= ood;
    w *= ood;
    const valueT u = 1. - v - w;
    const coordT intersect = vertices[0] * u + vertices[1] * v + vertices[2] * w;

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

inline Real3 reflect_plane(const Real3& begin, const Real3& end, const Triangle& face)
{
    return reflect_plane(begin, end, face.normal(), face.vertex_at(0));
}

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
