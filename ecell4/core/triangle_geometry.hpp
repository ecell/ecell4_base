#ifndef ECELL4_TRIANGLE_GEOMETRY
#define ECELL4_TRIANGLE_GEOMETRY
#include <ecell4/core/geometry.hpp>
#include <ecell4/core/Triangle.hpp>
#include <ecell4/core/TriangleView.hpp>

namespace ecell4
{

inline Real3 centroid(const Triangle& tri)
{
    return (tri.vertices()[0] + tri.vertices()[1] + tri.vertices()[2]) * (1. / 3.);
}
inline Real3 centroid(const TriangleView& tri)
{
    return (tri.vertices(0) + tri.vertices(1) + tri.vertices(2)) * (1. / 3.);
}
inline Real3 centroid(const TriangleConstView& tri)
{
    return (tri.vertices(0) + tri.vertices(1) + tri.vertices(2)) * (1. / 3.);
}

inline Real3 incenter(const Triangle& tri)
{
    const Real a = tri.lengths_of_edges()[1];
    const Real b = tri.lengths_of_edges()[2];
    const Real c = tri.lengths_of_edges()[0];
    const Real rabc = 1. / (a + b + c);
    return (tri.vertices()[0] * a + tri.vertices()[1] * b + tri.vertices()[2] * c) * rabc;
}

inline Real3 incenter(const TriangleView& tri)
{
    const Real a = length(tri.vertices(0));
    const Real b = length(tri.vertices(1));
    const Real c = length(tri.vertices(2));
    const Real rabc = 1. / (a + b + c);
    return (tri.vertices(0) * a + tri.vertices(1) * b + tri.vertices(2) * c) * rabc;
}

inline Real3 incenter(const TriangleConstView& tri)
{
    const Real a = length(tri.vertices(0));
    const Real b = length(tri.vertices(1));
    const Real c = length(tri.vertices(2));
    const Real rabc = 1. / (a + b + c);
    return (tri.vertices(0) * a + tri.vertices(1) * b + tri.vertices(2) * c) * rabc;
}

Real3 closest_point(const Real3&, const Triangle&);
Real3 closest_point(const Real3&, const TriangleView&);
Real3 closest_point(const Real3&, const TriangleConstView&);

inline Real distance_sq(const Real3& pos, const Triangle& tri)
{
    return length_sq(closest_point(pos, tri) - pos);
}

inline Real distance(const Real3& pos, const Triangle& tri)
{
    return std::sqrt(distance_sq(pos, tri));
}

inline Real distance_sq(const Real3& pos, const TriangleView& tri)
{
    return length_sq(closest_point(pos, tri) - pos);
}

inline Real distance(const Real3& pos, const TriangleView& tri)
{
    return std::sqrt(distance_sq(pos, tri));
}

inline Real distance_sq(const Real3& pos, const TriangleConstView& tri)
{
    return length_sq(closest_point(pos, tri) - pos);
}

inline Real distance(const Real3& pos, const TriangleConstView& tri)
{
    return std::sqrt(distance_sq(pos, tri));
}

} // ecell4
#endif// ECELL4_TRIANGLE_GEOMETRY
