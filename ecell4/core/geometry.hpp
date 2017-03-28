#ifndef ECELL4_GEOMETRY
#define ECELL4_GEOMETRY
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Triangle.hpp>

namespace ecell4
{

Real3 rotate(const Real angle, const Real3& axis, const Real3& target);
Real  angle(const Real3& lhs, const Real3& rhs);

inline Real3 centroid(const Triangle& tri)
{
    return (tri.vertices()[0] + tri.vertices()[1] + tri.vertices()[2]) / 3.0;
}

inline Real3 incenter(const Triangle& tri)
{
    const valueT a = tri.length_of_edge()[1];
    const valueT b = tri.length_of_edge()[2];
    const valueT c = tri.length_of_edge()[0];
    const valueT abc = a + b + c;
    return (tri.vertices()[0] * a + tri.vertices()[1] * b + tri.vertices()[2] * c)
            * (1e0 / abc);
}

Real3 closest_point(const Real3&, const Triangle&);
Real  distance_sq(const Real3&, const Triangle&);
Real  distance(const Real3&, const Triangle&);

} // ecell4
#endif// ECELL4_GEOMETRY
