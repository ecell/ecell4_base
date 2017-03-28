#ifndef ECELL4_GEOMETRY
#define ECELL4_GEOMETRY
#include <ecell4/core/Real3.hpp>

namespace ecell4
{

Real3 rotate(const Real angle, const Real3& axis, const Real3& target);
Real  angle(const Real3& lhs, const Real3& rhs);

} // ecell4
#endif// ECELL4_GEOMETRY
