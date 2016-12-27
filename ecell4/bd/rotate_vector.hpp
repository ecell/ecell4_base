#ifndef ECELL4_VECTOR_ROTATION
#define ECELL4_VECTOR_ROTATION
#include <ecell4/core/Real3.hpp>

namespace ecell4
{

Real3 rotate(const Real angle, const Real3& axis, const Real3& target);
inline Real angle(const Real3& lhs, const Real3& rhs)
{
    return acos(dot_product(lhs, rhs) / std::sqrt(length_sq(lhs) * length_sq(rhs)));
}


}// ecell4
#endif /* ECELL_VECTOR_ROTATION */
