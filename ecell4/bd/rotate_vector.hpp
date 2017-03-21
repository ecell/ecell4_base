#ifndef ECELL4_VECTOR_ROTATION
#define ECELL4_VECTOR_ROTATION
#include <ecell4/core/Real3.hpp>
#include <cmath>

namespace ecell4
{

Real3 rotate(const Real angle, const Real3& axis, const Real3& target);
inline Real angle(const Real3& lhs, const Real3& rhs)
{
    const Real costheta = dot_product(lhs, rhs) / std::sqrt(length_sq(lhs) * length_sq(rhs));
    if(1.0 < costheta)
    {
        if(std::abs(costheta - 1.0) < 1e-8)
            return std::acos(1.0);
        else
            throw std::out_of_range("acos");
    }
    else if(costheta < -1.0)
    {
        if(std::abs(costheta + 1.0) < 1e-8)
            return std::acos(-1.0);
        else
            throw std::out_of_range("acos");
    }
    else
        return std::acos(costheta);
}


}// ecell4
#endif /* ECELL_VECTOR_ROTATION */
