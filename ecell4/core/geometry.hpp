#ifndef ECELL4_GEOMETRY
#define ECELL4_GEOMETRY
#include <ecell4/core/Real3.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>

namespace ecell4
{

inline Real angle(const Real3& lhs, const Real3& rhs)
{
    const Real cosine  = dot_product(lhs, rhs) /
                         std::sqrt(length_sq(lhs) * length_sq(rhs));
         if(cosine < -1.){return boost::math::constants::pi<Real>();}
    else if(cosine >  1.){return 0;}
    else                 {return std::acos(cosine);}
}

Real3 rotate(const Real angle, const Real3& axis, const Real3& target);

} // ecell4
#endif// ECELL4_GEOMETRY
