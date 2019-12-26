#ifndef ECELL4_GEOMETRY_HPP
#define ECELL4_GEOMETRY_HPP
#include <ecell4/core/Real3.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/algorithm/clamp.hpp>
#include <cmath>

namespace ecell4
{

inline Real calc_angle(const Real3& lhs, const Real3& rhs)
{
    const Real cosine = dot_product(lhs, rhs) /
                        std::sqrt(length_sq(lhs) * length_sq(rhs));
    return std::acos(boost::algorithm::clamp(cosine, -1.0, 1.0));
}

Real3 rotate(const Real angle, const Real3& axis, const Real3& target);

} // ecell4
#endif// ECELL4_GEOMETRY_HPP
