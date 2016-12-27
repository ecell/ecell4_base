#include "rotate_vector.hpp"
#include <boost/math/quaternion.hpp>

namespace ecell4
{

Real3 rotate(const Real angle, const Real3& axis, const Real3& target)
{
    typedef boost::math::quaternion<Real> Quaternion;
    const Real cos_t(cos(angle / 2));
    const Real sin_t(sin(angle / 2));
    const Real sin_normalize(sin_t / length(axis));

    const Quaternion Q(cos_t, axis[0] * sin_normalize, 
                              axis[1] * sin_normalize,
                              axis[2] * sin_normalize);
    const Quaternion P(0e0, target[0], target[1], target[2]);
    const Quaternion S(Q * P * conj(Q));

    return Real3(S.R_component_2(), S.R_component_3(), S.R_component_4());
}

} // ecell4
