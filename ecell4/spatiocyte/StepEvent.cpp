#include "SpatiocyteEvent.hpp"

namespace ecell4
{

namespace spatiocyte
{

template <>
const Real calc_dt<3>(const Real R, const Real D)
{
    return 2 * R * R / 3 / D;
}

template <>
const Real calc_dt<2>(const Real R, const Real D)
{
    return R * R / D;
}

} // namespace spatiocyte

} // namespace ecell4
