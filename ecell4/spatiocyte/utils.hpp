#ifndef ECELL4_SPATIOCYTE_UTILS_HPP
#define ECELL4_SPATIOCYTE_UTILS_HPP

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

const Real calculate_dimensional_factor(
    const VoxelPool* mt0, const VoxelPool* mt1,
    const boost::shared_ptr<const SpatiocyteWorld>& world);

const Real calculate_alpha(
    const ReactionRule& rr, const boost::shared_ptr<SpatiocyteWorld>& world);

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_UTILS_HPP */
