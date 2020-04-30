#ifndef ECELL4_SPATIOCYTE_UTILS_HPP
#define ECELL4_SPATIOCYTE_UTILS_HPP

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

const Real calculate_dimensional_factor(
    std::shared_ptr<const VoxelPool> mt0, const Real D_A,
    std::shared_ptr<const VoxelPool> mt1, const Real D_B,
    std::shared_ptr<SpatiocyteWorld> world);

const Real calculate_alpha(const ReactionRule &rr,
                           const std::shared_ptr<SpatiocyteWorld> &world);

} // namespace spatiocyte

} // namespace ecell4

#endif /* ECELL4_SPATIOCYTE_UTILS_HPP */
