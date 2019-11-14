#ifndef ECELL4_SPATIOCYTE_UTILS_HPP
#define ECELL4_SPATIOCYTE_UTILS_HPP

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

const Real
calculate_dimensional_factor(boost::shared_ptr<const VoxelPool> mt0,
                             boost::shared_ptr<const VoxelPool> mt1,
                             boost::shared_ptr<SpatiocyteWorld> world);

const Real calculate_alpha(const ReactionRule &rr,
                           const boost::shared_ptr<SpatiocyteWorld> &world);

} // namespace spatiocyte

} // namespace ecell4

#endif /* ECELL4_SPATIOCYTE_UTILS_HPP */
