#ifndef ECELL4_SPATIOCYTE_OFFLATTICE_HPP
#define ECELL4_SPATIOCYTE_OFFLATTICE_HPP

#include <limits>
#include <ecell4/core/OffLatticeSpace.hpp>

namespace ecell4
{

namespace spatiocyte {

class OffLattice {
public:
    using positions_type = OffLatticeSpace::position_container;
    using adjoining_pairs_type = OffLatticeSpace::coordinate_pair_list_type;

    OffLattice(const Real voxel_radius, const positions_type positions, const adjoining_pairs_type adjoining_pairs)
        : voxel_radius_(voxel_radius), positions_(positions), adjoining_pairs_(adjoining_pairs) {}

    OffLattice(const Real voxel_radius, const positions_type positions)
        : voxel_radius_(voxel_radius), positions_(positions)
    {
        constexpr Real epsilon = std::numeric_limits<Real>::epsilon();

        adjoining_pairs_.clear();
        const std::size_t size = positions_.size();
        for (std::size_t i(0); i < size; ++i)
            for (std::size_t j(i+1); j < size; ++j)
                if (length(positions_[j] - positions_[i]) <= voxel_radius_ + epsilon)
                {
                    adjoining_pairs_.push_back(std::make_pair(i, j));
                }
    }

    inline const Real& voxel_radius() const
    {
        return voxel_radius_;
    }

    inline const positions_type& positions() const
    {
        return positions_;
    }

    inline const adjoining_pairs_type& adjoining_pairs() const
    {
        return adjoining_pairs_;
    }

    std::unique_ptr<OffLatticeSpace> generate_space(const Species& species) const
    {
        return std::unique_ptr<OffLatticeSpace>(
                new OffLatticeSpace(voxel_radius_, species, positions_, adjoining_pairs_));
    }

protected:

    Real voxel_radius_;
    positions_type positions_;
    adjoining_pairs_type adjoining_pairs_;
};

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_OFFLATTICE_HPP */
