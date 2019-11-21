#ifndef ECELL4_VOXELVIEW_HPP
#define ECELL4_VOXELVIEW_HPP

#include "Identifier.hpp"
#include "Species.hpp"
#include "types.hpp"

namespace ecell4
{

template <typename T>
struct ParticleBase
{
    ParticleID pid;
    const Species &species;
    T voxel;

    ParticleBase(ParticleID pid, const Species &species, T voxel)
        : pid(pid), species(species), voxel(voxel)
    {
    }
};

using VoxelView = ParticleBase<Integer>;

} // namespace ecell4

#endif /* ECELL4_VOXELVIEW_HPP */
