#ifndef ECELL4_SPATIOCYTE_VOXEL_HPP
#define ECELL4_SPATIOCYTE_VOXEL_HPP

#include <boost/weak_ptr.hpp>
#include <ecell4/core/VoxelSpaceBase.hpp>
#include <ecell4/core/types.hpp>

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteWorld;

struct Voxel
{
    typedef Integer coordinate_type;

    Voxel(const SpatiocyteWorld *world, boost::weak_ptr<VoxelSpaceBase> space,
          coordinate_type coordinate)
        : world(world), space(space), coordinate(coordinate)
    {
    }

    const SpatiocyteWorld *world;
    boost::weak_ptr<VoxelSpaceBase> space;
    coordinate_type coordinate;

public:
    const Real3 position() const
    {
        return space.lock()->coordinate2position(coordinate);
    }

    bool clear() const { return space.lock()->remove_voxel(coordinate); }

    boost::shared_ptr<VoxelPool> get_voxel_pool() const
    {
        return space.lock()->get_voxel_pool_at(coordinate);
    }

    Integer num_neighbors() const
    {
        return space.lock()->num_neighbors(coordinate);
    }

    Voxel get_neighbor(Integer nrand) const
    {
        return Voxel(world, space,
                     space.lock()->get_neighbor(coordinate, nrand));
    }

    Voxel get_neighbor_randomly(
        const boost::shared_ptr<RandomNumberGenerator> &rng) const
    {
        const Integer idx(rng->uniform_int(0, num_neighbors() - 1));
        return get_neighbor(idx);
    }

    Voxel
    get_neighbor_randomly(const boost::shared_ptr<RandomNumberGenerator> &rng,
                          Shape::dimension_kind dimension) const;
};

} // namespace spatiocyte

} // namespace ecell4

#endif
