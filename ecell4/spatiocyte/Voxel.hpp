#ifndef ECELL4_SPATIOCYTE_VOXEL_HPP
#define ECELL4_SPATIOCYTE_VOXEL_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/VoxelSpaceBase.hpp>
#include <boost/weak_ptr.hpp>

namespace ecell4
{

namespace spatiocyte
{

struct Voxel
{
    typedef Integer coordinate_type;

    Voxel(boost::weak_ptr<VoxelSpaceBase> space, coordinate_type coordinate)
        : space(space), coordinate(coordinate)
    {
    }

    boost::weak_ptr<VoxelSpaceBase> space;
    coordinate_type coordinate;

public:

    const Real3 position() const
    {
        return space.lock()->coordinate2position(coordinate);
    }

    bool clear() const
    {
        return space.lock()->remove_voxel(coordinate);
    }

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
        return Voxel(space, space.lock()->get_neighbor(coordinate, nrand));
    }
};

}

}

#endif
