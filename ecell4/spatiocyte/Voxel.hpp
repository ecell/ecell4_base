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
};

}

}

#endif
