#ifndef ECELL4_SPATIOCYTE_VOXEL_HPP
#define ECELL4_SPATIOCYTE_VOXEL_HPP

#include <ecell4/core/types.hpp>

namespace ecell4
{

namespace spatiocyte
{

struct Voxel
{
    typedef Integer coordinate_type;

    Voxel(coordinate_type coordinate)
        : coordinate(coordinate)
    {
    }

    coordinate_type coordinate;
};

}

}

#endif
