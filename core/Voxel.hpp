#ifndef __ECELL4_VOXEL_HPP
#define __ECELL4_VOXEL_HPP

#include "MolecularType.hpp" //XXX: needed only for the Coord definition.

namespace ecell4
{

class Voxel
{
public:

    Voxel()
    {
        ;
    }

    Voxel(const Species& sp, const Coord& coord, const Real& D)
        : species_(sp), coordinate_(coord), D_(D) {}

    const Species& species() const
    {
        return species_;
    }

    Species& species()
    {
        return species_;
    }

    const Coord& coordinate() const
    {
        return coordinate_;
    }

    Coord& coordinate()
    {
        return coordinate_;
    }

    const Real& D() const
    {
        return D_;
    }

    Real& D()
    {
        return D_;
    }

private:

    Species species_;
    Coord coordinate_;
    Real D_;
};

}

#endif
