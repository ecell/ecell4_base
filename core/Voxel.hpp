#ifndef __ECELL4_VOXEL_HPP
#define __ECELL4_VOXEL_HPP

namespace ecell4
{

class Voxel
{
public:

    typedef Integer coordinate_type;
    // typedef LatticeSpace::coordinate_type coordinate_type;

public:

    Voxel()
    {
        ;
    }

    Voxel(const Species& sp, const coordinate_type& coord, const Real& D)
        : species_(sp), coordinate_(coord), D_(D) {}

    const Species& species() const
    {
        return species_;
    }

    Species& species()
    {
        return species_;
    }

    const coordinate_type& coordinate() const
    {
        return coordinate_;
    }

    coordinate_type& coordinate()
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
    coordinate_type coordinate_;
    Real D_;
};

}

#endif
