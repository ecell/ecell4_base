#ifndef ECELL4_VOXEL_HPP
#define ECELL4_VOXEL_HPP

namespace ecell4
{

struct ParticleVoxel
{
    typedef Integer coordinate_type;

    ParticleVoxel()
    {
        ;
    }

    ParticleVoxel(const Species& sp,
                  const coordinate_type& coord,
                  const Real& radius,
                  const Real& D,
                  const std::string& loc = "")
        : species(sp), coordinate(coord), radius(radius), D(D), loc(loc)
    {}

    Species species;
    coordinate_type coordinate;
    Real radius;
    Real D;
    std::string loc;
};

}

#endif
