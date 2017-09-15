#ifndef ECELL4_VOXEL_HPP
#define ECELL4_VOXEL_HPP

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

    Voxel(const Species& sp, const coordinate_type& coord,
        const Real& radius, const Real& D, const std::string& loc = "")
        : species_(sp), coordinate_(coord), radius_(radius), D_(D), loc_(loc) {}

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

    const Real& radius() const
    {
        return radius_;
    }

    Real& radius()
    {
        return radius_;
    }

    const std::string& loc() const
    {
        return loc_;
    }

    std::string& loc()
    {
        return loc_;
    }

private:

    Species species_;
    coordinate_type coordinate_;
    Real radius_;
    Real D_;
    std::string loc_;
};

}

#endif
