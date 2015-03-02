#ifndef __ECELL4_STRUCTURE_TYPE_HPP
#define __ECELL4_STRUCTURE_TYPE_HPP

#include "MolecularType.hpp"


namespace ecell4
{

class StructureType
    : public MolecularType
{
public:

    typedef MolecularType base_type;
    typedef base_type::particle_info particle_info;
    typedef base_type::private_coordinate_type private_coordinate_type;
    typedef base_type::container_type container_type;
    typedef base_type::iterator iterator;
    typedef base_type::const_iterator const_iterator;

public:

    StructureType(
        const Species& species, MolecularTypeBase* location,
        const Real& radius = 0.0)
        : base_type(species, location, radius, 0)
        // : base_type(species, &(VacantType::getInstance()), radius, 0)
    {
        ;
    }

    virtual ~StructureType()
    {
        ;
    }

    // bool is_vacant() const
    // {
    //     return false;
    // }

    bool with_voxels() const
    {
        return false;
    }

    virtual void add_voxel_without_checking(const particle_info& info)
    {
        if (info.second != ParticleID())
        {
            throw NotSupported("No ParticleID is allowed.");
        }

        ; // do nothing
    }

    virtual void replace_voxel(
        const private_coordinate_type& from_coord,
        const particle_info& to_info)
    {
        ; // do nothing
    }

    virtual void replace_voxel(
        const private_coordinate_type& from_coord,
        const private_coordinate_type& to_coord)
    {
        ; // do nothing
    }

    virtual particle_info pop(const private_coordinate_type& coord)
    {
        return particle_info(coord, ParticleID());
    }

    virtual bool remove_voxel_if_exists(const private_coordinate_type& coord)
    {
        return true;
    }
};

} //ecell4

#endif /* __ECELL4_STRUCTURE_TYPE_HPP */
