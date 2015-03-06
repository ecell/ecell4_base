#ifndef __ECELL4_VACANT_TYPE_HPP
#define __ECELL4_VACANT_TYPE_HPP

#include "MolecularTypeBase.hpp"

namespace ecell4
{

class VacantType
    : public MolecularTypeBase
{
public:

    typedef MolecularTypeBase base_type;
    typedef base_type::particle_info particle_info;
    typedef base_type::private_coordinate_type private_coordinate_type;
    typedef base_type::container_type container_type;
    typedef base_type::iterator iterator;
    typedef base_type::const_iterator const_iterator;

public:

    ~VacantType()
    {
        ; // do nothing
    }

    bool is_vacant() const
    {
        return true;
    }

    bool with_voxels() const
    {
        return false;
    }

    static VacantType& getInstance()
    {
        static VacantType instance;
        return instance;
    }

    virtual void add_voxel_without_checking(const particle_info& info)
    {
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

    bool remove_voxel_if_exists(const private_coordinate_type& coord)
    {
        return true; // just return true
    }

private:

    VacantType()
        : base_type(Species("VACANT", "0", "0"), NULL, 0, 0)
    {
        ;
    }
};

} // ecell4

#endif /* __ECELL4_VACANT_TYPE_HPP */
