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

public:

    ~VacantType()
    {
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
        private_coordinate_type from_coord,
        const particle_info& to_info)
    {
        ; // do nothing
    }

    virtual void remove_voxel(const container_type::iterator& position)
    {
        ; // do nothing
    }

    virtual void swap(
        const container_type::iterator& a, const container_type::iterator& b)
    {
        ; // do nothing
    }

    void addVoxel(particle_info info)
    {
        ; // do nothing
    }

    bool removeVoxel(private_coordinate_type coord)
    {
        return true; // just return true
    }

    bool is_vacant() const
    {
        return true;
    }

private:

    VacantType()
        : MolecularTypeBase(Species("VACANT", "0"), NULL, 0, 0)
    {
        ;
    }
};

} // ecell4

#endif /* __ECELL4_VACANT_TYPE_HPP */
