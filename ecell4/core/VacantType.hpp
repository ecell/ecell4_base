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
    typedef base_type::coord_id_pair coord_id_pair;
    typedef base_type::coordinate_type coordinate_type;
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

    virtual void add_voxel_without_checking(const coord_id_pair& info)
    {
        ; // do nothing
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord,
        const coord_id_pair& to_info)
    {
        ; // do nothing
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord,
        const coordinate_type& to_coord)
    {
        ; // do nothing
    }

    bool remove_voxel_if_exists(const coordinate_type& coord)
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
