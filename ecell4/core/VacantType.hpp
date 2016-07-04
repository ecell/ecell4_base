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
    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::coordinate_type coordinate_type;
    typedef base_type::container_type container_type;
    typedef base_type::iterator iterator;
    typedef base_type::const_iterator const_iterator;
    typedef base_type::voxel_type_type voxel_type_type;

public:

    ~VacantType()
    {
        ; // do nothing
    }

    virtual voxel_type_type const voxel_type() const
    {
        return VACANT;
    }

    // bool is_vacant() const
    // {
    //     return true;
    // }

    bool with_voxels() const
    {
        return false;
    }

    static VacantType& getInstance()
    {
        static VacantType instance;
        return instance;
    }

    const Shape::dimension_kind get_dimension() const
    {
        return Shape::THREE;
    }

    virtual void add_voxel_without_checking(const coordinate_id_pair_type& info)
    {
        ; // do nothing
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord,
        const coordinate_id_pair_type& to_info)
    {
        ; // do nothing
    }

    virtual void replace_voxel(
        const coordinate_type& from_coord,
        const coordinate_type& to_coord,
        const std::size_t candidate=0)
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
    }
};

} // ecell4

#endif /* __ECELL4_VACANT_TYPE_HPP */
