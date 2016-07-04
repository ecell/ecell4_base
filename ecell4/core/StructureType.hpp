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
    typedef base_type::coord_id_pair coord_id_pair;
    typedef base_type::private_coordinate_type private_coordinate_type;
    typedef base_type::container_type container_type;
    typedef base_type::iterator iterator;
    typedef base_type::const_iterator const_iterator;
    typedef base_type::voxel_type_type voxel_type_type;

public:

    StructureType(
        const Species& species, MolecularTypeBase* location,
        const Real& radius = 0.0, const Shape::dimension_kind& dimension=Shape::UNDEF)
        : base_type(species, location, radius, 0),
        dimension_(std::min(dimension, location->get_dimension()))
        // : base_type(species, &(VacantType::getInstance()), radius, 0)
    {
        ;
    }

    virtual ~StructureType()
    {
        ;
    }

    virtual voxel_type_type const voxel_type() const
    {
        return STRUCTURE;
    }

    // bool is_vacant() const
    // {
    //     return false;
    // }

    bool with_voxels() const
    {
        return false;
    }

    // bool is_structure() const
    // {
    //     return true;
    // }

    const Shape::dimension_kind get_dimension() const
    {
        return dimension_;
    }

    virtual void add_voxel_without_checking(const coord_id_pair& info)
    {
        if (info.second != ParticleID())
        {
            throw NotSupported("No ParticleID is allowed.");
        }

        ; // do nothing
    }

    virtual void replace_voxel(
        const private_coordinate_type& from_coord,
        const coord_id_pair& to_info)
    {
        ; // do nothing
    }

    virtual void replace_voxel(
        const private_coordinate_type& from_coord,
        const private_coordinate_type& to_coord,
        const std::size_t candidate=0)
    {
        ; // do nothing
    }

    virtual coord_id_pair pop(const private_coordinate_type& coord)
    {
        return coord_id_pair(coord, ParticleID());
    }

    virtual bool remove_voxel_if_exists(const private_coordinate_type& coord)
    {
        return true;
    }

private:
    const Shape::dimension_kind dimension_;

};

} //ecell4

#endif /* __ECELL4_STRUCTURE_TYPE_HPP */
