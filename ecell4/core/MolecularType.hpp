#ifndef ECELL4_MOLECULAR_TYPE_HPP
#define ECELL4_MOLECULAR_TYPE_HPP

#include "VoxelPool.hpp"
#include "VacantType.hpp"

namespace ecell4
{

class MolecularType
    : public MoleculePool
{
public:

    typedef MoleculePool base_type;
    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::coordinate_type coordinate_type;
    typedef base_type::container_type container_type;
    typedef base_type::iterator iterator;
    typedef base_type::const_iterator const_iterator;
    typedef base_type::voxel_type_type voxel_type_type;

public:

    MolecularType(const std::string& name = "")
        : base_type(Species(name), &(VacantType::getInstance()), 0, 0)
    {
        ;
    }

    MolecularType(const Species& species, const Real& radius = 0.0, const Real& D = 0.0)
        : base_type(species, &(VacantType::getInstance()), radius, D)
    {
        ;
    }

    MolecularType(const Species& species, VoxelPool* location,
                  const Real& radius = 0.0, const Real& D = 0.0)
        : base_type(species, location, radius, D)
    {
        ;
    }

    ~MolecularType()
    {
        ;
    }

    voxel_type_type const voxel_type() const
    {
        return DEFAULT;
    }

    virtual const Shape::dimension_kind get_dimension() const
    {
        return location()->get_dimension();
    }
};

} // ecell4

#endif /* ECELL4_MOLECULAR_TYPE_HPP */
