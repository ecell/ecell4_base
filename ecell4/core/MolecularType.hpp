#ifndef __ECELL4_MOLECULAR_TYPE_HPP
#define __ECELL4_MOLECULAR_TYPE_HPP

#include "MolecularTypeBase.hpp"
#include "VacantType.hpp"

namespace ecell4
{

class MolecularType
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

    MolecularType(const std::string& name = "")
        : base_type(Species(name), false, &(VacantType::getInstance()), 0, 0)
    {
        ;
    }

    MolecularType(const Species& species, const Real& radius = 0.0,
            const Real& D = 0.0)
        : base_type(species, false, &(VacantType::getInstance()), radius, D)
    {
        ;
    }

    MolecularType(const Species& species, MolecularTypeBase* location,
            const Real& radius = 0.0, const Real& D = 0.0)
        : base_type(species, false,  location, radius, D)
    {
        ;
    }

    MolecularType(const Species& species, bool is_structure,  MolecularTypeBase* location,
            const Real& radius = 0.0, const Real& D = 0.0)
        : base_type(species, is_structure,  location, radius, D)
    {
        ;
    }


    ~MolecularType()
    {
        ;
    }

    bool is_vacant() const
    {
        return false;
    }

    virtual const Shape::dimension_kind get_dimension() const
    {
        return location()->get_dimension();
    }

};

} // ecell4

#endif /* __ECELL4_MOLECULAR_TYPE_HPP */
