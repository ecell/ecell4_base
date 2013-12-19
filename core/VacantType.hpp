#ifndef __ECELL4_VACANT_TYPE_HPP
#define __ECELL4_VACANT_TYPE_HPP

#include "MolecularTypeBase.hpp"

namespace ecell4
{

class VacantType
    : public MolecularTypeBase
{

public:

    ~VacantType()
    {
    }

    bool is_vacant() const
    {
        return true;
    }
};

} // ecell4

#endif /* __ECELL4_VACANT_TYPE_HPP */
