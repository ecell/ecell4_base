#ifndef __ECELL4_CONTEXT_HPP
#define __ECELL4_CONTEXT_HPP

#include "Species.hpp"


namespace ecell4
{

class Context
{
public:

    Context()
    {
        ;
    }

    virtual ~Context()
    {
        ;
    }

protected:
};

bool uspmatch(const UnitSpecies& pttrn, const Species& sp);
bool spmatch(const Species& pttrn, const Species& sp);

} // ecell4

#endif /* __ECELL4_CONTEXT_HPP */
