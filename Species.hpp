#ifndef __SPECIES_HPP
#define __SPECIES_HPP

#include <vector>

#include "types.hpp"


namespace ecell4
{

class Species
{
public:

    bool operator==(Species const& rhs) const;
    bool operator<(Species const& rhs) const;
    bool operator>(Species const& rhs) const;
};

typedef std::vector<Species> SpeciesVector;

} // ecell4

#endif /* __SPECIES_HPP */
