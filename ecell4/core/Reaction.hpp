#ifndef __ECELL4_REACTION_HPP
#define __ECELL4_REACTION_HPP

#include <vector>

#include "ReactionRule.hpp"
#include "Identifier.hpp"

namespace ecell4
{

template <typename T>
struct Reaction
{
    typedef std::pair<ParticleID, T> particle_type;

    std::vector<particle_type> reactants;
    std::vector<particle_type> products;
    ReactionRule rule;
};

} // ecell4

#endif /* __ECELL4_REACTION_HPP */
