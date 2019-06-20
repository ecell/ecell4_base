#ifndef ECELL4_SPACE_HPP
#define ECELL4_SPACE_HPP

namespace ecell4
{

struct Space
{
    typedef enum
    {
        ELSE,
        PARTICLE,
        LATTICE,
        COMPARTMENT,
        SUBVOLUME
    } space_kind;
};

} // ecell4

#endif /* ECELL4_SPACE_HPP */
