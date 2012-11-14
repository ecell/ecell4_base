#ifndef __SPACE_HPP
#define __SPACE_HPP

#include "types.hpp"


namespace ecell4
{

class Space
{
public:

    Space()
        : t_(0)
    {
        ;
    }

    Real const& t() const
    {
        return t_;
    }

protected:

    Real t_;
};

} // ecell4

#endif /* __SPACE_HPP */
