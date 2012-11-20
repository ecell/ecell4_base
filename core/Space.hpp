#ifndef __SPACE_HPP
#define __SPACE_HPP

#include <stdexcept>

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

    void set_t(Real const& t)
    {
        if (t < 0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

protected:

    Real t_;
};

} // ecell4

#endif /* __SPACE_HPP */
