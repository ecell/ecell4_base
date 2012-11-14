#ifndef __COMPARTMENT_SPACE_HPP
#define __COMPARTMENT_SPACE_HPP

#include "types.hpp"
#include "Space.hpp"


namespace ecell4
{

class CompartmentSpace
    : public Space
{
public:

    virtual Real const& volume() const = 0;
    virtual void set_volume(Real volume) = 0;
};

class CompartmentSpaceVectorImpl
    : public CompartmentSpace
{
public:

    CompartmentSpaceVectorImpl(Real const& volume)
        : volume_(1)
    {
        set_volume(volume);
    }

    Real const& volume() const;
    void set_volume(Real volume);

protected:

    Real volume_;
};

} // ecell4

#endif /* __COMPARTMENT_SPACE_HPP */
