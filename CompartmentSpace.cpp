#include "CompartmentSpace.hpp"


namespace ecell4
{

Real const& CompartmentSpaceVectorImpl::volume() const
{
    return volume_;
}

void CompartmentSpaceVectorImpl::set_volume(Real volume)
{
    if (volume >= 0)
    {
        volume_ = volume;
    }
}

} // ecell4
