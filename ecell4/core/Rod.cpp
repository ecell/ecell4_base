#include "Rod.hpp"

namespace ecell4
{

Rod::Rod()
    : length_(0.5e-6), radius_(2.0e-6)
{
    ;
}

Rod::Rod(const Real& length, const Real& radius)
    : length_(length), radius_(radius)
{
    //assert(length_ > 0);
    //assert(radius_ > 0)
}

Rod::Rod(const Rod& rhs)
    : length_(rhs.length_), radius_(rhs.radius_)
{
    ;
}

const Real& Rod::lengthX() const
{
    return length_;
}

const Real& Rod::radius() const
{
    return radius_;
}

Real Rod::is_inside(const Real3& pos) const
{
    return distance(pos);
}

Real Rod::distance(const Real3& pos) const
{
    if (pos[0] > length_/2)
        return length(pos - Real3(length_/2, 0, 0)) - radius_;
    if (pos[0] < -length_/2)
        return length(pos - Real3(-length_/2, 0, 0)) - radius_;
    return length(pos - Real3(pos[0], 0, 0)) - radius_;
}

Real3 Rod::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    // The Cylinder Part
    if (rng->uniform(-4*radius_, 3*length_) >= 0)
    {
        const Real x(rng->uniform(-length_/2, length_/2));
        const Real theta(rng->uniform(0, M_PI*2));
        const Real r(sqrtf(rng->uniform(0, radius_*radius_)));
        return Real3(x, r*cos(theta), r*sin(theta));
    }

    // The Terminal Part
    const Real theta(rng->uniform(0, M_PI));
    const Real phi(rng->uniform(0, M_PI));
    const Real r(pow(rng->uniform(0, pow(radius_, 3.0)), 1.0/3.0));
    const Real l(r*sin(phi));
    if (rng->uniform(-1, 1) < 0)
        return Real3(length_/2+l*sin(theta), r*cos(theta), r*cos(phi));
    else
        return Real3(length_/2-l*sin(theta), r*cos(theta), r*cos(phi));
}

RodSurface Rod::surface() const
{
    return RodSurface(length_, radius_);
}

RodSurface::RodSurface()
    : length_(0.5e-6), radius_(2.0e-6)
{
    ;
}

RodSurface::RodSurface(const Real& length, const Real& radius)
    : length_(length), radius_(radius)
{
    //assert(length_ > 0);
    //assert(radius_ > 0)
}

RodSurface::RodSurface(const RodSurface& rhs)
    : length_(rhs.length_), radius_(rhs.radius_)
{
    ;
}

const Real& RodSurface::lengthX() const
{
    return length_;
}

const Real& RodSurface::radius() const
{
    return radius_;
}

Real RodSurface::is_inside(const Real3& pos) const
{
    return distance(pos);
}

Real RodSurface::distance(const Real3& pos) const
{
    if (pos[0] > length_/2)
        return length(pos - Real3(length_/2, 0, 0)) - radius_;
    if (pos[0] < -length_/2)
        return length(pos - Real3(-length_/2, 0, 0)) - radius_;
    return length(pos - Real3(pos[0], 0, 0)) - radius_;
}

Real3 RodSurface::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    // The Cylinder Part
    if (rng->uniform(-2*radius_, length_) >= 0)
    {
        const Real x(rng->uniform(-length_/2, length_/2));
        const Real theta(rng->uniform(0, M_PI*2));
        return Real3(x, radius_*sin(theta), radius_*cos(theta));
    }

    // The Terminal Part
    const Real theta(rng->uniform(0, M_PI));
    const Real phi(rng->uniform(0, M_PI));
    const Real l(radius_*sin(theta));
    if (rng->uniform(-1, 1) < 0)
        return Real3(length_/2+l*sin(theta), l*cos(theta), radius_*cos(phi));
    else
        return Real3(-length_/2-l*sin(theta), l*cos(theta), radius_*cos(phi));
}

Rod RodSurface::inside() const
{
    return Rod(length_, radius_);
}

} // ecell4
