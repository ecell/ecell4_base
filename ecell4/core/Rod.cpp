#include "Rod.hpp"
#include "AABB.hpp"
#include "exceptions.hpp"
#include "collision.hpp"


namespace ecell4
{

Rod::Rod()
    : length_(0.5e-6), radius_(2.0e-6), origin_()
{
    ;
}

Rod::Rod(const Real& length, const Real& radius)
    : length_(length), radius_(radius), origin_()
{
    //assert(length_ > 0);
    //assert(radius_ > 0)
}

Rod::Rod(const Real& length, const Real& radius, const Real3& origin)
    : length_(length), radius_(radius), origin_(origin)
{
    ;
}

Rod::Rod(const Rod& rhs)
    : length_(rhs.length_), radius_(rhs.radius_), origin_(rhs.origin_)
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

const Real3& Rod::origin() const
{
    return origin_;
}

void Rod::shift(const Real3& vec)
{
    origin_ += vec;
}

Real Rod::is_inside(const Real3& pos) const
{
    return distance(pos);
}

Real Rod::distance(const Real3& pos) const
{
    return collision::distance_point_capsule(pos, *this);
}

Real3 Rod::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    // The Cylinder Part
    if (rng->uniform(-4*radius_, 3*length_) >= 0)
    {
        const Real x(rng->uniform(-length_/2, length_/2));
        const Real theta(rng->uniform(0, M_PI*2));
        const Real r(sqrt(rng->uniform(0, pow(radius_, 2.0))));
        return origin_ + Real3(x, r*cos(theta), r*sin(theta));
    }

    // The Terminal Part
    const Real theta(rng->uniform(0, M_PI));
    const Real phi(rng->uniform(0, M_PI));
    const Real r(pow(rng->uniform(0, pow(radius_, 3.0)), 1.0/3.0));
    const Real l(r*sin(phi));

    const Integer sign(2*Integer(rng->uniform(0,2))-1);
    return origin_ + Real3(sign*(length_/2+l*sin(theta)), l*cos(theta), r*cos(phi));
}

RodSurface Rod::surface() const
{
    return RodSurface(length_, radius_, origin_);
}

bool Rod::test_AABB(const Real3& lower, const Real3& upper) const
{
    const Real3 axis(1.0, 0.0, 0.0); //XXX: DEFAULT
    const Real3 d(axis * length_);
    const Real3 p0(origin_ - axis * (length_ * 0.5));

    Real t;
    return collision::intersect_moving_sphere_AABB(
        Sphere(p0, radius_), d, AABB(lower, upper), t);
}

RodSurface::RodSurface()
    : length_(0.5e-6), radius_(2.0e-6), origin_()
{
    ;
}

RodSurface::RodSurface(const Real& length, const Real& radius)
    : length_(length), radius_(radius), origin_()
{
    //assert(length_ > 0);
    //assert(radius_ > 0)
}

RodSurface::RodSurface(const Real& length, const Real& radius, const Real3& origin)
    : length_(length), radius_(radius), origin_(origin)
{
    ;
}

RodSurface::RodSurface(const RodSurface& rhs)
    : length_(rhs.length_), radius_(rhs.radius_), origin_(rhs.origin_)
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

const Real3& RodSurface::origin() const
{
    return origin_;
}

void RodSurface::shift(const Real3& vec)
{
    origin_ += vec;
}

Real RodSurface::is_inside(const Real3& pos) const
{
    return distance(pos);
}

Real RodSurface::distance(const Real3& pos) const
{
    return collision::distance_point_capsule(pos, this->inside());
}

Real3 RodSurface::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    // The Cylinder Part
    if (rng->uniform(-2*radius_, length_) >= 0)
    {
        const Real x(rng->uniform(-length_/2, length_/2));
        const Real theta(rng->uniform(0, M_PI*2));
        return origin_ + Real3(x, radius_*sin(theta), radius_*cos(theta));
    }

    // The Terminal Part
    const Real theta(rng->uniform(0, M_PI));
    const Real phi(rng->uniform(0, M_PI));
    const Real l(radius_*sin(phi));

    const Integer sign(2*Integer(rng->uniform(0,2))-1);
    return origin_ + Real3(sign*(length_/2+l*sin(theta)), l*cos(theta), radius_*cos(phi));
}

Rod RodSurface::inside() const
{
    return Rod(length_, radius_, origin_);
}

bool RodSurface::test_AABB(const Real3& lower, const Real3& upper) const
{
    throw NotImplemented("not implemented yet.");
}

} // ecell4
