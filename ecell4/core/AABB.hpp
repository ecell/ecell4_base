#ifndef __ECELL4_AABB_HPP
#define __ECELL4_AABB_HPP

#include "Shape.hpp"

namespace ecell4
{

struct AABB
    : public Shape
{
    AABB()
        : lower_(), upper_()
    {
        ;
    }

    AABB(const Real3& lower, const Real3& upper)
        : lower_(lower), upper_(upper)
    {
        ;
    }

    AABB(const AABB& rhs)
        : lower_(rhs.lower()), upper_(rhs.upper())
    {
        ;
    }

    const Real3& lower() const
    {
        return lower_;
    }

    const Real3& upper() const
    {
        return upper_;
    }

    Real distance(const Real3& pos) const
    {
        const Real3 closest(
            std::min(std::max(pos[0], lower_[0]), upper_[0]),
            std::min(std::max(pos[1], lower_[1]), upper_[1]),
            std::min(std::max(pos[2], lower_[2]), upper_[2]));
        const Real distance_plus(length(closest - pos));

        if (distance_plus > 0)
        {
            return distance_plus;
        }
        else
        {
            std::vector<Real> tmp(6);
            tmp[0] = lower_[0] - pos[0];
            tmp[1] = lower_[1] - pos[1];
            tmp[2] = lower_[2] - pos[2];
            tmp[3] = pos[0] - upper_[0];
            tmp[4] = pos[1] - upper_[1];
            tmp[5] = pos[2] - upper_[2];
            return *std::max_element(tmp.begin(), tmp.end());
        }
    }

    Real is_inside(const Real3& coord) const
    {
        return distance(coord);
    }

    Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const
    {
        const Real3 pos(
            rng->uniform(lower_[0], upper_[0]),
            rng->uniform(lower_[1], upper_[1]),
            rng->uniform(lower_[2], upper_[2]));
        return pos;
    }

protected:

    Real3 lower_, upper_;
};

}// ecell4

#endif /* __ECELL4_AABB_HPP */
