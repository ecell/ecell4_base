#ifndef __ECELL4_PLANAR_SURFACE_HPP
#define __ECELL4_PLANAR_SURFACE_HPP

#include "Shape.hpp"

#include "exceptions.hpp"

namespace ecell4
{

struct PlanarSurface
    : public Shape
{

    PlanarSurface();
    PlanarSurface(const Real3& origin, const Real3& e0, const Real3& e1);
    PlanarSurface(const PlanarSurface& rhs);
    Real is_inside(const Real3& coord) const;
    Real3 draw_position(
        boost::shared_ptr<RandomNumberGenerator>& rng) const;
    bool test_AABB(const Real3& lower, const Real3& upper) const;
    void bounding_box(
        const Real3& edge_lengths, Real3& lower, Real3& u) const;

    const Real3& origin() const
    {
        return origin_;
    }

    const Real3& e0() const
    {
        return e0_;
    }

    const Real3& e1() const
    {
        return e1_;
    }

    const Real3& normal() const
    {
        return n_;
    }

    dimension_kind dimension() const
    {
        return TWO;
    }

    Real3 reflection(const Real3& from, const Real3& displacement) const
    {
        Real3 provisional(from + displacement);
        Real3 new_pos;
        Real is_inside_from(this->is_inside(from));
        Real is_inside_provisional(this->is_inside(provisional));

        if ( (0 < is_inside_from && 0 < is_inside_provisional) || 
                (is_inside_from < 0 && is_inside_provisional < 0) )
        {
            // The same signs means No refrection
            new_pos = provisional;
        }
        else if( 0 < is_inside_from  )
        {
            // inside -> (Refrection) -> inside
            ecell4::Real distance_from_surface(std::abs(is_inside_provisional));
            new_pos = provisional + multiply(this->normal(), (-2.0) * distance_from_surface);
        } else if ( is_inside_from < 0 )
        {
            // outside -> (Refrection) -> outside
            ecell4::Real distance_from_surface(std::abs(is_inside_provisional));
            new_pos = provisional + multiply(this->normal(), 2.0 * distance_from_surface);
        }
        return new_pos;
    }

    bool cross(const Real3 &from, const Real3& displacement) const
    {
        Real is_inside = this->is_inside(from);
        Real is_inside_dest = this->is_inside(from + displacement);
        if (0. < (is_inside * is_inside_dest))
        {
            // Same signs means the same side between from and after
            return false;
        }
        else
        {
            return true;
        }
    }
    Real3 intrusion(const Real3& from, const Real3& displacement) const
    {
        Real is_inside = this->is_inside(from);
        Real is_inside_dest = this->is_inside(from + displacement);
        if (0. < (is_inside * is_inside_dest))
        {
            throw IllegalState("No Intrusion");
        }
        else
        {
            Real ratio(std::abs(is_inside) / (std::abs(is_inside) + std::abs(is_inside_dest)));
            return from + multiply(displacement, ratio);
        }
    }

protected:

    Real3 origin_, e0_, e1_, n_;
    Real d_;
};

}

#endif /* __ECELL4_PLANAR_SURFACE_HPP */
