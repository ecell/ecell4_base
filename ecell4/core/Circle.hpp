#ifndef ECELL4_CIRCLE_HPP
#define ECELL4_CIRCLE_HPP
#include "Shape.hpp"
#include <ostream>
#include <cmath>

namespace ecell4
{

struct Circumference;

struct Circle : public Shape
{
public:

    Circle(){}
    ~Circle(){}

    Circle(const Real radius, const Real3& center, const Real3 normal)
        : radius_(radius), center_(center), normal_(normal)
    {}
    Circle(const Circle& rhs)
        : radius_(rhs.radius_), center_(rhs.center_), normal_(rhs.normal_)
    {}
    Circle& operator=(const Circle& rhs)
    {
        radius_ = rhs.radius_;
        center_ = rhs.center_;
        normal_ = rhs.normal_;
        return *this;
    }

    Real  const& radius() const {return radius_;}
    Real3 const& center() const {return center_;}
    Real3 const& normal() const {return normal_;}

    Real3 const& position() const {return center_;}
    Real3&       position()       {return center_;}
    Real const&  size() const {return radius_;}
    Real&        size()       {return radius_;}
    dimension_kind dimension() const {return TWO;}

    Real is_inside(const Real3& coord) const
    {
        return this->distance(coord);
    }
    Real distance(const Real3& pos) const
    {
        const std::pair<bool, Real> result = this->distance_sq_(pos);
        return result.first ? std::sqrt(result.second) : -std::sqrt(result.second);
    }
    Real distance_sq(const Real3& pos) const
    {
        return this->distance_sq_(pos).second;
    }

    Circumference surface() const;

    Real3 draw_position(std::shared_ptr<RandomNumberGenerator>& rng) const
    {
        throw NotImplemented("Circle::draw_position");
    }

    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("Circle::test_AABB");
    }

protected:

    std::pair<bool, Real> distance_sq_(const Real3& pos) const
    {
        const Real dot = dot_product(center_ - pos, normal_);
        const bool sign = (dot < 0); // true means + side
        const Real3 projection(normal_ * dot);
        const Real dist_on_plane2 = length_sq(pos + projection - center_);
        if(dist_on_plane2 > radius_ * radius_)
        {
            const Real dr = std::sqrt(dist_on_plane2) - radius_;
            return std::make_pair(sign, length_sq(projection) + dr * dr);
        }
        else
        {
            return std::make_pair(sign, length_sq(projection));
        }
    }

protected:

    Real radius_;
    Real3 center_;
    Real3 normal_;
};

struct Circumference : public Shape
{
    Circumference(){}
    ~Circumference(){}

    Circumference(const Real radius, const Real3& center, const Real3& normal)
        : radius_(radius), center_(center), normal_(normal)
    {}
    Circumference(const Circumference& rhs)
        : radius_(rhs.radius()), center_(rhs.center()), normal_(rhs.normal())
    {}

    Real  const& radius() const {return radius_;}
    Real3 const& center() const {return center_;}
    Real3 const& normal() const {return normal_;}

    Real is_inside(const Real3& coord) const
    {
        return this->distance(coord);
    }
    Real distance(const Real3& pos) const
    {
        const std::pair<bool, Real> result = this->distance_sq_(pos);
        return result.first ? std::sqrt(result.second) : -std::sqrt(result.second);
    }
    Real distance_sq(const Real3& pos) const
    {
        return this->distance_sq_(pos).second;
    }

    Circle inside() const {return Circle(radius_, center_, normal_);}

    Real3 draw_position(std::shared_ptr<RandomNumberGenerator>& rng) const
    {
        throw NotImplemented("Circumference::draw_position");
    }
    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("Circumference::test_AABB");
    }

    dimension_kind dimension() const
    {
        return ONE;
    }

protected:

    std::pair<bool, Real> distance_sq_(const Real3& pos) const
    {
        const Real dot = dot_product(center_ - pos, normal_);
        const bool sign = (dot < 0); // true means + side
        const Real3 projection(normal_ * dot);
        const Real dist_on_plane2 = length_sq(pos + projection - center_);
        if(dist_on_plane2 > radius_ * radius_)
        {
            const Real dr = std::sqrt(dist_on_plane2) - radius_;
            return std::make_pair(sign, length_sq(projection) + dr * dr);
        }
        else
        {
            return std::make_pair(sign, length_sq(projection));
        }
    }

protected:

    Real radius_;
    Real3 center_;
    Real3 normal_;
};

inline Circumference Circle::surface() const
{
    return Circumference(radius_, center_, normal_);
}


template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Circle& c)
{
    os << "Circle(r=" << c.radius() << ", o=" << c.center()
       << ", n=" << c.normal() << ")";
    return os;
}


template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Circumference& c)
{
    os << "Circumference(r=" << c.radius() << ", o=" << c.center()
       << ", n=" << c.normal() << ")";
    return os;
}


} // ecell4
#endif //ECELL4_CIRCLE_HPP
