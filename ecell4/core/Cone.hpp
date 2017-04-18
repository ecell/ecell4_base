#ifndef ECELL4_CORE_CONE
#define ECELL4_CORE_CONE
#include <ostream>
#include "Shape.hpp"

namespace ecell4
{
struct ConicalSurface;

struct Cone : public Shape
{
public:
    typedef Real3 position_type;
    typedef position_type::value_type length_type;
    typedef position_type::value_type value_type;

public:

    Cone(){}
    ~Cone(){}
    Cone(const Real3 apex)/* infinite cone */
        : radius_(std::numeric_limit<Real>::max()), apex_(apex),
          bottom_(std::numeric_limit<Real>::max(),
                  std::numeric_limit<Real>::max(),
                  std::numeric_limit<Real>::max())
    {}
    Cone(const Cone& rhs)
        : radius_(rhs.radius_), apex_(rhs.apex_), bottom_(rhs.bottom_)
    {}
    Cone& operator=(const Cone& rhs)
    {
        radius_ = rhs.radius_;
        apex_   = rhs.apex_;
        bottom_ = rhs.bottom_;
    }

    Real  const& radius() const {return radius_;}
    Real3 const& apex()   const {return apex_;}
    Real3 const& bottom() const {return bottom_;}

    Real3 const& position() const {return bottom_;} // XXX ?
    Real3&       position()       {return bottom_;}
    Real const&  size() const {return radius_;} // XXX ?
    Real&        size()       {return radius_;}
    dimension_kind dimension() const {return THREE;}

    Real is_inside(const Real3& coord) const
    {
        throw NotImplemented("Cone::is_inside");
    }
    Real distance(const Real3& pos) const
    {
        throw NotImplemented("Cone::distance");
    }
    Real distance_sq(const Real3& pos) const
    {
        throw NotImplemented("Cone::distance_sq");
    }

    ConicalSurface surface() const {return ConicalSurface(center_, radius_, normal);}

    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
    {
        throw NotImplemented("Cone::draw_position");
    }

    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("Cone::test_AABB");
    }

protected:

    Real radius_;// bottom radius
    Real3 apex_;
    Real3 bottom_;

};

struct ConicalSurface : public Shape
{
public:
    typedef Real3 position_type;
    typedef position_type::value_type length_type;
    typedef position_type::value_type value_type;

public:

    ConicalSurface(){}
    ~ConicalSurface(){}
    ConicalSurface(const Real3 apex)
        : radius_(std::numeric_limit<Real>::max()), apex_(apex)
    {}
    ConicalSurface(const ConicalSurface& rhs)
        : radius_(rhs.radius_), apex_(rhs.apex_), bottom_(rhs.bottom_)
    {}
    ConicalSurface& operator=(const ConicalSurface& rhs)
    {
        radius_ = rhs.radius_;
        apex_   = rhs.apex_;
        bottom_ = rhs.bottom_;
    }

    Real  const& radius() const {return radius_;}
    Real3 const& apex()   const {return apex_;}
    Real3 const& bottom() const {return bottom_;}

    Real3 const& position() const {return bottom_;} // XXX ?
    Real3&       position()       {return bottom_;}
    Real const&  size() const {return radius_;} // XXX ?
    Real&        size()       {return radius_;}
    dimension_kind dimension() const {return THREE;}

    Real is_inside(const Real3& coord) const
    {
        return distance(coord);
    }
    Real distance(const Real3& pos) const
    {
        const std::pair<bool, Real> result = this->distance_sq_(pos);
        return result.first ? result.second : -(result.second);
    }
    Real distance_sq(const Real3& pos) const
    {
        return this->distance_sq_(pos).second;
    }

    ConicalSurface surface() const {return ConicalSurface(center_, radius_, normal);}

    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
    {
        throw NotImplemented("ConicalSurface::draw_position");
    }

    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("ConicalSurface::test_AABB");
    }

protected:
    // pairof(is_outside, dist)
    std::pair<bool, Real> distance_sq_(const Real3& pos) const
    {
        throw NotImplemented("protected: ConicalSurface::distance_sq_");
    }

protected:

    Real radius_;// bottom radius
    Real3 apex_;
    Real3 bottom_;

};

template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Cone& c)
{
    os << "Cone(apex=" << c.apex() << ", bottom=" << c.bottom()
       << ", radius=" << c.radius() << ")";
    return os;
}


template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const ConicalSurface& c)
{
    os << "ConicalSurface(apex=" << c.apex() << ", bottom=" << c.bottom()
       << ", radius=" << c.radius() << ")";
    return os;
}


} // ecell4
#endif //ECELL4_CORE_CONE
