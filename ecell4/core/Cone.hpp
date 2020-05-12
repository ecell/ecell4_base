#ifndef ECELL4_CORE_CONE
#define ECELL4_CORE_CONE
#include <ostream>
#include <limits>
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
    Cone(const Real3 apex, const Real apex_angl)/* infinite cone */
        : slant_height_(std::numeric_limits<Real>::max()),
          apex_angle_(apex_angl), apex_(apex)
    {}
    Cone(const Real3 apex, const Real apex_angl, const Real slant)
        : slant_height_(slant), apex_angle_(apex_angl), apex_(apex)
    {}
    Cone(const Cone& rhs)
        : slant_height_(rhs.slant_height_), apex_angle_(rhs.apex_angle_),
          apex_(rhs.apex_)
    {}
    Cone& operator=(const Cone& rhs)
    {
        slant_height_ = rhs.slant_height_;
        apex_angle_   = rhs.apex_angle_;
        apex_         = rhs.apex_;
        return *this;
    }

    Real  const& slant_height() const {return slant_height_;}
    Real  const& apex_angle()   const {return apex_angle_;}
    Real3 const& apex()         const {return apex_;}

    Real3 const& position() const {return apex_;}
    Real3&       position()       {return apex_;}
    Real const&  size() const {return slant_height_;}
    Real&        size()       {return slant_height_;}
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

    ConicalSurface surface() const;

    Real3 draw_position(std::shared_ptr<RandomNumberGenerator>& rng) const
    {
        throw NotImplemented("Cone::draw_position");
    }

    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("Cone::test_AABB");
    }

protected:

    Real slant_height_;//
    Real apex_angle_;
    Real3 apex_;

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
    ConicalSurface(const Real3 apex, const Real apex_angl)
        : slant_height_(std::numeric_limits<Real>::max()),
          apex_angle_(apex_angl), apex_(apex)
    {}
    ConicalSurface(const Real3 apex, const Real apex_angl, const Real slant)
        : slant_height_(slant), apex_angle_(apex_angl), apex_(apex)
    {}

    ConicalSurface(const ConicalSurface& rhs)
        : slant_height_(rhs.slant_height_), apex_angle_(rhs.apex_angle_),
          apex_(rhs.apex_)
    {}
    ConicalSurface& operator=(const ConicalSurface& rhs)
    {
        slant_height_     = rhs.slant_height_;
        apex_angle_ = rhs.apex_angle_;
        apex_       = rhs.apex_;
        return *this;
    }

    Real  const& slant_height()     const {return slant_height_;}
    Real  const& apex_angle() const {return apex_angle_;}
    Real3 const& apex()   const {return apex_;}

    Real3 const& position() const {return apex_;}
    Real3&       position()       {return apex_;}
    Real const&  size() const {return slant_height_;}
    Real&        size()       {return slant_height_;}
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

    Cone inside() const {return Cone(apex_, apex_angle_, slant_height_);}

    Real3 draw_position(std::shared_ptr<RandomNumberGenerator>& rng) const
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

    Real slant_height_;
    Real apex_angle_;
    Real3 apex_;
};

inline ConicalSurface Cone::surface() const
{
    return ConicalSurface(apex_, apex_angle_, slant_height_);
}

template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Cone& c)
{
    os << "Cone(apex=" << c.apex()  << ", apex_angle = " << c.apex_angle()
       << ", slant_height=" << c.slant_height() << ")";
    return os;
}


template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const ConicalSurface& c)
{
    os << "ConicalSurface(apex=" << c.apex() << ", apex_angle = " << c.apex_angle()
       << ", slant_height=" << c.slant_height() << ")";
    return os;
}


} // ecell4
#endif //ECELL4_CORE_CONE
