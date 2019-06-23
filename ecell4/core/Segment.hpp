#ifndef ECELL4_CORE_SEGMENT
#define ECELL4_CORE_SEGMENT
#include <ecell4/core/Shape.hpp>
#include <ecell4/core/exceptions.hpp>

namespace ecell4
{

struct Segment
    : public Shape
{
public:

    /** for epdp
     */
    typedef Real3 position_type;
    typedef position_type::value_type length_type;
    typedef position_type::value_type value_type;

public:

    Segment(){}
    Segment(const Real3& start, const Real3& stop): start_(start), stop_(stop){}
    Segment(const Segment& rhs) : start_(rhs.start_), stop_(rhs.stop_){}

    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
    {
        const Real r = rng->random();
        return start_ * r + stop_ * (1. - r);
    }
    Real is_inside(const Real3& coord) const
    {
        throw NotImplemented("Segment::is_inside");
    }
    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("Segment::is_inside");
    }

    Real3 const& start() const {return start_;}
    Real3&       start()       {return start_;}
    Real3 const& stop()  const {return stop_;}
    Real3&       stop()        {return stop_;}

    dimension_kind dimension() const
    {
        return ONE;
    }

protected:

    Real3 start_;
    Real3 stop_;
};

template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Segment& sgm)
{
    os << "Segment(" << sgm.start() << ", " << sgm.stop() << ')';
    return os;
}

} // ecell4
#endif// ECELL4_CORE_SEGMENT
