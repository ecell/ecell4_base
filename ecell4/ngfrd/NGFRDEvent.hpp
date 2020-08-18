#ifndef ECELL4_NGFRD_EVENT_HPP
#define ECELL4_NGFRD_EVENT_HPP

#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/ngfrd/DomainID.hpp>

namespace ecell4
{
namespace ngfrd
{

struct NGFRDEvent
{
    NGFRDEvent(const Real t, const DomainID& did) noexcept
        : time_(t), did_(did)
    {}

    Real      const& time()      const noexcept {return time_;}
    DomainID& const& domain_id() const noexcept {return did_;}

private:
    Real     time_;
    DomainID did_;
};


} // ngfrd
} // ecell4
#endif// ECELL4_NGFRD_EVENT_HPP
