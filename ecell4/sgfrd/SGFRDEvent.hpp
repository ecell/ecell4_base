#ifndef ECELL4_SGFRD_EVENT
#define ECELL4_SGFRD_EVENT
#include <ecell4/core/EventScheduler.hpp>
#include "DomainID.hpp"

namespace ecell4
{
namespace sgfrd
{

struct SGFRDEvent
{
public:

    SGFRDEvent(Real const& time, const DomainID& did)
        : time_(time), did_(did)
    {}

    Real     time()      const {return time_;}
    DomainID domain_id() const {return did_;}

private:

    Real time_;
    DomainID did_;
};

typedef ecell4::EventSchedulerBase<SGFRDEvent> SGFRDEventScheduler;
typedef SGFRDEventScheduler::identifier EventID;

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_EVENT
