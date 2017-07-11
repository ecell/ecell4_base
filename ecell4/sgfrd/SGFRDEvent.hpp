#ifndef ECELL4_SGFRD_EVENT
#define ECELL4_SGFRD_EVENT
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/sgfrd/Single.hpp>
#include <ecell4/sgfrd/Pair.hpp>
#include <ecell4/sgfrd/Multi.hpp>
#include <boost/variant.hpp>
#include <ostream>
#include <climits>

namespace ecell4
{
namespace sgfrd
{

struct SGFRDEvent
{
public:
    typedef boost::variant<Single, Pair, Multi> domain_type;

    enum domain_kind
    {
        single_domain = 0,
        pair_domain   = 1,
        multi_domain  = 2,
        invalid       = INT_MAX,
        // because std::numeric_limits<int>::max() is not constexpr in c++03,
        // there is no other way but to use macro defined in <climits>.
    };

public:

    template<typename domainT>
    SGFRDEvent(Real const& time, const domainT& dom)
        : time_(time), domain_(dom)
    {}

    Real const&        time()   const {return time_;}
    domain_type const& domain() const {return domain_;}
    domain_type &      domain()       {return domain_;}

    domain_kind which_domain() const
    {
        switch(this->domain_.which())
        {
            case 0: return single_domain;
            case 1: return pair_domain;
            case 2: return multi_domain;
            default: return invalid;
        }
    }

private:

    Real time_;
    domain_type domain_;
};

typedef ecell4::EventSchedulerBase<SGFRDEvent> SGFRDEventScheduler;
typedef SGFRDEventScheduler::identifier_type EventID;
typedef EventID DomainID; // XXX!

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_EVENT
