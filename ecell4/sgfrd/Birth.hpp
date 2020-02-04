#ifndef ECELL4_SGFRD_BIRTH_DOMAIN
#define ECELL4_SGFRD_BIRTH_DOMAIN
#include <ecell4/core/ReactionRule.hpp>

namespace ecell4
{
namespace sgfrd
{

class Birth
{
  public:

    typedef ReactionRule reaction_rule_type;

  public:

    Birth(const Real dt, const Real begin_time, const ReactionRule& rule)
        : dt_(dt), begin_time_(begin_time), rule_(rule)
    {}
    ~Birth(){}

    Real& dt()       noexcept {return dt_;}
    Real  dt() const noexcept {return dt_;}
    Real& begin_time()       noexcept {return begin_time_;}
    Real  begin_time() const noexcept {return begin_time_;}

    ReactionRule&       rule()       noexcept {return rule_;}
    ReactionRule const& rule() const noexcept {return rule_;}

    std::size_t num_shells()   const noexcept {return 1;}
    std::size_t multiplicity() const noexcept {return 1;}

  private:

    Real dt_;
    Real begin_time_;
    ReactionRule rule_;
};
} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SINGLE_DOMAIN */
