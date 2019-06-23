#ifndef ECELL4_SGFRD_BIRTH_DOMAIN
#define ECELL4_SGFRD_BIRTH_DOMAIN
#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>

namespace ecell4
{
namespace sgfrd
{

class Birth
{
  public:

    typedef Species      species_type;
    typedef ReactionRule reaction_rule_type;

  public:

    Birth(const Real dt, const Real begin_time,
          const Species& sp, const ReactionRule& rule)
        : dt_(dt), begin_time_(begin_time), species_(sp), rule_(rule)
    {}
    ~Birth(){}

    Real& dt()       throw() {return dt_;}
    Real  dt() const throw() {return dt_;}
    Real& begin_time()       throw() {return begin_time_;}
    Real  begin_time() const throw() {return begin_time_;}

    Species&       species()       throw() {return species_;}
    Species const& species() const throw() {return species_;}
    ReactionRule&       rule()       throw() {return rule_;}
    ReactionRule const& rule() const throw() {return rule_;}

    std::size_t num_shells() const {return 1;}
    std::size_t multiplicity() const {return 1;}

  private:

    Real dt_;
    Real begin_time_;
    Species      species_;
    ReactionRule rule_;
};


} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SINGLE_DOMAIN */
