#ifndef __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP
#define __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "GillespieWorld.hpp"


namespace ecell4
{

namespace gillespie
{

class GillespieSimulator
    : public Simulator<Model, GillespieWorld>
{
public:

    typedef Simulator<Model, GillespieWorld> base_type;

public:

    GillespieSimulator(
        boost::shared_ptr<Model> model,
        boost::shared_ptr<GillespieWorld> world)
        : base_type(model, world)
    {
        initialize();
    }

    GillespieSimulator(boost::shared_ptr<GillespieWorld> world)
        : base_type(world)
    {
        initialize();
    }

    // SimulatorTraits

    Real t(void) const;
    Real dt(void) const;

    void step(void) ;
    bool step(const Real & upto);

    // Optional members

    void set_t(const Real &t);
    std::vector<ReactionRule> last_reactions() const;

    /**
     * recalculate reaction propensities and draw the next time.
     */
    void initialize();

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    bool __draw_next_reaction(void);
    void draw_next_reaction(void);
    Integer num_molecules(const Species& sp);
    Integer num_molecules(const Species& sp1, const Species& sp2);
    ReactionRule draw_exact_reaction(const ReactionRule& rr);
    std::pair<ReactionRule::reactant_container_type, Integer>
        draw_exact_reactants(const Species& sp1);
    std::pair<ReactionRule::reactant_container_type, Integer>
        draw_exact_reactants(const Species& sp1, const Species& sp2);

    Real calculate_propensity(const ReactionRule& rr);

protected:

    Real dt_;
    ReactionRule next_reaction_;
    std::vector<ReactionRule> last_reactions_;
};

}

} // ecell4

#endif /* __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP */
