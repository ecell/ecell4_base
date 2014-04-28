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
    : public Simulator<NetworkModel, GillespieWorld>
{
public:

    typedef Simulator<NetworkModel, GillespieWorld> base_type;

public:

    GillespieSimulator(
        boost::shared_ptr<NetworkModel> model,
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

    /**
     * recalculate reaction propensities and draw the next time.
     */
    void initialize();

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    void draw_next_reaction(void);

protected:

    Real dt_;
    int next_reaction_num_; // the index of the next reaction.
};

}

} // ecell4

#endif /* __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP */
