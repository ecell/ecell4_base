#ifndef __GILLESPIESIMULATOR_HPP
#define __GILLESPIESIMULATOR_HPP

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
    : public Simulator
{
public:

    GillespieSimulator(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<GillespieWorld> world)
        : model_(model), world_(world), num_steps_(0)
    {
        this->initialize();
    }

    // SimulatorTraits

    Real t(void) const;
    Real dt(void) const;
    Integer num_steps(void) const;

    void step(void) ;
    bool step(Real const & upto);

    // Optional members

    void set_t(Real const &t);

    /**
     * recalculate reaction propensities and draw the next time.
     */
    void initialize(void);

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    void draw_next_reaction(void);

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<GillespieWorld> world_;
    Integer num_steps_;

    Real dt_;
    int next_reaction_num_; // the index of the next reaction.
};

}

} // ecell4

#endif //__GILLESPIESIMULATOR_HPP
