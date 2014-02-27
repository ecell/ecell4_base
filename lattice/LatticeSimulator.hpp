#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Simulator.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>

#include "LatticeWorld.hpp"
#include "EventScheduler.hpp"

namespace ecell4
{

namespace lattice
{

class LatticeSimulator
    : public Simulator
{
protected:
    struct Event : EventScheduler::Event
    {
        Event(LatticeSimulator* sim, const Species& species)
            : EventScheduler::Event(0.0), species_(species), sim_(sim)
        {
            const Real R(sim_->world_->normalized_voxel_radius());
            Real D = boost::lexical_cast<Real>(species.get_attribute("D"));
            if (D <= 0)
            {
                dt_ = 0;
            } else {
                dt_ = 4 * R / 6 / D;
            }

            time_ = dt_;
        }

        virtual ~Event()
        {
        }

        virtual void fire()
        {
            boost::shared_ptr<GSLRandomNumberGenerator> rng(sim_->world_->rng());

            std::vector<Coord> coords(sim_->world_->list_coords(species_));
            shuffle(*rng, coords);
            for (std::vector<Coord>::iterator itr(coords.begin());
                    itr != coords.end(); ++itr)
            {
                const Coord coord(*itr);
                const Integer nrnd(rng->uniform_int(0,11));
                std::pair<Coord, bool> retval(sim_->world_->move_to_neighbor(coord, nrnd));
                if (!retval.second)
                {
                    sim_->attempt_reaction_(coord, retval.first);
                }
            }

            time_ += dt_;
        }

    protected:
        Species species_;
        LatticeSimulator* sim_;
        MolecularTypeBase* mt_;
        Real dt_;
    };


public:

    LatticeSimulator(
            boost::shared_ptr<NetworkModel> model,
            boost::shared_ptr<LatticeWorld> world)
        : model_(model), world_(world), is_initialized_(false)
    {
    }

    virtual Real t() const
    {
        return world_->t();
    }

    virtual Real dt() const
    {
        return dt_;
    }

    Integer num_steps() const
    {
        return (Integer)(t() / dt());
    }

    void initialize();
    void step();
    bool step(const Real& upto);

protected:
    boost::shared_ptr<EventScheduler::Event> create_event(const Species& species);
    void attempt_reaction_(Coord from_coord, Coord to_coord);

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<LatticeWorld> world_;
    EventScheduler scheduler_;

    Real dt_;
    bool is_initialized_;

};

} // lattice

} // ecell4

#endif /* __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP */
