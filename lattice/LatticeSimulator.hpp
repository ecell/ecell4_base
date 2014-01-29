#ifndef __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP
#define __ECELL4_LATTICE_LATTICE_SIMULATOR_HPP

#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <ecell4/core/NetworkModel.hpp>
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
            : EventScheduler::Event(0.0), sim_(sim)
        {
            mt_ = sim_->world_->get_molecular_type(species);
            const Real R(sim_->world_->normalized_voxel_radius());
            Real D = boost::lexical_cast<Real>(species.get_attribute("D"));
            if (D <= 0)
            {
                dt_ = 0;
            } else {
                dt_ = 4 * R / 6 / D;
            }
        }

        virtual ~Event()
        {
        }

        virtual void fire()
        {
            std::cerr << "[" << time_ << "," << dt_ << "," <<
                mt_->species().name() << "," << mt_->species().get_attribute("D") << "]";
            boost::shared_ptr<GSLRandomNumberGenerator> rng(sim_->world_->rng());

            shuffle(*rng, mt_->voxels());
            for (MolecularType::container_type::iterator itr(mt_->begin());
                    itr != mt_->end(); ++itr)
            {
                Coord from_coord((*itr).first);
                const Integer rnd(rng->uniform_int(0,11));
                Coord to_coord(sim_->world_->get_neighbor(from_coord,
                            rnd));

                sim_->world_->move(from_coord, to_coord);
            }

            time_ += dt_;
        }

    protected:
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
        return (*world_).t();
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
