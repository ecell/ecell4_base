#include "LatticeSimulator.hpp"


namespace ecell4
{

namespace lattice
{

void LatticeSimulator::initialize()
{
    //(*world_).initialize();
    dt_ = 0.1;
    // dt = (4 * R^2) / (6 * D) for each species

    std::vector<Species> species(world_->list_species());
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        const boost::shared_ptr<EventScheduler::Event> event(create_event(*itr));
        scheduler_.add(event);
    }

    is_initialized_ = true;
}

boost::shared_ptr<EventScheduler::Event> LatticeSimulator::create_event(
        const Species& species)
{
    boost::shared_ptr<EventScheduler::Event> event(new Event(this, species));
    return event;
}

void LatticeSimulator::step()
{
    if (!is_initialized_)
    {
        initialize();
    }

    EventScheduler::value_type const& top(scheduler_.top());
    Real time(top.second->time());
    top.second->fire(); // top.second->time_ is updated in fire()
    (*world_).set_t(time);
    scheduler_.update(top);

    /*
    boost::shared_ptr<GSLRandomNumberGenerator> rng((*world_).rng());

    std::vector<Species> species((*world_).list_species());
    for (std::vector<Species>::iterator s_itr(species.begin());
            s_itr != species.end(); ++s_itr)
    {
        MolecularTypeBase* mt((*world_).get_molecular_type(*s_itr));
        //shuffle(*rng, mt->voxels());
        for (MolecularType::container_type::iterator itr(mt->begin());
                itr != mt->end(); ++itr)
        {
            Coord from_coord((*itr).first);
            Coord to_coord((*world_).get_neighbor(from_coord,
                        (*rng).uniform_int(0, 11)));
            (*world_).move(from_coord, to_coord);
        }
    }

    (*world_).set_t(t() + dt());
    */
}

bool LatticeSimulator::step(const Real& upto)
{
    if (!is_initialized_)
    {
        initialize();
    }

    /*
    Integer count((Integer)(upto / dt()) - num_steps());
    for (Integer i(0); i < count; ++i)
    {
        step();
    }
    */

    if (upto < t())
    {
        return false;
    }

    EventScheduler::value_type next_event;
    while ((next_event = scheduler_.top()).second->time() < upto)
    {
        Real time(next_event.second->time());
        next_event.second->fire(); // top.second->time_ is updated in fire()
        (*world_).set_t(time);
        scheduler_.update(next_event);
    }

    return true;
}

} // lattice

} // ecell4
