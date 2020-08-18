#include <ecell4/ngfrd/MultiDomain.hpp>
// #include <ecell4/ngfrd/NGFRDSimulator.hpp>

namespace ecell4
{
namespace ngfrd
{


void MultiDomain::step(const Model& model, NGFRDSimulator& sim, NGFRDWorld& world, const Real dt)
{
    assert(this->dt_              > 0.0);
    assert(this->reaction_length_ > 0.0);

    this->last_reactions_.clear();
    this->kind_ = EventKind::None;

    BDPropagator propagator(model, world, sim, sim.random_number_generator(),
            this->dt_, this->max_retry, this->particle_ids_, this->shell_ids_,
            last_reactions);

    while(propagator())
    {
        if(!last_reactions_.empty())
        {
            kind_ = REACTION;
        }
        if(propagator.vc().escaped())
        {
            kind_ = ESCAPE;
        }
    }
    return;
}

} // ngfrd
} // ecell4
