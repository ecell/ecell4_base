#include <ecell4/ngfrd/NGFRDSimulator.hpp>

namespace ecell4
{
namespace ngfrd
{
constexpr Real NGFRDSimulator::SAFETY;
constexpr Real NGFRDSimulator::SINGLE_SHELL_FACTOR;
constexpr Real NGFRDSimulator::MULTI_SHELL_FACTOR;
constexpr Real NGFRDSimulator::DEFAULT_DT_FACTOR;
constexpr Real NGFRDSimulator::CUTOFF_FACTOR;

void NGFRDSimulator::form_domain_2D(const ParticleID& pid, const Particle& p)
{
    // TODO
}
void NGFRDSimulator::form_domain_3D(const ParticleID& pid, const Particle& p)
{
    // TODO
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::fire_multi(const DomainID& did, MultiDomain dom)
{
    dom.step(*(this->model_), *this, *(this->world_));

    // XXX: If no (reaction, escapement) happens, we don't need to break it
    //      down into independent domains. For the efficiency, it re-inserts
    //      this domain into scheduler and returns nothing.
    //          Since nothing is returned, no domain will be formed at the end
    //      of this step.
    if(dom.eventkind() == MultiDomain::EventKind::None)
    {
        // add event with the same domain ID
        const auto evid = scheduler_.add(
                std::make_shared<event_type>(this->t(), did));

        // update begin_time and re-insert domain into domains_ container
        dom.begin_time() = this->t();
        domains_[did] = std::make_pair(evid, Domain(std::move(dom));

        return {/* All particles are re-inserted as Multi! */};
    }
    // something happens. remove multi and assign domains for each particles.
    return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>(
            dom.particles().begin(), dom.particles().end());
}

} // ngfrd
} // ecell4
