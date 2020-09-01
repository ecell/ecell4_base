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

void NGFRDSimulator::form_domain_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    // TODO: Currently we always form a multi domain.

    // -----------------------------------------------------------------------
    // form_multi
    const Real multi_radius = p.radius() * MULTI_SHELL_FACTOR;

    // list 2D domains within the multi shell
    std::vector<DomainID> intruders;
    for(const auto& item : shells_.list_shells_within_radius_2D(
                std::make_pair(p.position, fid), multi_radius))
    {
        const auto& shell = item.first.second;
        const auto  did   = shell.domain_id().get();
        if(std::find(intruders.begin(), intruders.end(), did) == intruders.end())
        {
            intruders.push_back(did);
        }
    }
    // list 3D multi domains that overlap with the shell.

    // XXX: todo
    // 1. find domains that overlaps with a face (it should be a multi)
    // 2. calculate the cross section of 3D domain with the plane on which the face belongs.
    // 3. check overlap between 2D domain and cross section of 3D domain.

    if(intruders.empty())
    {
        // form a new multi!
        const auto did = didgen_();
        const auto sid = sidgen_();

        Shell sh(CircularShell(Circle(multi_radius, p.position(),
                        polygon_->triangle_at(fid).normal()), fid), did);
        this->shells_.update_shell(sid, sh);

        MultiDomain dom(this->t());
        dom.add_particle(pid);
        dom.add_shell(sid);
        dom.determine_parameters(*(this->model()), *(this->world()));

        // add event with the same domain ID
        const auto evid = this->scheduler_.add(
                std::make_shared<event_type>(this->t(), did));

        // update begin_time and re-insert domain into domains_ container
        this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));
        return;
    }

    // XXX Currently all the domains are multi domains. merge all those multi domains.

    const auto host_id = intruders.back();
    intruders.pop_back();

    const auto sid = sidgen_();
    Shell sh(CircularShell(Circle(multi_radius, p.position(),
                    polygon_->triangle_at(fid).normal()), fid), did);
    this->shells_.update_shell(sid, sh);

    auto& host = domains_.at(host_id).second;
    assert(host.is_multi());
    host.add_particle(pid);
    host.add_shell(sid);
    for(const auto& did : intruders)
    {
        const auto dom_iter = domains_.find(did);

        const auto evid = dom_iter->second.first;
        scheduler_.remove(evid);

        const auto& dom = dom_iter->second.second;
        assert(dom.is_multi());
        for(const auto& pid : dom.particle_ids())
        {
            host.add_particle(pid);
        }
        for(const auto& sid : dom.shell_ids())
        {
            host.add_shell(sid);
        }
        domains_.erase(dom_iter);
    }
    dom.determine_parameters(*(this->model()), *(this->world()));

    return ;
}
void NGFRDSimulator::form_domain_3D(const ParticleID& pid, const Particle& p)
{
    // TODO: Currently we always form a multi domain.

    // -----------------------------------------------------------------------
    // form_multi
    const Real multi_radius = p.radius() * MULTI_SHELL_FACTOR;

    // list 3D domains within the multi shell
    std::vector<DomainID> intruders;
    for(const auto& item : shells_.list_shells_within_radius_3D(
                p.position, multi_radius))
    {
        const auto& shell = item.first.second;
        const auto  did   = shell.domain_id().get();
        if(std::find(intruders.begin(), intruders.end(), did) == intruders.end())
        {
            intruders.push_back(did);
        }
    }
    // list 3D multi domains that overlap with the shell.

    // XXX: todo
    // 0. check if multi domain overlaps with a face
    // 1. find domains are on the face
    // 2. calculate the cross section of 3D domain with the plane on which the face belongs.
    // 3. check overlap between 2D domain and cross section of 3D domain.

    if(intruders.empty())
    {
        // form a new multi!
        const auto did = didgen_();
        const auto sid = sidgen_();

        Shell sh(SphericalShell(Sphere(p.position(), multi_radius)), did);
        this->shells_.update_shell(sid, sh);

        MultiDomain dom(this->t());
        dom.add_particle(pid);
        dom.add_shell(sid);
        dom.determine_parameters(*(this->model()), *(this->world()));

        // add event with the same domain ID
        const auto evid = this->scheduler_.add(
                std::make_shared<event_type>(this->t(), did));

        // update begin_time and re-insert domain into domains_ container
        this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));
        return;
    }

    // XXX Currently all the domains are multi domains. merge all those multi domains.

    const auto host_id = intruders.back();
    intruders.pop_back();

    const auto sid = sidgen_();
    Shell sh(SphericalShell(Sphere(p.position(), multi_radius)), host_id);
    this->shells_.update_shell(sid, sh);

    auto& host = domains_.at(host_id).second;
    assert(host.is_multi());
    host.add_particle(pid);
    host.add_shell(sid);
    for(const auto& did : intruders)
    {
        const auto dom_iter = domains_.find(did);

        const auto evid = dom_iter->second.first;
        scheduler_.remove(evid);

        const auto& dom = dom_iter->second.second;
        assert(dom.is_multi());
        for(const auto& pid : dom.particle_ids())
        {
            host.add_particle(pid);
        }
        for(const auto& sid : dom.shell_ids())
        {
            host.add_shell(sid);
        }
        domains_.erase(dom_iter);
    }
    dom.determine_parameters(*(this->model()), *(this->world()));

    return;
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
        domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

        return {/* All particles are re-inserted as Multi! */};
    }
    // something happens. remove multi and assign domains for each particles.
    return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>(
            dom.particles().begin(), dom.particles().end());
}

} // ngfrd
} // ecell4
