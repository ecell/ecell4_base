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
//     std::cerr << "form_domain_2D: forming domain for particle " << pid << std::endl;

    // -----------------------------------------------------------------------
    // form_multi
    const Real multi_radius = p.radius() * MULTI_SHELL_FACTOR;

    // list 2D domains within the multi shell
    std::vector<DomainID> intruders;
    for(const auto& item : shells_.list_shells_within_radius_2D(
                std::make_pair(p.position(), fid), multi_radius))
    {
        const auto& shell = item.first.second;
        const auto  did   = shell.domain_id().get();
        if(std::find(intruders.begin(), intruders.end(), did) == intruders.end())
        {
//             std::cerr << "intruder found: Shell " << item.first.first << " in " << did << std::endl;
            intruders.push_back(did);
        }
    }

//     std::cerr << "form_domain_2D: " << intruders.size() << " intrusive shells found" << std::endl;
    // list 3D multi domains that overlap with the shell.

    {
        // TODO: consider overlap between 3D domains and 2D domains
        // 1. find domains that overlaps with a face (it should be a multi)
        // 2. calculate the cross section of 3D domain with the plane on which the face belongs.
        // 3. check overlap between 2D domain and cross section of 3D domain.
    }

    if(intruders.empty())
    {
        // form a new multi!
        const auto did = didgen_();
        const auto sid = sidgen_();

        Shell sh(CircularShell(Circle(multi_radius, p.position(),
                        this->polygon().triangle_at(fid).normal()), fid), did);

//         std::cerr << "shell created: " << sid << std::endl;
        this->shells_.update_shell(sid, sh);
//         std::cerr << "shell container updated" << std::endl;

        MultiDomain dom(this->t());
        dom.add_particle(pid);
        dom.add_shell(sid);
        dom.determine_parameters(*(this->model()), *(this->world()));

//         std::cerr << "multi domain created " << did << std::endl;

        // add event with the same domain ID
        const auto evid = this->scheduler_.add(
                std::make_shared<event_type>(this->t() + dom.dt(), did));

//         std::cerr << "event added" << std::endl;

        // update begin_time and re-insert domain into domains_ container
        this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

//         std::cerr << "domain container updated" << std::endl;
        return;
    }

//     std::cerr << "intruder found. merge all domains" << std::endl;
    // XXX Currently all the domains are multi domains. merge all those multi domains.
    assert(false);

    const auto host_id = intruders.back();
    intruders.pop_back();

    const auto sid = sidgen_();
    Shell sh(CircularShell(Circle(multi_radius, p.position(),
                    this->polygon().triangle_at(fid).normal()), fid), host_id);
//     std::cerr << "shell created " << sid << std::endl;
    this->shells_.update_shell(sid, sh);
//     std::cerr << "shell inserted" << std::endl;

    assert(domains_.at(host_id).second.is_multi());
    auto& host = domains_.at(host_id).second.as_multi();
    host.add_particle(pid);
    host.add_shell(sid);
    for(const auto& did : intruders)
    {
        const auto dom_iter = domains_.find(did);

        const auto evid = dom_iter->second.first;
        scheduler_.remove(evid);

        assert(dom_iter->second.second.is_multi());
        const auto& dom = dom_iter->second.second.as_multi();
        for(const auto& pid : dom.particle_ids())
        {
            host.add_particle(pid);
        }
        for(const auto& sid : dom.shell_ids())
        {
            this->shells_.at(sid).second.domain_id() = host_id;
            host.add_shell(sid);
        }
        domains_.erase(dom_iter);
    }
    host.determine_parameters(*(this->model()), *(this->world()));

    return ;
}
void NGFRDSimulator::form_domain_3D(const ParticleID& pid, const Particle& p)
{
    // TODO: Currently we always form a multi domain.
//     std::cerr << "form_domain_3D: forming domain for particle " << pid << std::endl;

    // -----------------------------------------------------------------------
    // form_multi
    const Real multi_radius = p.radius() * MULTI_SHELL_FACTOR;

    // list 3D domains within the multi shell
    std::vector<DomainID> intruders;
    for(const auto& item : shells_.list_shells_within_radius_3D(
                p.position(), multi_radius))
    {
        const auto& shell = item.first.second;
        const auto  did   = shell.domain_id().get();
        if(std::find(intruders.begin(), intruders.end(), did) == intruders.end())
        {
//             std::cerr << "intruder found: Shell " << item.first.first << " in " << did << std::endl;
            intruders.push_back(did);
        }
    }
    // list 3D multi domains that overlap with the shell.

    // XXX: todo
    {
        // 0. check if multi domain overlaps with a face
        // 1. find domains are on the face
        // 2. calculate the cross section of 3D domain with the plane on which the face belongs.
        // 3. check overlap between 2D domain and cross section of 3D domain.
    }

    if(intruders.empty())
    {
        // form a new multi!
        const auto did = didgen_();
        const auto sid = sidgen_();

        Shell sh(SphericalShell(Sphere(p.position(), multi_radius)), did);
//         std::cerr << "shell created: " << sid << std::endl;
        this->shells_.update_shell(sid, sh);
//         std::cerr << "shell container updated" << std::endl;

        MultiDomain dom(this->t());
        dom.add_particle(pid);
        dom.add_shell(sid);
        dom.determine_parameters(*(this->model()), *(this->world()));

//         std::cerr << "multi domain created " << did << std::endl;

        // add event with the same domain ID
        const auto evid = this->scheduler_.add(
                std::make_shared<event_type>(this->t() + dom.dt(), did));

//         std::cerr << "event added" << std::endl;

        // update begin_time and re-insert domain into domains_ container
        this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

//         std::cerr << "domain container updated" << std::endl;
//         std::cerr << "all done." << std::endl;
        return;
    }
    assert(false);

    // XXX Currently all the domains are multi domains. merge all those multi domains.
//     std::cerr << "intruder found. merge all domains" << std::endl;

    const auto host_id = intruders.back();
    intruders.pop_back();

    const auto sid = sidgen_();
    Shell sh(SphericalShell(Sphere(p.position(), multi_radius)), host_id);
//     std::cerr << "shell created " << sid << std::endl;
    this->shells_.update_shell(sid, sh);
//     std::cerr << "shell inserted" << std::endl;

    assert(domains_.at(host_id).second.is_multi());
    auto& host = domains_.at(host_id).second.as_multi();
    host.add_particle(pid);
    host.add_shell(sid);
    for(const auto& did : intruders)
    {
        const auto dom_iter = domains_.find(did);

        const auto evid = dom_iter->second.first;
        scheduler_.remove(evid);

        assert(dom_iter->second.second.is_multi());
        const auto& dom = dom_iter->second.second.as_multi();
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
    host.determine_parameters(*(this->model()), *(this->world()));

    return;
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::fire_multi(const DomainID& did, MultiDomain dom)
{
//     std::cerr << "firing multi: " << did << std::endl;
//     for(const auto& sid : dom.shell_ids())
//     {
//         std::cerr << "included shell = " << sid << std::endl;
//     }

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
                // it performed stepping. the next event is at t+dt
                std::make_shared<event_type>(this->t() + dom.dt(), did));

        // update begin_time and re-insert domain into domains_ container
        dom.begin_time() = this->t() + dom.dt();
        domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

        return {/* All particles are re-inserted as Multi! */};
    }
    // something happens. remove multi and assign domains for each particles.

    // remove shells
    for(const auto& sid : dom.shell_ids())
    {
//         std::cerr << "removing shell " << sid << std::endl;
        this->shells_.remove_shell(sid);
    }

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> retval;
    for(const auto& pid : dom.particle_ids())
    {
        retval.push_back(world_->get_particle(pid));
    }
    return retval;
}

} // ngfrd
} // ecell4
