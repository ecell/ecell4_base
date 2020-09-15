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
    ECELL4_NGFRD_LOG_FUNCTION();
    // TODO: Currently we always form a multi domain.
    ECELL4_NGFRD_LOG("form_domain_2D: forming domain for particle ", pid);

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
            ECELL4_NGFRD_LOG("intruder found: 2DShell ", item.first.first, " in ", did);
            intruders.push_back(did);
        }
    }
    ECELL4_NGFRD_LOG("form_domain_2D: ", intruders.size(), " intrusive shells found");

    // list 3D multi domains that overlap with the shell.
    // Here, we list 3D shells that are within the bounding sphere and overlap
    // with the polygon.
    for(const auto& item : shells_.list_shells_within_radius_3D(
                p.position(), multi_radius))
    {
        const auto& shell = item.first.second;
        const auto  did   = shell.domain_id().get();

        if(!this->domains_.at(did).second.is_multi())
        {
            // non-Multi domain never overlaps with the polygon (?)
            continue;
        }
        // Multi 3D Shell is always spherical.
        const auto& sh = shell.as_spherical();
        if(this->polygon().has_overlapping_faces(sh.shape().position(), sh.shape().radius()))
        {
            // This shell is within the bounding sphere of 2D shell and overlaps
            // with Polygon. Insert it to multi (that does not always mean that
            // the domain overlaps with 2D shell, but a nice approximation (sort of))
            if(std::find(intruders.begin(), intruders.end(), did) == intruders.end())
            {
                ECELL4_NGFRD_LOG("intruder found: 3DShell ", item.first.first, " in ", did);
                intruders.push_back(did);
            }
        }
    }

    if(intruders.empty())
    {
        // form a new multi!
        const auto did = didgen_();
        const auto sid = sidgen_();

        Shell sh(CircularShell(Circle(multi_radius, p.position(),
                        this->polygon().triangle_at(fid).normal()), fid), did);

        ECELL4_NGFRD_LOG("shell created: ", sid);
        this->shells_.update_shell(sid, sh);
        ECELL4_NGFRD_LOG("shell container updated");

        MultiDomain dom(this->t());
        dom.add_particle(pid);
        dom.add_shell(sid);
        dom.determine_parameters(*(this->model()), *(this->world()));

        ECELL4_NGFRD_LOG("multi domain created ", did);

        // add event with the same domain ID
        const auto evid = this->scheduler_.add(
                std::make_shared<event_type>(this->t() + dom.dt(), did));

        ECELL4_NGFRD_LOG("event added");

        // update begin_time and re-insert domain into domains_ container
        this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));
        return;
    }

    ECELL4_NGFRD_LOG("intruder found. merge all domains");
    // XXX Currently all the domains are multi domains. merge all those multi domains.
    // Later we need to burst domains and check if the resulting particle should
    // be in Multi or form Single

    const auto host_id = intruders.back();
    intruders.pop_back();

    const auto sid = sidgen_();
    Shell sh(CircularShell(Circle(multi_radius, p.position(),
                    this->polygon().triangle_at(fid).normal()), fid), host_id);
    ECELL4_NGFRD_LOG("shell created ", sid);
    this->shells_.update_shell(sid, sh);
    ECELL4_NGFRD_LOG("shell inserted");

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
    ECELL4_NGFRD_LOG_FUNCTION();
    // TODO: Currently we always form a multi domain.
    ECELL4_NGFRD_LOG("form_domain_3D: forming domain for particle ", pid);

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
            ECELL4_NGFRD_LOG("intruder found: Shell ", item.first.first, " in ", did);
            intruders.push_back(did);
        }
    }
    // list 3D multi domains that overlap with the shell.
    const auto pbc = this->world_->boundary();
    for(const auto fidpd : this->polygon().list_faces_within_radius(
                p.position(), multi_radius))
    {
        const auto fid = fidpd.first.first;

        // shells on this triangle
        if(const auto shids = shells_.shells_on(fid))
        {
            for(const auto shid : *shids)
            {
                // check overlap with the bounding sphere of 2D shell
                // as an approximation
                const auto sh = this->shells_.get_shell(shid).second;
                const auto bs = sh.bounding_sphere();
                const auto dist = length(
                        pbc.periodic_transpose(p.position(), bs.position()) -
                        bs.position());
                if(dist <= multi_radius + bs.radius())
                {
                    intruders.push_back(*sh.domain_id());
                }
            }
        }
        // shells on the vertices of the triangle
        for(const auto vid : this->polygon().vertices_of(fid))
        {
            if(const auto shids = shells_.shells_on(vid))
            {
                for(const auto shid : *shids)
                {
                    // check overlap with the bounding sphere of 2D shell
                    // as an approximation
                    const auto sh = this->shells_.get_shell(shid).second;
                    const auto bs = sh.bounding_sphere();
                    const auto dist = length(
                            pbc.periodic_transpose(p.position(), bs.position()) -
                            bs.position());
                    if(dist <= multi_radius + bs.radius())
                    {
                        intruders.push_back(*sh.domain_id());
                    }
                }
            }
        }
        // shells on the neighboring triangles
        for(const auto nfid : this->polygon().neighbor_faces_of(fid))
        {
            if(const auto shids = shells_.shells_on(nfid))
            {
                for(const auto shid : *shids)
                {
                    // check overlap with the bounding sphere of 2D shell
                    // as an approximation
                    const auto sh = this->shells_.get_shell(shid).second;
                    const auto bs = sh.bounding_sphere();
                    const auto dist = length(
                            pbc.periodic_transpose(p.position(), bs.position()) -
                            bs.position());
                    if(dist <= multi_radius + bs.radius())
                    {
                        intruders.push_back(*sh.domain_id());
                    }
                }
            }
        }
    }

    if(intruders.empty())
    {
        // form a new multi!
        const auto did = didgen_();
        const auto sid = sidgen_();

        Shell sh(SphericalShell(Sphere(p.position(), multi_radius)), did);
        ECELL4_NGFRD_LOG("shell created: ", sid);
        this->shells_.update_shell(sid, sh);
        ECELL4_NGFRD_LOG("shell inserted");

        MultiDomain dom(this->t());
        dom.add_particle(pid);
        dom.add_shell(sid);
        dom.determine_parameters(*(this->model()), *(this->world()));

        ECELL4_NGFRD_LOG("multi domain created ", did);

        // add event with the same domain ID
        const auto evid = this->scheduler_.add(
                std::make_shared<event_type>(this->t() + dom.dt(), did));
        ECELL4_NGFRD_LOG("event added: ", evid);

        // update begin_time and re-insert domain into domains_ container
        this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));
        ECELL4_NGFRD_LOG("domain container updated");
        return;
    }

    // XXX Currently all the domains are multi domains. merge all those multi domains.
    // Later we need to burst domains and check if the resulting particle should
    // be in Multi or form Single
    ECELL4_NGFRD_LOG("intruder found. merge all domains");

    const auto host_id = intruders.back();
    intruders.pop_back();

    const auto sid = sidgen_();
    Shell sh(SphericalShell(Sphere(p.position(), multi_radius)), host_id);
    ECELL4_NGFRD_LOG("shell created: ", sid);
    this->shells_.update_shell(sid, sh);
    ECELL4_NGFRD_LOG("shell inserted");

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
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("firing multi: ", did);
    ECELL4_NGFRD_LOG("included shells: ", dom.shell_ids());

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
        ECELL4_NGFRD_LOG("removing shell: ", sid);
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
