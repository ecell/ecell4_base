#include "SGFRDSimulator.hpp"

namespace ecell4
{
namespace sgfrd
{
const Real SGFRDSimulator::single_circular_shell_factor        = 1.5;
const Real SGFRDSimulator::single_circular_shell_mergin        = 1.0 - 1e-7;
const Real SGFRDSimulator::single_conical_surface_shell_factor = 1.5;
const Real SGFRDSimulator::single_conical_surface_shell_mergin = 1.0 - 1e-7;

void SGFRDSimulator::fire_single(const Single& dom, DomainID did)
{
    SGFRD_LOG(trace, "single fired");
    const ShellID sid(dom.shell_id());
    ParticleID pid; Particle p; FaceID fid;
    switch(dom.eventkind())
    {
    case Single::ESCAPE:
    {
        SGFRD_LOG(trace, "single escape");
        boost::tie(pid, p, fid) = boost::apply_visitor(make_visitor(
            resolve<const circular_shell_type&,
                    boost::tuple<ParticleID, Particle, FaceID> >(boost::bind(
                &self_type::escape_single<circular_shell_type>,
                this, _1, dom)),
            resolve<const conical_surface_shell_type&,
                    boost::tuple<ParticleID, Particle, FaceID> >(boost::bind(
                &self_type::escape_single<conical_surface_shell_type>,
                this, _1, dom))
            ), get_shell(sid));
        this->remove_shell(sid);
        this->create_event(pid, p, fid);
        return;
    }
    case Single::REACTION:
    {
        SGFRD_LOG(trace, "single reaction");
        BOOST_AUTO(results, boost::apply_visitor(make_visitor(
            resolve<const circular_shell_type&,
                    boost::container::static_vector<
                        boost::tuple<ParticleID, Particle, FaceID>, 2>
                    >(boost::bind(
                &self_type::reaction_single<circular_shell_type>,
                this, _1, dom, did)),
            resolve<const conical_surface_shell_type&,
                    boost::container::static_vector<
                        boost::tuple<ParticleID, Particle, FaceID>, 2>
                    >(boost::bind(
                &self_type::reaction_single<conical_surface_shell_type>,
                this, _1, dom, did))
            ), get_shell(sid)));
        this->remove_shell(sid);

        BOOST_FOREACH(boost::tie(pid, p, fid), results)
        {
            this->create_event(pid, p, fid);
        }
        return;
    }
    case Single::UNKNOWN:
        throw std::logic_error("when firing Single: event unspecified");
    default:
        throw std::logic_error("when firing Single: invalid enum value");
    }
}

bool SGFRDSimulator::burst_and_shrink_overlaps(
        const Particle& p, const face_id_type& fid)
{
    SGFRD_LOG(trace, "burst_and_shrink_overlaps");
    const Real tm = this->time();
    BOOST_AUTO(intruders, this->get_intrusive_domains(
               std::make_pair(p.position(), fid), p.radius()));

    bool no_overlap = true;
    DomainID did;
    BOOST_FOREACH(boost::tie(did, boost::tuples::ignore), intruders)
    {
        ParticleID pid_; Particle p_; FaceID fid_;
        BOOST_FOREACH(boost::tie(pid_, p_, fid_),
                      burst_event(std::make_pair(did, pickout_event(did)), tm))
        {
            const Real dist = distance(p.position(), p_.position());
            no_overlap = no_overlap && (dist > p.radius() + p_.radius());
            add_event(create_closely_fitted_domain(
                create_closely_fitted_shell(pid_, p_, fid_), pid_, p_));
        }
    }
    return no_overlap;
}

/*XXX assuming doms are multi or shrinked single domain and are sorted by dist*/
DomainID SGFRDSimulator::form_multi(
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const std::vector<std::pair<DomainID, Real> >& doms)
{
    SGFRD_LOG(trace, "form_multi called");
    // make Multi domain to join if the closest domain is no Multi.
    DomainID id_of_multi_to_join = doms.front().first;
    BOOST_AUTO(closest_dom, pickout_event(doms.front().first)->domain());
    if(closest_dom.which() != event_type::idx_multi)
    {
        // 1. closest domain is shrinked single. create empty multi
        Multi new_multi = create_empty_multi();
        // 2. then add a particle in the closest domain creating min_shell.
        Single& shrinked = boost::get<Single>(closest_dom);
        ParticleID pid_;
        boost::tie(pid_, boost::tuples::ignore) = shrinked.particle_id_pair();
        ShellID sid = shrinked.shell_id();
        new_multi.add_particle(pid);
        new_multi.add_shell(sid);
        // 3. and remove the event and shell of shrinked shell.
        remove_shell(sid);
        remove_event(doms.front().first);
        // 4. at last, assign the domain to event scheduler.
        id_of_multi_to_join = add_event(new_multi);
    }
    Multi& multi_to_join = boost::get<Multi>(
            pickout_event(id_of_multi_to_join)->domain());

    // create min shell for particle passed as an argument.
    BOOST_AUTO(sh, create_minimum_single_shell(pid, p, fid));
    const bool add_pt_result = multi_to_join.add_particle(pid);
    const bool add_sh_result = multi_to_join.add_shell(sh.first);
    assert(add_pt_result);
    assert(add_sh_result);

    // lookup rest of domains and if there are domains to join in to the multi,
    // add them recursively.
    const Real new_shell_size = sh.second.size();
    DomainID did; Real dist;
    BOOST_FOREACH(boost::tie(did, dist), doms)
    {
        if(dist < new_shell_size)
        {
            form_multi_recursive(multi_to_join, did);
        }
    }
    return id_of_multi_to_join;
}

void SGFRDSimulator::form_multi_recursive(Multi& multi_to_join, const DomainID did)
{
    SGFRD_LOG(trace, "form_multi_recursive called");
    event_type& ev = *(pickout_event(did));
    if(ev.which_domain() == event_type::idx_multi)
    {
        merge_multi(boost::get<Multi>(ev.domain()), multi_to_join);
        return;
    }
    else if(ev.which_domain() == event_type::idx_pair)
    {
        throw std::logic_error("shrinked pair domain exists!");
    }

    Single& dom = boost::get<Single>(ev.domain());
    ParticleID pid; Particle p;
    boost::tie(pid, p) = dom.particle_id_pair();
    FaceID  fid = this->get_face_id(pid);

    remove_shell(dom.shell_id());
    remove_event(did);

    const Real new_shell_size = calc_min_single_circular_shell_radius(p);

    BOOST_AUTO(sh, create_single_circular_shell(
                       std::make_pair(p.position(), fid), new_shell_size));
    const bool add_pt_result = multi_to_join.add_particle(pid);
    const bool add_sh_result = multi_to_join.add_shell(sh.first);
    assert(add_pt_result);
    assert(add_sh_result);

    const std::vector<std::pair<DomainID, Real> > intrusive_domains(
            get_intrusive_domains(std::make_pair(p.position(), fid), new_shell_size));

    const std::vector<std::pair<DomainID, Real> > bursted_domains(
        burst_and_shrink_non_multis(pid, p, fid, intrusive_domains));

    DomainID did_; Real dist;
    BOOST_FOREACH(boost::tie(did_, dist), bursted_domains)
    {
        if(dist < new_shell_size)
        {
            form_multi_recursive(multi_to_join, did_);
        }
    }
    return;
}

DomainID SGFRDSimulator::create_event(
            const ParticleID& pid, const Particle& p, const face_id_type fid)
{
    SGFRD_LOG(trace, "create_event called");
    const std::pair<Real3, face_id_type> pos = std::make_pair(p.position(), fid);

    const Real min_circle_size = p.radius() * single_circular_shell_factor;
          Real max_circle_size = get_max_circle_size(pos);

    SGFRD_LOG(debug, boost::format("min_circle_size = %1%, max_circle_size = %2%")
              % min_circle_size % max_circle_size);

    if(max_circle_size < min_circle_size)// draw conical shell
    {
        SGFRD_LOG(trace, "drawing conical shell");
        const std::vector<std::pair<vertex_id_type, Real> > intrusive_vertices(
                get_intrusive_vertices(pos, std::numeric_limits<Real>::infinity()));

        const vertex_id_type& vid = intrusive_vertices.front().first;
        SGFRD_LOG(debug, boost::format("vertex id = %1%, distance = %2%")
                  % vid % intrusive_vertices.front().second);

        const Real min_cone_size =
            (p.radius() + intrusive_vertices.front().second) *
            single_conical_surface_shell_factor;
        const Real max_cone_size = get_max_cone_size(vid);
        SGFRD_LOG(debug, boost::format("min_cone_size = %1%, max_cone_size = %2%")
                  % min_cone_size % max_cone_size);

        const std::vector<std::pair<DomainID, Real> > intrusive_domains(
                get_intrusive_domains(vid, max_cone_size));
        SGFRD_LOG(debug, boost::format("intrusive_domain_size = %1%")
                  % intrusive_domains.size());

        if(intrusive_domains.empty())
        {
            return add_event(create_single(create_single_conical_surface_shell(
                                           vid, max_cone_size), pid, p));
        }

        if(intrusive_domains.front().second > min_cone_size)
        {
            SGFRD_LOG(trace, "intrusive domains exist but enough distant");
            return add_event(create_single(create_single_conical_surface_shell(
                vid, intrusive_domains.front().second *
                single_conical_surface_shell_mergin), pid, p));
        }

        // burst intruder_domains and get new positions of particles
        std::vector<std::pair<DomainID, Real> > shrinked_or_multi =
            burst_and_shrink_non_multis(pid, p, fid, intrusive_domains);
        SGFRD_LOG(trace, "close domains are bursted.");

        if(shrinked_or_multi.front().second > min_cone_size)
        {
            SGFRD_LOG(trace, "after burst, no intruders exist");
            return add_event(create_single(create_single_conical_surface_shell(
                vid, shrinked_or_multi.front().second *
                     single_conical_surface_shell_mergin),
                pid, p));
        }

        SGFRD_LOG(trace, "forming multi");
        return form_multi(pid, p, fid, shrinked_or_multi);
    }

    SGFRD_LOG(trace, "drawing circular shell");

    const std::vector<std::pair<DomainID, Real> > intrusive_domains(
            get_intrusive_domains(pos, max_circle_size));
    SGFRD_LOG(debug, boost::format("intrusive_domain_size = %1%")
              % intrusive_domains.size());

    if(intrusive_domains.empty())
    {
        SGFRD_LOG(trace, "no intrusive domains exists");
        return add_event(create_single(create_single_circular_shell(
                         pos, max_circle_size), pid, p));
    }

    if(intrusive_domains.front().second > min_circle_size)
    {
        SGFRD_LOG(trace, "intrusive domains exists but enough distant");
        return add_event(create_single(create_single_circular_shell(
            pos, intrusive_domains.front().second * single_circular_shell_mergin),
            pid, p));
    }

    std::vector<std::pair<DomainID, Real> > shrinked_or_multi =
        burst_and_shrink_non_multis(pid, p, fid, intrusive_domains);
    SGFRD_LOG(trace, "intruder domains are bursted");

    if(shrinked_or_multi.front().second > min_circle_size)
    {
        SGFRD_LOG(trace, "after burst, no intruders exist");
        return add_event(create_single(create_single_circular_shell(
            pos, shrinked_or_multi.front().second * single_circular_shell_mergin),
            pid, p));
    }

    SGFRD_LOG(trace, "forming multi");
    return form_multi(pid, p, fid, shrinked_or_multi);
}

} // sgfrd
} // ecell4
