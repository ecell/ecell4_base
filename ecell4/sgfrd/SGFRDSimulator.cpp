#include "SGFRDSimulator.hpp"

namespace ecell4
{
namespace sgfrd
{
const Real SGFRDSimulator::single_circular_shell_factor        = 1.5;
const Real SGFRDSimulator::single_circular_shell_mergin        = 1.0 - 1e-7;
const Real SGFRDSimulator::single_conical_surface_shell_factor = 1.5;
const Real SGFRDSimulator::single_conical_surface_shell_mergin = 1.0 - 1e-7;
const Real SGFRDSimulator::reaction_length                     = 1e-5;

void SGFRDSimulator::fire_single(const Single& dom, DomainID did)
{
    SGFRD_SCOPE(us, fire_single, tracer_);
    SGFRD_TRACE(tracer_.write("fire single domain %1%", did))

    const ShellID sid(dom.shell_id());
    ParticleID pid; Particle p; FaceID fid;
    switch(dom.eventkind())
    {
    case Single::ESCAPE:
    {
        SGFRD_SCOPE(us, single_escape, tracer_);
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
        SGFRD_TRACE(tracer_.write("shell %1% removed", sid))

        SGFRD_TRACE(tracer_.write("adding next event for %1%", pid))
        this->create_event(pid, p, fid);
        return;
    }
    case Single::REACTION:
    {
        SGFRD_SCOPE(us, single_reaction, tracer_);

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
            SGFRD_TRACE(tracer_.write("adding next event for %1%", pid))
            add_event(create_closely_fitted_domain(create_closely_fitted_shell(
                      pid, p, fid), pid, p));
//             this->create_event(pid, p, fid); to consider two domains evenly
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
        const Particle& p, const face_id_type& fid, const DomainID& did)
{
    SGFRD_SCOPE(us, burst_and_shrink_overlaps, tracer_);
    const Real tm = this->time();
    BOOST_AUTO(intruders, this->get_intrusive_domains(
               std::make_pair(p.position(), fid), p.radius()));
    SGFRD_TRACE(tracer_.write("there are %1% intruders", intruders.size()))

    bool no_overlap = true;
    DomainID did_;
    BOOST_FOREACH(boost::tie(did_, boost::tuples::ignore), intruders)
    {
        SGFRD_SCOPE(ms, intruders, tracer_)
        SGFRD_TRACE(tracer_.write("burst domain %1%", did_))

        if(did == did_)
        {
            SGFRD_TRACE(tracer_.write("domain %1% was ignored", did_))
            continue;
        }

        ParticleID pid_; Particle p_; FaceID fid_;
        BOOST_FOREACH(boost::tie(pid_, p_, fid_),
                      burst_event(std::make_pair(did_, get_event(did_)), tm))
        {
            const Real dist = this->polygon().distance(
                std::make_pair(p.position(), fid),
                std::make_pair(p_.position(), fid_));
            no_overlap = no_overlap && (dist > p.radius() + p_.radius());
            add_event(create_closely_fitted_domain(
                create_closely_fitted_shell(pid_, p_, fid_), pid_, p_));
            scheduler_.remove(did_);
        }
    }
    return no_overlap;
}

/*XXX assuming doms are multi or shrinked single domain and are sorted by dist*/
DomainID SGFRDSimulator::form_multi(
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const std::vector<std::pair<DomainID, Real> >& doms)
{
    SGFRD_SCOPE(us, form_multi, tracer_);
    bool skip_first = false;
    DomainID formed_multi_id;
    if(get_event(doms.front().first)->which_domain() == event_type::multi_domain)
    {
        SGFRD_TRACE(tracer_.write("closest intruder is a multi domain. add all to this"))
        skip_first = true;
        formed_multi_id = doms.front().first;
    }
    else
    {
        SGFRD_TRACE(tracer_.write("closest intruder is not a multi. make empty multi"))

        BOOST_AUTO(new_multi, create_empty_multi());
        formed_multi_id = add_event(new_multi);
        SGFRD_TRACE(tracer_.write("new multi(%1%) created", formed_multi_id))
    }
    const domain_id_setter didset(formed_multi_id);
    Multi& formed_multi = boost::get<Multi>(get_event(formed_multi_id)->domain());

    BOOST_AUTO(minsh, create_minimum_single_shell(pid, p, fid));
    const Real new_shell_radius = minsh.second.size();
    formed_multi.add_particle(pid);
    formed_multi.add_shell(minsh.first);

    SGFRD_TRACE(tracer_.write("particle (%1%) and shell (%2%) is added to multi",
                pid, minsh.first));

    DomainID did; Real dist;
    BOOST_FOREACH(boost::tie(did, dist), doms)
    {
        SGFRD_SCOPE(us, intruder_domain, tracer_);
        SGFRD_TRACE(tracer_.write("for domain %1%", did));
        if(skip_first){skip_first = false; continue;}

        if(dist < new_shell_radius) // add the domain to new multi
        {
            BOOST_AUTO(ev, get_event(did));
            if(ev->which_domain() == event_type::multi_domain)
            {
                SGFRD_TRACE(tracer_.write("domain (%1%) is multi. merging it...", did))
                merge_multi(boost::get<Multi>(ev->domain()), formed_multi);
            }
            else if(ev->which_domain() == event_type::pair_domain)
            {
                throw std::logic_error("pair event cannot join to multi");
            }
            else
            {
                SGFRD_TRACE(tracer_.write("domain (%1%) is single. adding it...", did))

                // update shell with min_single_circular_shell!
                ParticleID pid_; Particle p_;
                boost::tie(pid_, p_) =
                    boost::get<Single>(ev->domain()).particle_id_pair();
                const ShellID sid = boost::get<Single>(ev->domain()).shell_id();
                SGFRD_TRACE(tracer_.write("domain (%1%) has particle(%2%), shell(%3%)",
                            did, pid_, sid));

                // edit shell size to be min_shell_radius.
                circular_shell_type clsh =
                    boost::get<circular_shell_type>(get_shell(sid));
                clsh.shape().size() =
                    calc_min_single_circular_shell_radius(p_);
                update_shell(sid, clsh, clsh.structure_id());
                SGFRD_TRACE(tracer_.write("shell(%1%) size updated to %2%.",
                            sid, clsh.shape().size()));

                formed_multi.add_particle(pid_);
                formed_multi.add_shell(sid);

                remove_event(did);
            }
        }
    }
    mut_sh_vis_applier(didset, formed_multi);

    // search intruders on the new multi, burst them all and add to multi if needed
    add_to_multi_recursive(formed_multi);
    return formed_multi_id;
}

void SGFRDSimulator::add_to_multi_recursive(Multi& multi_to_join)
{
    SGFRD_SCOPE(us, add_to_multi_recursive, tracer_);

    const Real tm = this->time();
    bool multi_enlarged = false;
    const DomainID multi_to_join_id = get_domain_id(multi_to_join);
    const domain_id_setter didset(multi_to_join_id);

    BOOST_FOREACH(ShellID sid, multi_to_join.shell_ids())
    {
        // assuming multi has only a circular_shell...
        BOOST_AUTO(const& sh, boost::get<circular_shell_type>(get_shell(sid)));
        BOOST_AUTO(sh_pos, std::make_pair(sh.position(), sh.structure_id()));
        BOOST_AUTO(const intruder, get_intrusive_domains(
                   std::make_pair(sh.position(), sh.structure_id()), sh.size()));
        SGFRD_TRACE(tracer_.write("intrusive domains on shell(%1%) are collected(size = %2%)",
                    sid, intruder.size()));

        DomainID did;
        BOOST_FOREACH(boost::tie(did, boost::tuples::ignore), intruder)
        {
            if(did == multi_to_join_id) continue;

            SGFRD_TRACE(tracer_.write("bursting domain(%1%)", did));
            BOOST_AUTO(ev, get_event(did));

            if(ev->which_domain() == event_type::multi_domain)
            {
                merge_multi(boost::get<Multi>(ev->domain()), multi_to_join);
                multi_enlarged = true;
            }
            else
            {
                ParticleID pid; Particle p; FaceID fid;
                BOOST_FOREACH(boost::tie(pid, p, fid),
                              burst_event(std::make_pair(did, ev), tm))
                {
                    const Real dist = this->polygon().distance(
                            sh_pos, std::make_pair(p.position(), fid)) - sh.size();
                    const Real min_shell_radius =
                        calc_min_single_circular_shell_radius(p);
                    if(dist < min_shell_radius)
                    {
                        // add the particle into multi
                        BOOST_AUTO(minsh, create_minimum_single_shell(pid, p, fid));
                        multi_to_join.add_particle(pid);
                        multi_to_join.add_shell(minsh.first);
                        multi_enlarged = true;
                    }
                    else // enough distant. add closely-fitted shell
                    {
                        add_event(
                            create_closely_fitted_domain(
                                create_closely_fitted_shell(
                                    pid, p, this->get_face_id(pid)
                                    ), pid, p
                                )
                            );
                    }
                }
                remove_event(did);
            }
        }
    }
    if(multi_enlarged)
    {
        mut_sh_vis_applier(didset, multi_to_join);
        add_to_multi_recursive(multi_to_join);
    }

    return;
}

DomainID SGFRDSimulator::create_event(
            const ParticleID& pid, const Particle& p, const face_id_type fid)
{
    SGFRD_SCOPE(us, create_event, tracer_);
    const std::pair<Real3, face_id_type> pos = std::make_pair(p.position(), fid);

    const Real min_circle_size = p.radius() * single_circular_shell_factor;
          Real max_circle_size = get_max_circle_size(pos);

    SGFRD_TRACE(tracer_.write("min_circle_size = %1%, max_circle_size = %2%",
                min_circle_size, max_circle_size));

    if(max_circle_size < min_circle_size)// draw conical shell
    {
        SGFRD_SCOPE(us, draw_conical_shell, tracer_);
        const std::vector<std::pair<vertex_id_type, Real> > intrusive_vertices(
                get_intrusive_vertices(pos, std::numeric_limits<Real>::infinity()));

        const vertex_id_type& vid = intrusive_vertices.front().first;
        SGFRD_TRACE(tracer_.write("vertex id = %1%, distance = %2%",
                    vid, intrusive_vertices.front().second));

        const Real min_cone_size =
            (p.radius() + intrusive_vertices.front().second) *
            single_conical_surface_shell_factor;
        const Real max_cone_size = get_max_cone_size(vid);
        SGFRD_TRACE(tracer_.write("min_cone_size = %1%, max_cone_size = %2%",
                    min_cone_size, max_cone_size));

        const std::vector<std::pair<DomainID, Real> > intrusive_domains(
                get_intrusive_domains(vid, max_cone_size));
        SGFRD_TRACE(tracer_.write("intrusive_domain_size = %1%", intrusive_domains.size()));

        if(intrusive_domains.empty())
        {
            return add_event(create_single(create_single_conical_surface_shell(
                                           vid, max_cone_size), pid, p));
        }

        std::vector<std::pair<DomainID, Real> > min_shell_intruder;
        for(typename std::vector<std::pair<DomainID, Real> >::const_iterator
                iter = intrusive_domains.begin(), end = intrusive_domains.end();
                iter != end; ++iter)
            if(iter->second <= min_cone_size)
                min_shell_intruder.push_back(*iter);

        if(min_shell_intruder.empty())
        {
            SGFRD_TRACE(tracer_.write("intrusive domains exist but enough distant"));
#ifndef ECELL4_SGFRD_NO_TRACE
            for(std::size_t i=0; i<intrusive_domains.size(); ++i)
            {
                tracer_.write("domain %1%; dist = %2%;",
                    intrusive_domains[i].first, intrusive_domains[i].second);
            }
#endif//ECELL4_SGFRD_NO_TRACE

            const Real shell_size =
                std::min(max_cone_size, intrusive_domains.front().second) *
                single_conical_surface_shell_mergin;

            return add_event(create_single(
                        create_single_conical_surface_shell(vid, shell_size),
                        pid, p));
        }

        // burst intruder_domains and get new positions of particles
        std::vector<std::pair<DomainID, Real> > shrinked_or_multi =
            burst_and_shrink_non_multis(vid, min_shell_intruder);
        SGFRD_TRACE(tracer_.write("close domains are bursted."));

        if(shrinked_or_multi.front().second > min_cone_size)
        {
            SGFRD_TRACE(tracer_.write(
                "after burst, no intruders exist in the min-range %1%",
                min_cone_size));

#ifndef ECELL4_SGFRD_NO_TRACE
            for(std::size_t i=0; i<shrinked_or_multi.size(); ++i)
            {
                tracer_.write("domain %1%; dist = %2%;",
                    shrinked_or_multi[i].first, shrinked_or_multi[i].second);
            }
#endif//ECELL4_SGFRD_NO_TRACE
            const Real shell_size =
                std::min(max_cone_size, shrinked_or_multi.front().second) *
                single_conical_surface_shell_mergin;

            return add_event(create_single(
                        create_single_conical_surface_shell(vid, shell_size),
                        pid, p));
        }

        SGFRD_TRACE(tracer_.write("forming multi"));
        return form_multi(pid, p, fid, shrinked_or_multi);
    }
    SGFRD_SCOPE(us, draw_circular_shell, tracer_);

    const std::vector<std::pair<DomainID, Real> > intrusive_domains(
            get_intrusive_domains(pos, max_circle_size));
    SGFRD_TRACE(tracer_.write("intrusive_domain_size = %1%", intrusive_domains.size()))

    if(intrusive_domains.empty())
    {
        SGFRD_TRACE(tracer_.write("no intrusive domains exists"))
        return add_event(create_single(create_single_circular_shell(
                         pos, max_circle_size), pid, p));
    }

    std::vector<std::pair<DomainID, Real> > min_shell_intruder;
    for(typename std::vector<std::pair<DomainID, Real> >::const_iterator
            iter = intrusive_domains.begin(), end = intrusive_domains.end();
            iter != end; ++iter)
        if(iter->second <= min_circle_size)
            min_shell_intruder.push_back(*iter);

    if(min_shell_intruder.empty())
    {
        SGFRD_TRACE(tracer_.write("intrusive domains exists but enough distant"))
#ifndef ECELL4_SGFRD_NO_TRACE
        for(std::size_t i=0; i<intrusive_domains.size(); ++i)
        {
            tracer_.write("domain %1%; dist = %2%;",
                intrusive_domains[i].first, intrusive_domains[i].second);
        }
#endif//ECELL4_SGFRD_NO_TRACE

        const Real shell_size =
            std::min(max_circle_size, intrusive_domains.front().second) *
            single_circular_shell_mergin;

        return add_event(create_single(
                    create_single_circular_shell(pos, shell_size), pid, p));
    }

    std::vector<std::pair<DomainID, Real> > shrinked_or_multi =
        burst_and_shrink_non_multis(pid, p, fid, min_shell_intruder);
    SGFRD_TRACE(tracer_.write("min_shell_intruder domains are bursted"))

    if(shrinked_or_multi.front().second > min_circle_size)
    {
        SGFRD_TRACE(tracer_.write("after burst, no intruders exist"))
#ifndef ECELL4_SGFRD_NO_TRACE
        for(std::size_t i=0; i<shrinked_or_multi.size(); ++i)
        {
            tracer_.write("domain %1%; dist = %2%;",
                shrinked_or_multi[i].first, shrinked_or_multi[i].second);
        }
#endif//ECELL4_SGFRD_NO_TRACE

        const Real shell_size =
            std::min(max_circle_size, shrinked_or_multi.front().second) *
            single_circular_shell_mergin;

        return add_event(create_single(
                    create_single_circular_shell(pos, shell_size), pid, p));
    }

    SGFRD_TRACE(tracer_.write("forming multi"))
    return form_multi(pid, p, fid, shrinked_or_multi);
}

} // sgfrd
} // ecell4
