#include "SGFRDSimulator.hpp"

namespace ecell4
{
namespace sgfrd
{
const Real SGFRDSimulator::single_circular_shell_factor        = 1.5;
const Real SGFRDSimulator::single_circular_shell_mergin        = 1.0 - 1e-7;
const Real SGFRDSimulator::single_conical_surface_shell_factor = 1.5;
const Real SGFRDSimulator::single_conical_surface_shell_mergin = 1.0 - 1e-7;

boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::propagate_single_circular(
        const circular_shell_type& sh, const Single& dom, const Real tm)
{
    SGFRD_SCOPE(us, propagate_single_circular, tracer_);

    Particle   p   = dom.particle();
    ParticleID pid = dom.particle_id();

    greens_functions::GreensFunction2DAbsSym gf(p.D(), sh.size() - p.radius());

    const Real del_t = tm - dom.begin_time();

    SGFRD_TRACE(tracer_.write("delta t for domain having shell %1% is %2%",
                del_t, dom.shell_id()));
    SGFRD_TRACE(tracer_.write("its own dt = %1%, and begin_time = %2%",
                dom.dt(), dom.begin_time()));

    const Real r     = gf.drawR(this->uniform_real(), del_t);
    const Real theta = this->uniform_real() * 2 * M_PI;

    SGFRD_TRACE(tracer_.write("r = %1%, theta = %2%", r, theta));

    const FaceID         fid  = this->get_face_id(pid);
    const triangle_type& face = this->polygon().triangle_at(fid);
    const Real3 direction = rotate(theta, face.normal(), face.represent());
    const Real  len_direction = length(direction);

    SGFRD_TRACE(tracer_.write("direction = %1%, len = %2%", direction, len_direction));

    std::pair<std::pair<Real3, FaceID>, Real3> state =
        std::make_pair(std::make_pair(p.position(), fid),
                       direction * r / len_direction);

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%", state.first.first, state.first.second));

    std::size_t continue_count =
        ecell4::polygon::travel(this->polygon(), state.first, state.second, 2);
    if(continue_count == 0)
    {
        SGFRD_TRACE(tracer_.write("moving on face: precision lost"))
    }

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%, dsp = %3%, count = %4%",
                state.first.first, state.first.second, state.second, continue_count))

    p.position() = state.first.first;
    this->update_particle(pid, p, state.first.second);
    SGFRD_TRACE(tracer_.write("particle updated"))

    return boost::make_tuple(pid, p, state.first.second);
}

boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::propagate_single_conical(
    const conical_surface_shell_type& sh, const Single& dom, const Real tm)
{
    SGFRD_SCOPE(us, propagate_single_conical, tracer_);

    Particle         p   = dom.particle();
    const ParticleID pid = dom.particle_id();
    const FaceID     fid = this->get_face_id(pid);

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%", p.position(), fid));

    const Real r_max = sh.size() - p.radius();
    greens_functions::GreensFunction2DRefWedgeAbs
        gf(/* D   = */ p.D(),
           /* r0  = */ length(p.position() - sh.position()),
           /* a   = */ r_max,
           /* phi = */ sh.shape().apex_angle());

    const Real del_t = tm - dom.begin_time();
    const Real r     = gf.drawR(this->uniform_real(), del_t);
    const Real theta = gf.drawTheta(this->uniform_real(), r, del_t);

    SGFRD_TRACE(tracer_.write("r = %1%, theta = %2%", r, theta));

    const std::pair<Real3, FaceID> state =
        ecell4::polygon::roll(this->polygon(), std::make_pair(p.position(), fid),
                      sh.structure_id(), r, theta);
    SGFRD_TRACE(tracer_.write("propagateed : pos = %1%, fid = %2%",
                state.first, state.second));

    p.position() = state.first;
    this->update_particle(pid, p, state.second);
    SGFRD_TRACE(tracer_.write("particle updated"))

    return boost::make_tuple(pid, p, state.second);
}

boost::container::static_vector<
    boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>, 2>
SGFRDSimulator::reaction_single(
        const shell_type& sh, const Single& dom, const DomainID did)
{
    // considering shell remains
    SGFRD_SCOPE(us, reaction_single, tracer_)
    ParticleID pid; Particle p; FaceID fid;
    boost::tie(pid, p, fid) = this->propagate_single(sh, dom, this->time());
    return attempt_reaction_single(sh, did, dom, pid, p, fid);
}

boost::container::static_vector<SGFRDSimulator::pid_p_fid_tuple_type, 2>
SGFRDSimulator::attempt_reaction_single(
        const shell_type& sh,  const DomainID did, const Single& dom,
        const ParticleID& pid, const Particle&  p, const FaceID& fid)
{
    SGFRD_SCOPE(us, attempt_reaction_single, tracer_)

    BOOST_AUTO(const rules, this->model_->query_reaction_rules(p.species()));
    if(rules.empty())
    {
        SGFRD_TRACE(tracer_.write("rule is empty. return particle kept intact"))
        return boost::container::static_vector<pid_p_fid_tuple_type, 2>(1,
                boost::make_tuple(pid, p, fid));
    }
    const ReactionRule& rule = determine_reaction_rule(rules);

    switch(rule.products().size())
    {
        case 0:
        {
            SGFRD_TRACE(tracer_.write("degradation reaction occurs."))
            this->remove_particle(pid, fid);
            last_reactions_.push_back(std::make_pair(rule,
                make_degradation_reaction_info(this->time(), pid, p)));
            return boost::container::static_vector<pid_p_fid_tuple_type, 2>(0);
        }
        case 1:
        {
            SGFRD_TRACE(tracer_.write("attempting 1 to 1 reaction."))
            return attempt_reaction_1_to_1(rule, sh, did, dom, pid, p, fid);
        }
        case 2:
        {
            SGFRD_TRACE(tracer_.write("attempting 1 to 2 reaction."))
            return attempt_reaction_1_to_2(rule, sh, did, dom, pid, p, fid);
        }
        default: throw NotImplemented("SGFRD Single Reaction:"
            "more than two products from one reactant are not allowed");
    }
}

boost::container::static_vector<SGFRDSimulator::pid_p_fid_tuple_type, 2>
SGFRDSimulator::attempt_reaction_1_to_1(const ReactionRule& rule,
        const shell_type& sh,  const DomainID did, const Single& dom,
        const ParticleID& pid, const Particle&  p, const FaceID& fid)
{
    SGFRD_SCOPE(us, attempt_reaction_1_to_1, tracer_)

    const Species species_new =
        this->model_->apply_species_attributes(rule.products().front());
    const molecule_info_type mol_info =
        this->world_->get_molecule_info(species_new);
    const Real radius_new = mol_info.radius;
    const Real D_new      = mol_info.D;
    const Particle p_new(species_new, p.position(), radius_new, D_new);

    inside_checker is_inside_of(p.position(), p.radius(), fid, this->polygon());
    if(!boost::apply_visitor(is_inside_of, sh))
    {
        SGFRD_SCOPE(us, particle_goes_outside, tracer_)
        // particle goes outside of the shell. must clear the volume.
        const bool no_overlap = this->burst_and_shrink_overlaps(p_new, fid, did);
        SGFRD_TRACE(tracer_.write("no_overlap = %1%", no_overlap))
        if(!no_overlap)
        {// cannot avoid overlapping... reject the reaction.
            SGFRD_TRACE(tracer_.write("reject the reaction because of no space"))
            return boost::container::static_vector<pid_p_fid_tuple_type, 2>(
                    1, boost::make_tuple(pid, p, fid));
        }
    }

    SGFRD_TRACE(tracer_.write("reaction is accepted."))
    // reaction occurs. record this reaction to `last_reactions`.
    last_reactions_.push_back(std::make_pair(rule,
        make_unimolecular_reaction_info(this->time(), pid, p, pid, p_new)));

    // update particle and return resulting particles.
    this->update_particle(pid, p_new, fid);
    SGFRD_TRACE(tracer_.write("particle updated."))
    return boost::container::static_vector<pid_p_fid_tuple_type, 2>(
            1, boost::make_tuple(pid, p_new, fid));
}

boost::container::static_vector<SGFRDSimulator::pid_p_fid_tuple_type, 2>
SGFRDSimulator::attempt_reaction_1_to_2(const ReactionRule& rule,
        const shell_type& sh,  const DomainID did, const Single& dom,
        const ParticleID& pid, const Particle&  p, const FaceID& fid)
{
    SGFRD_SCOPE(us, attempt_reaction_1_to_2, tracer_)

    const Species sp1 =
        this->model_->apply_species_attributes(rule.products().at(0));
    const Species sp2 =
        this->model_->apply_species_attributes(rule.products().at(1));

    const molecule_info_type mol1 = world_->get_molecule_info(sp1);
    const molecule_info_type mol2 = world_->get_molecule_info(sp2);

    const Real D1(mol1.D),      D2(mol2.D);
    const Real r1(mol1.radius), r2(mol2.radius), r12(mol1.radius + mol2.radius);

    boost::array<std::pair<Real3, face_id_type>, 2> newpfs;
    newpfs[0] = std::make_pair(p.position(), fid);
    newpfs[1] = std::make_pair(p.position(), fid);

    boost::array<Particle, 2> particles_new;
    particles_new[0] = Particle(sp1, newpfs[0].first, r1, D1);
    particles_new[1] = Particle(sp2, newpfs[1].first, r2, D2);

    bool rejected = false;
    Real separation_factor = r12 * 1e-7;
    std::size_t separation_count = 10;
    while(separation_count != 0)
    {
        --separation_count;
        SGFRD_SCOPE(us, try_to_split, tracer_)

        const Real3 ipv(random_circular_uniform(r12 + separation_factor, fid));
        SGFRD_TRACE(tracer_.write("length of ipv drawn now is %1%", length(ipv)));

        Real3 disp1(ipv * ( r1 / r12)), disp2(ipv * (-r2 / r12));

        if(0 == ecell4::polygon::travel(this->polygon(), newpfs[0], disp1, 100))
        {
            std::cerr << "[WARNING] moving on face by BD: precision lost\n";
        }
        if(0 == ecell4::polygon::travel(this->polygon(), newpfs[1], disp2, 100))
        {
            std::cerr << "[WARNING] moving on face by BD: precision lost\n";
        }

        // if two particle overlaps...
        const Real dist =
            ecell4::polygon::distance(this->polygon(), newpfs[0], newpfs[1]);
        if(dist <= r12)
        {
            separation_factor *= 2.0;
            continue;
        }

        // check whether new particles are inside of the shell, reject the move.
        {
            particles_new[0].position() = newpfs[0].first;
            inside_checker
                is_inside_of(newpfs[0].first, r1, newpfs[0].second, this->polygon());
            if(!boost::apply_visitor(is_inside_of, sh))
            {
                const bool no_overlap = this->burst_and_shrink_overlaps(
                    particles_new[0], newpfs[0].second, did);
                if(!no_overlap)
                {
                    rejected = true;
                    break;
                }
            }
        }
        {
            particles_new[1].position() = newpfs[1].first;
            inside_checker
                is_inside_of(newpfs[1].first, r2, newpfs[1].second, this->polygon());
            if(!boost::apply_visitor(is_inside_of, sh))
            {
                const bool no_overlap = this->burst_and_shrink_overlaps(
                    particles_new[1], newpfs[1].second, did);
                if(!no_overlap)
                {
                    rejected = true;
                    break;
                }
            }
        }
        break;
    }
    if(rejected)
    {
        SGFRD_TRACE(tracer_.write("reaction is rejected because there are no space."))
        return boost::container::static_vector<pid_p_fid_tuple_type, 2>(
                1, boost::make_tuple(pid, p, fid));
    }

    SGFRD_TRACE(tracer_.write("single reaction 1(%1%) -> 2 occurs", pid))

    //------------------------ update particles -----------------------------
    this->update_particle(pid, particles_new[0], newpfs[0].second);
    std::pair<std::pair<ParticleID, Particle>, bool>
        pp2 = this->create_particle(particles_new[1], newpfs[1].second);
    assert(pp2.second);
    const ParticleID pid2(pp2.first.first);

    //------------------------ record reaction -----------------------------
    last_reactions_.push_back(std::make_pair(rule, make_unbinding_reaction_info(
        this->time(), pid, p, pid, particles_new[0], pid2, particles_new[1])));

    boost::container::static_vector<pid_p_fid_tuple_type, 2> retval(2);
    retval[0] = boost::make_tuple(pid , particles_new[0], newpfs[0].second);
    retval[1] = boost::make_tuple(pid2, particles_new[1], newpfs[1].second);
    return retval;
}

boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::escape_single_circular(
        const circular_shell_type& sh, const Single& dom)
{
    SGFRD_SCOPE(us, escape_single_circular, tracer_);

    if(sh.size() == dom.particle().radius())
    {
        SGFRD_TRACE(tracer_.write("closely fitted shell. didnot move."));
        return boost::make_tuple(dom.particle_id(), dom.particle(),
                                 this->get_face_id(dom.particle_id()));
    }

    Particle   p   = dom.particle();
    ParticleID pid = dom.particle_id();

    const Real r   = sh.size() - p.radius();
    const Real theta = this->uniform_real() * 2.0 * M_PI;

    SGFRD_TRACE(tracer_.write("r = %1%, theta = %2%", r, theta))

    const FaceID         fid  = this->get_face_id(pid);
    const triangle_type& face = this->polygon().triangle_at(fid);
    const Real3 direction = rotate(theta, face.normal(), face.represent());

    SGFRD_TRACE(tracer_.write("dir = %1%, len = %2%", direction, length(direction)))

    std::pair<std::pair<Real3, FaceID>, Real3> state =
        std::make_pair(/*position = */std::make_pair(p.position(), fid),
                   /*displacement = */direction * r / length(direction));

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%", state.first.first, state.first.second))

    const std::size_t continue_count =
        ecell4::polygon::travel(this->polygon(), state.first, state.second, 2);
    if(continue_count == 0)
    {
        SGFRD_TRACE(tracer_.write("moving on face: precision lost"))
    }
    SGFRD_TRACE(tracer_.write("escaped"))
    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%",
                              state.first.first, state.first.second))

    p.position() = state.first.first;
    this->update_particle(pid, p, state.first.second);
    SGFRD_TRACE(tracer_.write("particle updated"))
    return boost::make_tuple(pid, p, state.first.second);
}

boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::escape_single_conical(
        const conical_surface_shell_type& sh, const Single& dom)
{
    SGFRD_SCOPE(us, escape_single_conical, tracer_);

    Particle           p   = dom.particle();
    const ParticleID   pid = dom.particle_id();
    const FaceID       fid = this->get_face_id(pid);

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%", p.position(), fid))

    const Real r     = sh.size() - p.radius();
    greens_functions::GreensFunction2DRefWedgeAbs
        gf(p.D(), length(p.position() - sh.position()),
           r,     sh.shape().apex_angle());
    const Real theta = gf.drawTheta(this->uniform_real(), r, dom.dt());

    SGFRD_TRACE(tracer_.write("r = %1%, theta = %2%", r, theta))

    const std::pair<Real3, FaceID> state =
        ecell4::polygon::roll(this->polygon(), std::make_pair(p.position(), fid),
                      sh.structure_id(), r, theta);
    SGFRD_TRACE(tracer_.write("escaped. pos = %1%, fid = %2%", p.position(), fid));

    p.position() = state.first;
    this->update_particle(pid, p, state.second);
    SGFRD_TRACE(tracer_.write("particle updated"))
    return boost::make_tuple(pid, p, state.second);
}

SGFRDSimulator::bursted_type
SGFRDSimulator::burst_single(const Single& dom, const Real tm)
{
    SGFRD_SCOPE(us, burst_single, tracer_);
    const ShellID sid(dom.shell_id());
    SGFRD_TRACE(tracer_.write("shell id = %1%", sid));
    bursted_type results;
    results.push_back(this->propagate_single(this->get_shell(sid), dom, tm));
    this->remove_shell(sid);
    SGFRD_TRACE(tracer_.write("shell removed"));
    return results;
}

void SGFRDSimulator::fire_single(const Single& dom, DomainID did)
{
    SGFRD_SCOPE(us, fire_single, tracer_);
    SGFRD_TRACE(tracer_.write("fire single domain %1%", did))

    const ShellID sid(dom.shell_id());
    switch(dom.eventkind())
    {
        case Single::ESCAPE:
        {
            SGFRD_SCOPE(us, single_escape, tracer_);

            ParticleID pid; Particle p; FaceID fid;
            boost::tie(pid, p, fid) =
                this->escape_single(this->get_shell(sid), dom);

            this->remove_shell(sid);
            SGFRD_TRACE(tracer_.write("shell %1% removed", sid))

            SGFRD_TRACE(tracer_.write("adding next event for %1%", pid))
            this->create_event(pid, p, fid);
            return;
        }
        case Single::REACTION:
        {
            SGFRD_SCOPE(us, single_reaction, tracer_);
            BOOST_AUTO(results,
                       this->reaction_single(this->get_shell(sid), dom, did));
            this->remove_shell(sid);

            ParticleID pid; Particle p; FaceID fid;
            BOOST_FOREACH(boost::tie(pid, p, fid), results)
            {
                SGFRD_TRACE(tracer_.write("adding next event for %1%", pid))
                // here, by calling create_event, the first domain might include
                // the other particle that is one of the result of 1->2 reaction
                add_event(create_closely_fitted_domain(
                            create_closely_fitted_shell(pid, p, fid), pid, p));
            }
            return;
        }
        case Single::UNKNOWN:
        {
            throw std::logic_error("when firing Single: event unspecified");
        }
        default:
        {
            throw std::logic_error("when firing Single: invalid enum value");
        }
    }// switch
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
        }
        remove_event(did_);
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
                clsh.domain_id() = formed_multi_id;
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

    SGFRD_TRACE(tracer_.write("add domain to multi %1% ", multi_to_join_id));

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
            if(did == multi_to_join_id){continue;}

            SGFRD_TRACE(tracer_.write("bursting domain(%1%)", did));
            BOOST_AUTO(ev, get_event(did));

            if(ev->which_domain() == event_type::multi_domain)
            {
                SGFRD_TRACE(tracer_.write("intruder is multi. merge."))
                merge_multi(boost::get<Multi>(ev->domain()), multi_to_join);
                multi_enlarged = true;
            }
            else
            {
                SGFRD_TRACE(tracer_.write("intruder is not multi. burst."))
                ParticleID pid; Particle p; FaceID fid;
                BOOST_FOREACH(boost::tie(pid, p, fid),
                              burst_event(std::make_pair(did, ev), tm))
                {
                    const Real dist = this->polygon().distance(
                            sh_pos, std::make_pair(p.position(), fid)
                            ) - sh.size() - p.radius();
                    const Real min_shell_radius =
                        calc_min_single_circular_shell_radius(p);
                    if(dist < min_shell_radius)
                    {
                        SGFRD_TRACE(tracer_.write("add the particle to multi"))

                        // assign new shell to shell_container and return its ID.
                        BOOST_AUTO(minsh, create_minimum_single_shell(pid, p, fid));

                        multi_to_join.add_particle(pid);
                        multi_to_join.add_shell(minsh.first);

                        // In the next loop, next shell may find this shell.
                        // and if so, the domain_id would not be initialized.
                        mut_sh_vis_applier(didset, multi_to_join);
//                         boost::get<circular_shell_type>(this->get_shell(minsh.first)).domain_id() = multi_to_join_id;

                        multi_enlarged = true;
                    }
                    else // enough distant. add closely-fitted shell
                    {
                        SGFRD_TRACE(tracer_.write("add tight shell to the particle"))
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

        Real dist_to_max_shell_intruder = max_cone_size;
        std::vector<std::pair<DomainID, Real> > min_shell_intruder;
        for(typename std::vector<std::pair<DomainID, Real> >::const_iterator
                iter = intrusive_domains.begin(), end = intrusive_domains.end();
                iter != end; ++iter)
        {
            if(iter->second <= min_cone_size)
            {
                min_shell_intruder.push_back(*iter);
            }
            else
            {
                // XXX because `intrusive_domains` are sorted by comparing the
                // distance to them, once we found the element is far away, all
                // the successors are much further.
                dist_to_max_shell_intruder =
                    std::min(iter->second, dist_to_max_shell_intruder);
                break;
            }
        }

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

            const Real shell_size = dist_to_max_shell_intruder *
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
                std::min(dist_to_max_shell_intruder,
                         shrinked_or_multi.front().second) *
                single_conical_surface_shell_mergin;

            /* assertion */{
            if(false == this->shell_container_.list_shells_within_radius(
                std::make_pair(this->polygon().vertex_at(vid).position, vid),
                shell_size).empty())
            {
                std::cout << "nearest shell: "
                    << this->shell_container_.list_shells_within_radius(
                        std::make_pair(this->polygon().vertex_at(vid).position, vid),
                        shell_size).front().first.first << std::endl;
                std::cout << "distance: "
                    << this->shell_container_.list_shells_within_radius(
                        std::make_pair(this->polygon().vertex_at(vid).position, vid),
                        shell_size).front().second << std::endl;
                std::cout << "shell size " << shell_size << std::endl;
                assert(false);
            }
            }

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

    Real dist_to_max_shell_intruder = max_circle_size;
    std::vector<std::pair<DomainID, Real> > min_shell_intruder;
    for(typename std::vector<std::pair<DomainID, Real> >::const_iterator
            iter = intrusive_domains.begin(), end = intrusive_domains.end();
            iter != end; ++iter)
    {
        SGFRD_TRACE(tracer_.write("check domain %1%: distance = %2%",
                    iter->first, iter->second));
        if(iter->second <= min_circle_size)
        {
            SGFRD_TRACE(tracer_.write("%1% is inside of minimum circle size",
                        iter->first));
            min_shell_intruder.push_back(*iter);
        }
        else
        {
            SGFRD_TRACE(tracer_.write(
                "%1% does not intersect with minimum circle", iter->first));

            // XXX because `intrusive_domains` are sorted by comparing the
            // distance to them, once we found the element is far away, all
            // the successors are much further.
            dist_to_max_shell_intruder =
                std::min(iter->second, dist_to_max_shell_intruder);
            break;
        }
    }

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
            dist_to_max_shell_intruder * single_circular_shell_mergin;

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
            std::min(dist_to_max_shell_intruder, shrinked_or_multi.front().second) *
            single_circular_shell_mergin;

        return add_event(create_single(
                    create_single_circular_shell(pos, shell_size), pid, p));
    }

    SGFRD_TRACE(tracer_.write("forming multi"))
    return form_multi(pid, p, fid, shrinked_or_multi);
}

} // sgfrd
} // ecell4
