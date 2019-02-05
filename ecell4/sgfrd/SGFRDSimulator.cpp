#include "SGFRDSimulator.hpp"

namespace ecell4
{
namespace sgfrd
{
const Real SGFRDSimulator::single_circular_shell_factor        = 1.1;
const Real SGFRDSimulator::single_circular_shell_mergin        = 1.0 - 1e-7;
const Real SGFRDSimulator::single_conical_surface_shell_factor = 1.1;
const Real SGFRDSimulator::single_conical_surface_shell_mergin = 1.0 - 1e-7;
const Real SGFRDSimulator::minimum_separation_factor           = 1e-7;

boost::tuple<ParticleID, Particle, SGFRDSimulator::FaceID>
SGFRDSimulator::propagate_single_circular(
        const circular_shell_type& sh, const Single& dom, const Real tm)
{
    SGFRD_SCOPE(us, propagate_single_circular, tracer_);

    Particle     p   = dom.particle();
    ParticleID   pid = dom.particle_id();
    const FaceID fid = this->get_face_id(pid);

    if(p.D() == 0.0)
    {
        // it never moves.
        return boost::make_tuple(pid, p, fid);
    }

    greens_functions::GreensFunction2DAbsSym gf(p.D(), sh.size() - p.radius());

    const Real del_t = tm - dom.begin_time();

    SGFRD_TRACE(tracer_.write("delta t for domain having shell %1% is %2%",
                del_t, dom.shell_id()));
    SGFRD_TRACE(tracer_.write("its own dt = %1%, and begin_time = %2%",
                dom.dt(), dom.begin_time()));

    const Real r     = gf.drawR(this->uniform_real(), del_t);
    const Real theta = this->uniform_real() * 2 * M_PI;

    SGFRD_TRACE(tracer_.write("r = %1%, theta = %2%", r, theta));

    const Triangle& face = this->polygon().triangle_at(fid);
    const Real3 direction = rotate(theta, face.normal(), face.represent());
    const Real  len_direction = length(direction);

    SGFRD_TRACE(tracer_.write("direction = %1%, len = %2%", direction, len_direction));

    std::pair<std::pair<Real3, FaceID>, Real3> state =
        std::make_pair(std::make_pair(p.position(), fid),
                       direction * r / len_direction);

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%", state.first.first, state.first.second));

    state.first = ecell4::polygon::travel(this->polygon(), state.first, state.second, 2);

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%, dsp = %3%",
                state.first.first, state.first.second, state.second))

    p.position() = state.first.first;
    this->update_particle(pid, p, state.first.second);
    SGFRD_TRACE(tracer_.write("particle updated"))

    assert(// check particle is still inside of the shell after this propagation
        p.radius() + ecell4::polygon::distance(this->polygon(),
                state.first, std::make_pair(sh.position(), sh.structure_id()))
        <= sh.size() * (1.0+1e-8) /* <- tolerance */);

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

    if(p.D() == 0.0)
    {
        // it never moves.
        return boost::make_tuple(pid, p, fid);
    }

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%", p.position(), fid));

    const Real D   = p.D();
    const Real a   = sh.size() - p.radius();
    const Real phi = sh.shape().apex_angle();
    const Real r0  = length(
        this->polygon().periodic_transpose(p.position(), sh.shape().apex()) -
        sh.shape().apex());

    greens_functions::GreensFunction2DRefWedgeAbs gf(D, r0, a, phi);

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

    assert(// check particle is still inside of the shell after this propagation
        length(polygon().periodic_transpose(p.position(), sh.shape().apex()) -
               sh.shape().apex()) + p.radius() <= sh.size() * (1.0+1e-8)
    );

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
    const Real radius_new = species_new.get_attribute_as<Real>("radius");
    const Real D_new      = species_new.get_attribute_as<Real>("D");
    const Particle p_new(species_new, p.position(), radius_new, D_new);

    inside_checker is_inside_of(
            p_new.position(), p_new.radius(), fid, this->polygon());
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

    const Real D1 = sp1.get_attribute_as<Real>("D");
    const Real D2 = sp2.get_attribute_as<Real>("D");
//     const Real D12 = D1 + D2
    const Real r1 = sp1.get_attribute_as<Real>("radius");
    const Real r2 = sp2.get_attribute_as<Real>("radius");
    const Real r12 = r1 + r2;

    SGFRD_TRACE(tracer_.write("products has D1(%1%), D2(%2%), r1(%3%), r2(%4%)",
                D1, D2, r1, r2));

    boost::array<std::pair<Real3, FaceID>, 2> newpfs;
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
        SGFRD_SCOPE(us, try_to_split, tracer_);

        SGFRD_TRACE(tracer_.write("separation count = %1%", separation_count));

        const Real3 ipv(random_circular_uniform(r12 + separation_factor, fid));
        SGFRD_TRACE(tracer_.write("length of ipv drawn now is %1%", length(ipv)));

        Real3 disp1(ipv * ( r1 / r12)), disp2(ipv * (-r2 / r12));

        newpfs[0] = ecell4::polygon::travel(this->polygon(), newpfs[0], disp1);
        newpfs[1] = ecell4::polygon::travel(this->polygon(), newpfs[1], disp2);

        // if two particle overlaps...
        const Real dist =
            ecell4::polygon::distance(this->polygon(), newpfs[0], newpfs[1]);
        if(dist <= r12)
        {
            newpfs[0] = std::make_pair(p.position(), fid); // rollback
            newpfs[1] = std::make_pair(p.position(), fid);
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
    SGFRD_TRACE(tracer_.write("ext particle %1% is updateed", pid))
    SGFRD_TRACE(tracer_.write("new particle %1% is assigned", pid2))

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

    const FaceID    fid  = this->get_face_id(pid);
    const Triangle& face = this->polygon().triangle_at(fid);
    const Real3 direction = rotate(theta, face.normal(), face.represent());

    SGFRD_TRACE(tracer_.write("dir = %1%, len = %2%", direction, length(direction)))

    const Real3 displacement = direction * r / length(direction);
    SGFRD_TRACE(tracer_.write("displacement = %1%", displacement))

    std::pair<std::pair<Real3, FaceID>, Real3> state =
        std::make_pair(std::make_pair(p.position(), fid), displacement);

    SGFRD_TRACE(tracer_.write("pos = %1%, fid = %2%", state.first.first, state.first.second))

    state.first = ecell4::polygon::travel(this->polygon(), state.first, state.second, 2);
    SGFRD_TRACE(tracer_.write("escaped. pos = %1%, fid = %2%",
                              state.first.first, state.first.second))

    p.position() = state.first.first;
    this->update_particle(pid, p, state.first.second);
    SGFRD_TRACE(tracer_.write("particle updated"))

    assert(// check particle is still inside of the shell after this propagation
        p.radius() + ecell4::polygon::distance(this->polygon(),
                state.first, std::make_pair(sh.position(), sh.structure_id()))
        <= sh.size() * (1.0+1e-8) /* <- tolerance */);

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
    SGFRD_TRACE(tracer_.write("escaped. pos = %1%, fid = %2%",
                              state.first, state.second));

    p.position() = state.first;
    this->update_particle(pid, p, state.second);
    SGFRD_TRACE(tracer_.write("particle updated"))

    assert(// check particle is still inside of the shell after this propagation
        length(polygon().periodic_transpose(p.position(), sh.shape().apex()) -
               sh.shape().apex()) + p.radius() <= sh.size() * (1.0+1e-8)
    );

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
    STAT(if(this->get_shell(dom.shell_id()).which() ==
            shell_container_type::circular_shell){
            stat_fired_events.add_count(FireSingleCircular);
        } else {
            stat_fired_events.add_count(FireSingleConical);
        })

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
                add_event(create_tight_domain(
                            create_tight_shell(pid, p, fid), pid, p));
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
        const Particle& p, const FaceID& fid, const DomainID& did)
{
    // Here, particle is being outside of a shell, by reaction or BDstep.
    // But before bursting, the particle is not assigned to World.
    // So it has no ID.
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
        if(!(this->event_exists(did_)))
        {
            SGFRD_TRACE(tracer_.write("domain %1% does not exist. it may have "
                        "already been fired.", did_));
            assert(false);
        }

        boost::shared_ptr<event_type> ev_(get_event(did_));

        if(ev_->which_domain() == SGFRDEvent::multi_domain)
        {
            // Multi overlaps with Multi. it sometimes causes recursive bursting.
            // in this case, it is considered as a collision.
            no_overlap = false;
            continue;
        }

        ParticleID pid_; Particle p_; FaceID fid_;
        BOOST_FOREACH(boost::tie(pid_, p_, fid_),
                      burst_event(std::make_pair(did_, ev_), tm))
        {
            const Real dist = ecell4::polygon::distance(this->polygon(),
                std::make_pair(p.position(), fid),
                std::make_pair(p_.position(), fid_));
            no_overlap = no_overlap && (dist > p.radius() + p_.radius());
            add_event(create_tight_domain(
                create_tight_shell(pid_, p_, fid_), pid_, p_));
        }
        remove_event(did_);
    }
    return no_overlap;
}

boost::optional<DomainID>
SGFRDSimulator::form_pair(
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const std::vector<std::pair<DomainID, Real> >& intruders)
{
    SGFRD_SCOPE(us, form_pair, tracer_);

    // the first (nearest) domain in the intruders is the partner to form pair.
    const boost::shared_ptr<event_type> nearest =
        this->get_event(intruders.front().first);
    if(nearest->which_domain() != event_type::single_domain)
    {
        // if the partner is Pair of Multi, pair cannot be formed.
        SGFRD_TRACE(tracer_.write(
                    "nearest intruder is not single. can't form pair."))
        return boost::none;
    }

    // the domain that will be consumed to form a pair.
    const Single& sgl = boost::get<Single>(nearest->domain());
    const ParticleID partner_id  = sgl.particle_id();
    const Particle   partner     = sgl.particle();

    if(p.D() == 0.0 && partner.D() == 0.0)
    {
        return boost::none;
    }
    if(sgl.dt() != 0.0)
    {
        SGFRD_TRACE(tracer_.write(
                    "nearest intruder is not tight. first we need to burst it"));
        return boost::none;
    }


    const FaceID   partner_fid = this->get_face_id(partner_id);
    const ShellID  partner_sid = sgl.shell_id();
    const DomainID partner_did = intruders.front().first;

    SGFRD_TRACE(tracer_.write(
                "try to form pair domain for %1% and %2%", pid, partner_id))

    const Real r1(p.radius()), r2(partner.radius());
    const Real D1(p.D()), D2(partner.D());
    const Real D12 = D1 + D2;
    const Real r12 = r1 + r2;

    const Real3 ipv = ecell4::polygon::direction(this->polygon(),
        std::make_pair(p.position(), fid), // ->
        std::make_pair(partner.position(), partner_fid));
    const Real len_ipv = length(ipv);
    const Real sh_minim =
        std::max(len_ipv * D1 / D12 + r1, len_ipv * D2 / D12 + r2) * 3;
    // XXX this `3` is just a parameter. it should be tuned later.

    if(len_ipv < r12)
    {
        throw std::runtime_error((boost::format(
            "form_pair: particle %1% and %2% already collides!\n"
            "len_ipv  = %3%\n"
            "r1       = %4%\n"
            "r2       = %5%\n"
            "r12      = %6%\n") % pid % partner_id % len_ipv % r1 % r2 % r12
            ).str());
    }

    const std::pair<Real3, FaceID> pos_com = ecell4::polygon::travel(
        this->polygon(), std::make_pair(p.position(), fid), ipv * D1 / D12, 2);

    // check other shells are not in the range...
    Real max_dist = get_max_circle_size(pos_com);
    for(std::vector<std::pair<DomainID, Real> >::const_iterator
        iter(intruders.begin()+1), iend(intruders.end()); iter != iend; ++iter)
    {
        const boost::shared_ptr<event_type> intruder_ev =
            this->get_event(iter->first);
        if(intruder_ev->which_domain() != event_type::single_domain)
        {
            continue;
        }
        const Single&    intruder_dom = boost::get<Single>(nearest->domain());
        const ParticleID intruder_pid = intruder_dom.particle_id();
        const Particle&  intruder_p   = intruder_dom.particle();
        const FaceID     intruder_fid = this->get_face_id(intruder_pid);

        const Real d_to_sh = ecell4::polygon::distance(this->polygon(),
            pos_com, std::make_pair(intruder_p.position(), intruder_fid)) -
            calc_min_single_circular_shell_radius(intruder_p);

        if(d_to_sh < sh_minim)
        {
            SGFRD_TRACE(tracer_.write(
                        "intrusive domains exists. multi should be formed"))
            // other shells overlap to the pair. multi should be formed.
            return boost::none;
        }
        max_dist = std::min(max_dist, d_to_sh);
    }

    std::vector<std::pair<std::pair<ShellID, shell_type>, Real> >
        other_shells(this->shell_container_.list_shells_within_radius(
                     pos_com, max_dist));
    Real pair_shell_size = max_dist;
    for(std::vector<std::pair<std::pair<ShellID, shell_type>, Real> >::iterator
            iter = other_shells.begin(), iend = other_shells.end();
            iter != iend; ++iter)
    {
        if(iter->first.first == partner_sid)
        {
            continue;
        }
        pair_shell_size = std::min(pair_shell_size, iter->second);
    }

    if(pair_shell_size >= sh_minim)
    {
        SGFRD_TRACE(tracer_.write("pair shell size is larger than the minimum. "
                                  "pair can be formed"))

        this->remove_shell(partner_sid);
        SGFRD_TRACE(tracer_.write("remove partner's shell, %1%", partner_sid))

        this->remove_event(partner_did);
        SGFRD_TRACE(tracer_.write("remove partner's domain, %1%", partner_did))

        // check that the pair shell does not overlap with other particles
        assert(get_intrusive_domains(pos_com, pair_shell_size).empty());

        const ShellID shid(shell_id_gen());
        const circle_type pair_circle(
                pair_shell_size * single_circular_shell_mergin, pos_com.first,
                this->polygon().triangle_at(pos_com.second).normal());
        const circular_shell_type pair_shell(pair_circle, pos_com.second);
        shell_container_.check_add_shell(shid, pair_shell, pos_com.second,
                                         "create pair shell");

        return add_event(create_pair(
                    std::make_pair(shid, pair_circle),
                    pid, p, partner_id, partner, ipv, len_ipv));
    }
    SGFRD_TRACE(tracer_.write("pair shell size = %1%, min_shell_size = %2%",
                pair_shell_size, sh_minim);)
    SGFRD_TRACE(tracer_.write("min-pair-intruder exists. multi should be formed"))
    return boost::none;
}

DomainID SGFRDSimulator::form_multi(
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const std::vector<std::pair<DomainID, Real> >& doms)
{
    SGFRD_SCOPE(us, form_multi, tracer_);

    Multi new_multi(create_empty_multi());
    const DomainID formed_multi_id = this->add_event(new_multi);
    Multi& formed_multi =
        boost::get<Multi>(get_event(formed_multi_id)->domain());
    SGFRD_TRACE(tracer_.write("new multi(%1%) created", formed_multi_id))

    const domain_id_setter didset(formed_multi_id);

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
                this->update_shell(sid, clsh, clsh.structure_id());
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
    formed_multi.determine_reaction_length();
    formed_multi.determine_delta_t();

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
                    const Real dist = ecell4::polygon::distance(this->polygon(),
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
                        add_event(create_tight_domain(
                            create_tight_shell(pid, p, this->get_face_id(pid)),
                            pid, p));
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

expected<DomainID, std::vector<std::pair<DomainID, Real> > >
SGFRDSimulator::form_single_conical_event(
        const ParticleID& pid, const Particle& p, const FaceID fid)
{
    SGFRD_SCOPE(us, form_single_conical_event, tracer_);

    const std::pair<Real3, FaceID> pos = std::make_pair(p.position(), fid);
    const std::vector<std::pair<VertexID, Real> > intrusive_vertices(
            get_intrusive_vertices(pos, std::numeric_limits<Real>::infinity()));

    const VertexID& vid  = intrusive_vertices.front().first;
    const Real dist_to_v = intrusive_vertices.front().second;
    SGFRD_TRACE(tracer_.write("vertex id = %1%, distance = %2%", vid, dist_to_v));

    if(p.D() == 0.0)
    {
        SGFRD_TRACE(tracer_.write("diffusion coefficient is 0."))
        SGFRD_TRACE(tracer_.write("creating tight conical single event"))

        const Real shell_size = (p.radius() + dist_to_v) *
                                (1.0 + minimum_separation_factor);

        // check particle does not overlap with any others.
        assert(get_intrusive_domains(vid, shell_size).empty());

        return ok(add_event(create_single(
            create_single_conical_surface_shell(vid, shell_size), pid, p)));
    }

    const Real min_cone_size = (p.radius() + dist_to_v) *
                               single_conical_surface_shell_factor;
    const Real max_cone_size = get_max_cone_size(vid);
    SGFRD_TRACE(tracer_.write("min_cone_size = %1%, max_cone_size = %2%",
                min_cone_size, max_cone_size));

    if(min_cone_size > max_cone_size)
    {
        // cannot form cone nor circle. use multi.
        return err(std::vector<std::pair<DomainID, Real> >(0));
    }

    const std::vector<std::pair<DomainID, Real> > intrusive_domains(
            get_intrusive_domains(vid, max_cone_size));
    SGFRD_TRACE(tracer_.write("intrusive_domain_size = %1%", intrusive_domains.size()));

    if(intrusive_domains.empty())
    {
        return ok(add_event(create_single(
            create_single_conical_surface_shell(vid, max_cone_size), pid, p)));
    }

    Real dist_to_max_shell_intruder = max_cone_size;
    std::vector<std::pair<DomainID, Real> > min_shell_intruder;
    for(std::vector<std::pair<DomainID, Real> >::const_iterator
        iter = intrusive_domains.begin(), iend = intrusive_domains.end();
        iter != iend; ++iter)
    {
        // here, it calculates the distance between domain edges and the vertex.
        // the value `iter->second` is not the distance between particles,
        // it is a distance between domain and vertex.
        //
        // So it is not a problem that the value `iter->second` is shorter than
        // the raidus of the particle.

        if(iter->second <= min_cone_size)
        {
            SGFRD_TRACE(tracer_.write("%1% <= %2%, min shell intruder found!",
                        iter->second, min_cone_size));
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
        SGFRD_TRACE(
            tracer_.write("intrusive domains exist but enough distant");
            for(std::size_t i=0; i<intrusive_domains.size(); ++i)
            {
                tracer_.write("domain %1%; dist = %2%;",
                    intrusive_domains[i].first, intrusive_domains[i].second);
            }
        )

        const Real shell_size = dist_to_max_shell_intruder *
                                single_conical_surface_shell_mergin;

        return ok(add_event(create_single(
            create_single_conical_surface_shell(vid, shell_size), pid, p)));
    }

    // burst intruder_domains and get new positions of particles
    std::vector<std::pair<DomainID, Real> > shrinked_or_multi =
        burst_and_shrink_non_multis(vid, min_shell_intruder);
    SGFRD_TRACE(tracer_.write("close domains are bursted."));

    if(shrinked_or_multi.front().second > min_cone_size)
    {
        SGFRD_TRACE(
            tracer_.write("after burst, no intruder exist in the min-range %1%",
                          min_cone_size);
            for(std::size_t i=0; i<shrinked_or_multi.size(); ++i)
            {
                tracer_.write("domain %1%; dist = %2%;",
                    shrinked_or_multi[i].first, shrinked_or_multi[i].second);
            }
        )
        const Real shell_size =
            std::min(dist_to_max_shell_intruder,
                     shrinked_or_multi.front().second) *
            single_conical_surface_shell_mergin;

        /* assertion */{
        if(false == this->shell_container_.list_shells_within_radius(
            std::make_pair(this->polygon().position_at(vid), vid),
            shell_size).empty())
        {
            std::cout << "nearest shell: "
                << this->shell_container_.list_shells_within_radius(
                    std::make_pair(this->polygon().position_at(vid), vid),
                    shell_size).front().first.first << std::endl;
            std::cout << "distance: "
                << this->shell_container_.list_shells_within_radius(
                    std::make_pair(this->polygon().position_at(vid), vid),
                    shell_size).front().second << std::endl;
            std::cout << "shell size " << shell_size << std::endl;
            assert(false);
        }
        }

        return ok(add_event(create_single(
            create_single_conical_surface_shell(vid, shell_size), pid, p)));
    }

    return err(shrinked_or_multi);
}

expected<DomainID, std::vector<std::pair<DomainID, Real> > >
SGFRDSimulator::form_single_circular_event(
    const ParticleID& pid, const Particle& p, const FaceID fid,
    const Real max_circle_size)
{
    SGFRD_SCOPE(us, form_single_circular_event, tracer_);
    SGFRD_TRACE(tracer_.write("forming single domain for particle %1% r = %2%",
                pid, p.radius()));

    const Real min_circle_size = p.radius() * single_circular_shell_factor;
    const std::pair<Real3, FaceID> pos = std::make_pair(p.position(), fid);

    if(p.D() == 0.0)
    {
        SGFRD_TRACE(tracer_.write("diffusion coefficient is 0."))
        SGFRD_TRACE(tracer_.write("creating tight single event"))
        const Real shell_size = p.radius() * (1.0 + minimum_separation_factor);

        // check the particle does not overlap with any others.
        assert(get_intrusive_domains(pos, shell_size).empty());
        return ok(add_event(create_single(
                create_single_circular_shell(pos, shell_size), pid, p)));
    }

    /* XXX:TAKE CARE! the distance in the element of intrusive_domains, typed *
     * as `std::pair<DomainID, Real>` is not a distance between particle and  *
     * shell, but a distance between center of particle and shell surface.    */
    const std::vector<std::pair<DomainID, Real> > intrusive_domains(
            get_intrusive_domains(pos, max_circle_size));
    SGFRD_TRACE(tracer_.write(
                "intrusive_domain_size = %1%", intrusive_domains.size()))

    if(intrusive_domains.empty())
    {
        SGFRD_TRACE(tracer_.write("no intrusive domains exists."))
        SGFRD_TRACE(tracer_.write(
            "creating single event; shell size = %1%", max_circle_size))
        return ok(add_event(create_single(
            create_single_circular_shell(pos, max_circle_size), pid, p)));
    }

    Real distance_to_nearest = max_circle_size; //XXX nearest (but not intruder)
    std::vector<std::pair<DomainID, Real> > min_shell_intruder;
    for(std::vector<std::pair<DomainID, Real> >::const_iterator
            iter = intrusive_domains.begin(), end = intrusive_domains.end();
            iter != end; ++iter)
    {
        SGFRD_TRACE(tracer_.write("check domain %1%: distance = %2%",
                    iter->first, iter->second));
        if(iter->second <= min_circle_size)
        {
            SGFRD_TRACE(tracer_.write("%1% is inside of minimum circle size",
                        iter->first));
            if(iter->second < p.radius())
            {
                throw std::runtime_error((
                    boost::format("form_single_circular_event: nearest domain "
                        "%1% overlaps with particle %2%. distance from point = "
                        "%3%, radius = %4%") % iter->first % pid %
                        iter->second % p.radius()
                    ).str());
            }
            min_shell_intruder.push_back(*iter);
            // collect all the min-shell-intruders.
            continue;
        }
        else
        {
            SGFRD_TRACE(tracer_.write(
                "%1% does not intersect with minimum circle", iter->first));

            // calculate modest distance if this one is a single domain.
            boost::shared_ptr<event_type> ev(this->get_event(iter->first));
            if(ev->which_domain() == event_type::single_domain)
            {
                SGFRD_TRACE(tracer_.write("calculating modest r."))
                const Single&     sgl       = boost::get<Single>(ev->domain());
                const Particle&   nearest_p = sgl.particle();
                const shell_type& sh        = this->get_shell(sgl.shell_id());
                const Real sh_size = boost::apply_visitor(shell_size_getter(), sh);

                SGFRD_TRACE(tracer_.write("raw distance = %1%", iter->second))
                SGFRD_TRACE(tracer_.write("shell radius = %1%", sh_size))
                SGFRD_TRACE(tracer_.write("nearp radius = %1%", nearest_p.radius()))

                const Real modest_dist = calc_modest_shell_size(p, nearest_p,
                    iter->second + sh_size - p.radius() - nearest_p.radius());
                distance_to_nearest = std::min(iter->second, modest_dist);
            }
            else // nearest domain is not a single domain.
            {
                distance_to_nearest = iter->second;
            }
            SGFRD_TRACE(tracer_.write("distance_to_nearest = %1%",
                                      distance_to_nearest));

            // XXX because `intrusive_domains` are sorted by their distance,
            //     from nearest to distant, all the rests are far away.
            break;
        }
    }

    if(min_shell_intruder.empty())
    {
        SGFRD_TRACE(
            tracer_.write("intrusive domains exists but enough distant");
            for(std::size_t i=0; i<intrusive_domains.size(); ++i)
            {
                tracer_.write("domain %1%; dist = %2%;",
                    intrusive_domains[i].first, intrusive_domains[i].second);
            }
        )
        const Real shell_size =
            distance_to_nearest * single_circular_shell_mergin;

        assert(get_intrusive_domains(pos, shell_size).empty());

        SGFRD_TRACE(tracer_.write(
            "creating single event; shell size = %1%", shell_size))

        return ok(add_event(create_single(
            create_single_circular_shell(pos, shell_size), pid, p)));
    }

    std::vector<std::pair<DomainID, Real> > shrinked_or_multi =
        burst_and_shrink_non_multis(pid, p, fid, min_shell_intruder);
    SGFRD_TRACE(tracer_.write("min_shell_intruder domains are bursted"))

    if(shrinked_or_multi.front().second > min_circle_size)
    {
        SGFRD_TRACE(
            tracer_.write("after burst, no intruders exist");
            for(std::size_t i=0; i<shrinked_or_multi.size(); ++i)
            {
                tracer_.write("domain %1%; dist = %2%;",
                    shrinked_or_multi[i].first, shrinked_or_multi[i].second);
            }
        )

        const Real shell_size =
            std::min(distance_to_nearest, shrinked_or_multi.front().second) *
            single_circular_shell_mergin;

        assert(get_intrusive_domains(pos, shell_size).empty());

        SGFRD_TRACE(tracer_.write(
            "creating single event; shell size = %1%", shell_size))

        return ok(add_event(create_single(
            create_single_circular_shell(pos, shell_size), pid, p)));
    }

    SGFRD_TRACE(tracer_.write("failed to create single."))
    return err(shrinked_or_multi);
}

DomainID SGFRDSimulator::create_event(
            const ParticleID& pid, const Particle& p, const FaceID fid)
{
    SGFRD_SCOPE(us, create_event, tracer_);
    const std::pair<Real3, FaceID> pos = std::make_pair(p.position(), fid);

    const Real min_circle_size = p.radius() * single_circular_shell_factor;
    const Real max_circle_size = get_max_circle_size(pos);

    SGFRD_TRACE(tracer_.write("min_circle_size = %1%, max_circle_size = %2%",
                min_circle_size, max_circle_size));

    if(max_circle_size < min_circle_size)// draw conical shell
    {
        expected<DomainID, std::vector<std::pair<DomainID, Real> >
            > single_conical = this->form_single_conical_event(pid, p, fid);
        if(single_conical.is_ok())
        {
            SGFRD_TRACE(tracer_.write("single conical was successfully formed"))
            return single_conical.unwrap();
        }
        else
        {
            STAT(stat_multi_reason.add_count(SingleConicalFailed);)
            SGFRD_TRACE(tracer_.write("single conical could not be formed"))
            SGFRD_TRACE(tracer_.write("forming multi..."))
            return form_multi(pid, p, fid, single_conical.unwrap_error());
        }
    }
    else // draw circluar shell
    {
        expected<DomainID, std::vector<std::pair<DomainID, Real> >
            > single_circular =
                this->form_single_circular_event(pid, p, fid, max_circle_size);
        if(single_circular.is_ok())
        {
            SGFRD_TRACE(tracer_.write("single circular was successfully formed"))
            return single_circular.unwrap();
        }

        SGFRD_TRACE(tracer_.write("single circular could not be formed"))

        const std::vector<std::pair<DomainID, Real> >& intruders =
            single_circular.unwrap_error();

        boost::optional<DomainID> pair_ =
            this->form_pair(pid, p, fid, intruders);
        if(pair_)
        {
            SGFRD_TRACE(tracer_.write("pair circular was formed"))
            return *pair_;
        }
        STAT(stat_multi_reason.add_count(PairFailed);)
        SGFRD_TRACE(tracer_.write("forming multi..."))
        return form_multi(pid, p, fid, intruders);
    }
}

bool SGFRDSimulator::diagnosis() const
{//{{{
//     const boost::chrono::steady_clock::time_point start_ =
//         boost::chrono::high_resolution_clock::now();

    bool result = true;
    // 1. check overlap between particles
    // 2. check overlap between shells
    // 3. check all the particles are inside of its shell

    BOOST_AUTO(particles, this->world_->list_particles());
    BOOST_AUTO(shells,    this->shell_container_.list_shells());

    // 1.
    ParticleID pid; Particle p;
    BOOST_FOREACH(boost::tie(pid, p), particles)
    {
        const FaceID fid = this->get_face_id(pid);
        std::pair<Real3, FaceID> pos = std::make_pair(p.position(), fid);

        ParticleID _pid; Particle _p;
        BOOST_FOREACH(boost::tie(_pid, _p), particles)
        {
            if(pid == _pid) {continue;}
            const FaceID _fid = this->get_face_id(_pid);
            const Real dist = this->world_->distance(pos,
                                  std::make_pair(_p.position(), _fid));
            if(dist < p.radius() + _p.radius())
            {
                result = false;
                std::cerr << "ERROR: particle " << pid << " and " << _pid
                          << "overlaps!\n";
                std::cerr << "     : distance = " << dist << " < sum of radii = "
                          << p.radius() + _p.radius() << '\n';
                std::cerr << "     : particle " << pid << " has radius "
                          << p.radius() << " at " << p.position() << " on "
                          << fid << '\n';
                std::cerr << "     : particle " << _pid << " has radius "
                          << _p.radius() << " at " << _p.position() << " on "
                          << _fid << '\n';
            }
        }

        if(!this->polygon().is_inside_of_boundary(p.position()))
        {
            std::cerr << "ERROR: particle " << pid << " is outside of the boundary!\n";
            std::cerr << "     : position = " << p.position() << ", boundary = " << polygon().edge_lengths() << "\n";
            result = false;
        }

        const Triangle&   tri  = this->polygon().triangle_at(fid);
        const Barycentric bary = ::ecell4::to_barycentric(p.position(), tri);
        if(!is_inside(bary))
        {
            std::cerr << "ERROR: particle " << pid << " is not on the face " << fid << "\n";
            std::cerr << "     : position    = " << p.position() << ", face = " << tri << "\n";
            std::cerr << "     : barycentric = " << bary << ", face = " << tri << "\n";
            result = false;
        }
    }

    // 2.
    ShellID shid; shell_type sh;
    BOOST_FOREACH(boost::tie(shid, sh), shells)
    {
        ShellID _shid; shell_type _sh;
        switch(sh.which())
        {
            case shell_container_type::circular_shell:
            {
                circular_shell_type ccl = boost::get<circular_shell_type>(sh);
                distance_calculator_on_surface<FaceID>
                    dist_calc(ccl.get_surface_position(), this->polygon());

                BOOST_FOREACH(boost::tie(_shid, _sh), shells)
                {
                    if(_shid == shid){continue;}
                    const Real dist = boost::apply_visitor(dist_calc, _sh);
                    const DomainID _did = boost::apply_visitor(domain_id_getter(), _sh);
                    if(dist < ccl.size() && ccl.domain_id() != _did)
                    {
                        result = false;
                        std::cerr << "ERROR: circular shell " << shid
                                  << " and shell " << _shid << "overlaps\n";
                        std::cerr << "     : distance = " << dist - ccl.size()
                                  << '\n';
                        std::cerr <<  shid << " -> " <<  sh << '\n';
                        std::cerr << _shid << " -> " << _sh << '\n';
                    }
                }
                break;
            }
            case shell_container_type::conical_shell:
            {
                conical_surface_shell_type con =
                    boost::get<conical_surface_shell_type>(sh);
                distance_calculator_on_surface<VertexID>
                    dist_calc(con.get_surface_position(), this->polygon());

                BOOST_FOREACH(boost::tie(_shid, _sh), shells)
                {
                    if(_shid == shid){continue;}
                    const Real dist = boost::apply_visitor(dist_calc, _sh);
                    const DomainID _did = boost::apply_visitor(
                            domain_id_getter(), _sh);
                    if(dist < con.size() && con.domain_id() != _did)
                    {
                        result = false;
                        std::cerr << "ERROR: conical shell " << shid
                                  << " and shell " << _shid << "overlaps\n";
                        std::cerr << "     : distance = " << dist - con.size()
                                  << "\n";
                        std::cerr <<  shid << " -> " <<  sh << '\n';
                        std::cerr << _shid << " -> " << _sh << '\n';
                    }
                }
                break;
            }
            default:
            {
                result = false;
                std::cerr << "ERROR: shell " << shid
                          << " has invalid which() value " << sh.which() << '\n';
                break;
            }
        }
    }

    // 3.
    std::map<ParticleID, EventID> pid2evid;
    std::map<ShellID,    EventID> sid2evid;
    EventID evid; boost::shared_ptr<event_type> ev_ptr;
    BOOST_FOREACH(boost::tie(evid, ev_ptr), this->scheduler_.events())
    {
        SGFRDEvent::domain_type const& dom = ev_ptr->domain();
        switch(ev_ptr->which_domain())
        {
            case event_type::single_domain:
            {
                const Single&   sgl  = boost::get<Single>(dom);
                const ShellID   _shid = sgl.shell_id();
                const ParticleID _pid = sgl.particle_id();

                if(pid2evid.count(_pid)  == 0){pid2evid[_pid]  = evid;}
                if(sid2evid.count(_shid) == 0){sid2evid[_shid] = evid;}

                BOOST_AUTO(found_p, std::find_if(particles.begin(), particles.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                        ParticleID, Particle>(_pid)));
                if(found_p == particles.end())
                {
                    result = false;
                    std::cerr << "ERROR: particle might assigned to two "
                              << "different domains\n";
                    std::cerr << "     : Single domain " << evid << " has particle "
                              << _pid << " but the particle is already erased\n";
                    if(pid2evid.count(_pid) == 1)
                    {
                        std::cerr << "     : event " << pid2evid[_pid]
                                  << " has particle " << _pid << '\n';
                    }
                    break;
                }
                BOOST_AUTO(found_s, std::find_if(shells.begin(), shells.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                        ShellID, shell_type>(_shid)));
                if(found_s == shells.end())
                {
                    result = false;
                    std::cerr << "ERROR: shell might assigned to two"
                              << "different domains\n";
                    std::cerr << "     : Single domain " << evid << " has shell "
                              << _shid << " but the shell is already erased\n";
                    if(sid2evid.count(_shid) == 1)
                    {
                        std::cerr << "     : event " << sid2evid[_shid]
                                  << " has shell " << _shid << '\n';
                    }

                    break;
                }
                const FaceID fid_p = this->get_face_id(_pid);

                distance_calculator_on_surface<FaceID> dist_calc(
                    std::make_pair(found_p->second.position(), fid_p),
                    this->polygon());

                const Real dist = boost::apply_visitor(dist_calc, found_s->second) +
                                  found_p->second.radius();
                if(dist > 1e-8)
                {
                    result = false;
                    std::cerr << "ERROR: particle is not inside of the Single (ID="
                              << evid << ")\n";
                    std::cerr << "     : shell size   = "
                              << boost::apply_visitor(shell_size_getter(), found_s->second)
                              << '\n';
                    std::cerr << "     : shell pos    = "
                              << boost::apply_visitor(shell_position_getter(), found_s->second)
                              << '\n';
                    std::cerr << "     : particle pos = " << found_p->second.position()
                              << '\n';
                    std::cerr << "     : dist - r_shell + r_particle = " << dist
                              << '\n';
                }

                particles.erase(found_p);
                shells.erase(found_s);
                break;
            }
            case event_type::pair_domain:
            {
                const Pair&      pr   = boost::get<Pair>(dom);
                const ShellID    _shid = pr.shell_id();
                const ParticleID _pid0 = pr.particle_id_at(0);
                const ParticleID _pid1 = pr.particle_id_at(1);
                if(pid2evid.count(_pid0) == 0){pid2evid[_pid0] = evid;}
                if(pid2evid.count(_pid1) == 0){pid2evid[_pid1] = evid;}
                if(sid2evid.count(_shid) == 0){sid2evid[_shid] = evid;}

                BOOST_AUTO(found_p0, std::find_if(particles.begin(), particles.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                        ParticleID, Particle>(_pid0)));
                if(found_p0 == particles.end())
                {
                    result = false;
                    std::cerr << "ERROR: particle might assigned to two"
                              << "different domains\n";
                    std::cerr << "     : Pair domain " << evid << " has particle "
                              << _pid0 << " but the particle is already erased\n";
                    if(pid2evid.count(_pid0) == 1)
                    {
                        std::cerr << "     : event " << pid2evid[_pid0]
                                  << " has particle " << _pid0 << '\n';
                    }
                    break;
                }
                BOOST_AUTO(found_p1, std::find_if(particles.begin(), particles.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                        ParticleID, Particle>(_pid1)));
                if(found_p1 == particles.end())
                {
                    result = false;
                    std::cerr << "ERROR: particle might assigned to two"
                              << "different domains\n";
                    std::cerr << "     : Pair domain " << evid << " has particle "
                              << _pid1 << " but the particle is already erased\n";
                    if(pid2evid.count(_pid1) == 1)
                    {
                        std::cerr << "     : event " << pid2evid[_pid1]
                                  << " has particle " << _pid1 << '\n';
                    }
                    break;
                }
                BOOST_AUTO(found_s, std::find_if(shells.begin(), shells.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                        ShellID, shell_type>(_shid)));
                if(found_s == shells.end())
                {
                    result = false;
                    std::cerr << "ERROR: shell might assigned to two"
                              << "different domains\n";
                    std::cerr << "     : Pair domain " << evid << " has shell "
                              << _shid << " but the shell is already erased\n";
                    if(sid2evid.count(_shid) == 1)
                    {
                        std::cerr << "     : event " << sid2evid[_shid]
                                  << " has shell " << _shid << '\n';
                    }
                    break;
                }
                const FaceID fid_p0 = this->get_face_id(_pid0);
                const FaceID fid_p1 = this->get_face_id(_pid1);

                if(found_s->second.which() != shell_container_type::circular_shell)
                {
                    result = false;
                    std::cerr << "ERROR: currently, pair is only for circular\n";
                    std::cerr << "     : domain " << evid << "has shell " << _shid
                              << ", but it has invalid shape " << found_s->second.which()
                              << '\n';
                    break;
                }

                distance_calculator_on_surface<FaceID>
                    dist_calc0(std::make_pair(found_p0->second.position(), fid_p0),
                               this->polygon());
                distance_calculator_on_surface<FaceID>
                    dist_calc1(std::make_pair(found_p1->second.position(), fid_p1),
                               this->polygon());

                const Real dist0 = boost::apply_visitor(dist_calc0, found_s->second) +
                                   found_p0->second.radius();
                const Real dist1 = boost::apply_visitor(dist_calc1, found_s->second) +
                                   found_p1->second.radius();
                if(dist0 > 0)
                {
                    result = false;
                    std::cerr << "ERROR: particle " << _pid0
                              << " is not inside of the Pair " << evid << "\n";
                    std::cerr << "     : dist - r_shell + r_particle = " << dist0
                              << '\n';
                }
                if(dist1 > 0)
                {
                    result = false;
                    std::cerr << "ERROR: particle " << _pid1
                              << " is not inside of the Pair " << evid << "\n";
                    std::cerr << "     : dist - r_shell + r_particle = " << dist1
                              << '\n';
                }
                particles.erase(found_p0);
                // to avoid iterator break
                found_p1 = std::find_if(particles.begin(), particles.end(),
                    ecell4::utils::pair_first_element_unary_predicator<
                        ParticleID, Particle>(_pid1));
                particles.erase(found_p1);
                shells.erase(found_s);
                break;
            }
            case event_type::multi_domain:
            {
                const Multi& mul = boost::get<Multi>(dom);
                ShellID _shid;
                ParticleID _pid; Particle _p;
                BOOST_FOREACH(boost::tie(_pid, _p), mul.particles())
                {
                    if(pid2evid.count(_pid) == 0){pid2evid[_pid] = evid;}
                    bool within = false;
                    BOOST_AUTO(found_p, std::find_if(particles.begin(), particles.end(),
                        ecell4::utils::pair_first_element_unary_predicator<
                            ParticleID, Particle>(_pid)));
                    if(found_p == particles.end())
                    {
                        result = false;
                        std::cerr << "ERROR: particle might assigned to two "
                                  << "different domains\n";
                        std::cerr << "     : domain " << evid << " has particle "
                                  << _pid << " but the particle is already erased\n";
                        if(pid2evid.count(_pid) == 1)
                        {
                            std::cerr << "     : event " << pid2evid[_pid]
                                      << " has particle " << _pid << '\n';
                        }
                        continue;
                    }
                    const FaceID fid_p = this->get_face_id(_pid);
                    distance_calculator_on_surface<FaceID>
                        dist_calc(std::make_pair(found_p->second.position(), fid_p),
                                  this->polygon());

                    BOOST_FOREACH(_shid, mul.shell_ids())
                    {
                        BOOST_AUTO(found_s, std::find_if(shells.begin(), shells.end(),
                            ecell4::utils::pair_first_element_unary_predicator<
                                ShellID, shell_type>(_shid)));
                        if(found_s == shells.end())
                        {
                            result = false;
                            std::cerr << "ERROR: shell might assigned to two "
                                      << "different domains\n";
                            std::cerr << "     : Multi domain " << evid << " has shell "
                                      << _shid << " but the shell is already erased\n";
                            continue;
                        }
                        const Real dist =
                            boost::apply_visitor(dist_calc, found_s->second) +
                            found_p->second.radius();
                        if(dist < 0)
                        {
                            within = true;
                        }
                    }

                    if(!within)
                    {
                        result = false;
                        std::cerr << "ERROR: particle is not inside of any multi shell!\n";
                        std::cerr << "PID = " << _pid << ", DID = " << evid << '\n';
                    }
                    particles.erase(found_p);
                }
                BOOST_FOREACH(_shid, mul.shell_ids())
                {
                    if(sid2evid.count(_shid) == 0){sid2evid[_shid] = evid;}

                    BOOST_AUTO(found_s, std::find_if(shells.begin(), shells.end(),
                        ecell4::utils::pair_first_element_unary_predicator<
                            ShellID, shell_type>(_shid)));
                    if(found_s == shells.end())
                    {
                        result = false;
                        std::cerr << "ERROR: shell might assigned to two"
                                  << "different domains\n";
                        std::cerr << "     : Multi domain " << evid << " has shell "
                                  << _shid << " but the shell is already erased\n";
                        if(sid2evid.count(_shid) == 1)
                        {
                            std::cerr << "     : event " << sid2evid[_shid]
                                      << " has shell " << _shid << '\n';
                        }
                        continue;
                    }
                    shells.erase(found_s);
                }
                break;
            }
            case event_type::birth_domain:
            {
//                 std::cerr << "INFO : birth_domain is assigned\n";
                break;
            }
            default:
            {
                result = false;
                std::cerr << "ERROR: event " << evid
                          << " has invalid domain_kind " << ev_ptr->which_domain()
                          << '\n';
                break;
            }
        }
    }

    if(!particles.empty())
    {
        result = false;
        std::cerr << "ERROR: some of particles are not assigned to Domain\n";
        BOOST_FOREACH(boost::tie(pid, p), particles)
        {
            std::cerr << "     : particle id " << pid
                      << " is not assigned to any Domain\n";
        }
    }
    if(!shells.empty())
    {
        result = false;
        std::cerr << "ERROR: some of shells are not assigned to Domain\n";
        BOOST_FOREACH(boost::tie(shid, sh), shells)
        {
            std::cerr << "     : shell id " << shid
                      << " is not assigned to any Domain\n";
            const DomainID _did = boost::apply_visitor(domain_id_getter(), sh);
            std::cerr << "     : it should be contained by " << _did << '\n';
        }
    }
    if(result)
    {
//         const boost::chrono::steady_clock::time_point end_ =
//             boost::chrono::high_resolution_clock::now();
        std::cerr << "time = " << this->time() << " simulator is sanitized.\n";
//                   << "it took " << static_cast<double>(boost::chrono::duration_cast<
//                       boost::chrono::milliseconds>(end_ - start_).count()) / 1e3
//                   << " seconds.\n";
    }
    std::cerr << std::flush;
    return result;
}// }}}

} // sgfrd
} // ecell4
