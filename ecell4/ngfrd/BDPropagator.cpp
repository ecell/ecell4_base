#include <ecell4/ngfrd/BDPropagator.hpp>
#include <ciso646>

namespace ecell4
{
namespace ngfrd
{

// ============================================================================
// 2D propagation

void BDPropagator::propagate_2D_particle(
        const ParticleID& pid, Particle p, FaceID fid)
{
    if(attempt_single_reaction_2D(pid, p, fid))
    {
        return; // reaction happened. done.
    }
    if(p.D() == Real(0.0))
    {
        return; // particle does not move. done.
    }

    const auto prev_pos = std::make_pair(p.position(), fid);

    auto disp    = this->draw_2D_displacement(p, fid);
    auto new_pos = ecell4::polygon::travel(world_.polygon(), prev_pos, disp);

    if(not is_inside_of_shells_2D(new_pos, p.radius()))
    {
        // determine positions of particles in overlapping shells.
        sim_.determine_positions(new_pos, p.radius());

        if(world_.has_overlapping_particle(new_pos, p.radius(), /*ignore = */ pid))
        {
            // if overlap exists, the movement would be rejected.

            // XXX this includes the case where the particle collides with a
            //     particle in a multi shell, not a particle at outside...
            ++rejected_move_count_;
            return;
        }
    }

    // collect particles within reactive range
    auto overlapped = world_.list_particles_within_radius(
            new_pos, p.radius() + reaction_length_, /*ignore = */ pid);

    bool core_overlapped = false;
    for(const auto& pidpd : overlapped)
    {
        // {{ParticleID, Particle}, Real(distance)}
        if(pidpd.second < p.radius() + pidpd.first.second.radius())
        {
            core_overlapped = true;
            break;
        }
    }

    if(not core_overlapped)
    {
        // movement accepted.
        p.position() = new_pos.first;
        fid          = new_pos.second;
        world_.upadte_particle(pid, p, fid);
    }

    if(core_overlapped)
    {
        // reject this displacement. But still, it can attempt pair reaction
        // if there is a reactive partner around the original position.

        // Update overlapping particles.
        overlapped = world_.list_particles_within_radius(
            prev_pos, p.radius() + reaction_length_, /*ignore = */ pid);
    }

    if(overlapped.empty())
    {
        return; // No reactive partner, we don't need to attempt pair reaction.
    }

    attempt_pair_reaction_2D(pid, p, fid, overlapped);
    return ;
}

// ============================================================================
// 2D single reaction

// return true if reaction happens
bool BDPropagator::attempt_single_reaction_2D(const ParticleID& pid,
        const Particle& p, const FaceID& fid)
{
    const auto& rules = this->model_.query_reaction_rules(p.species());
    if(rules.empty())
    {
        return false;
    }

    Real probability = 0.0;
    const Real threshold = this->rng_.uniform(0.0, 1.0);
    for(const auto& rule : rules)
    {
        probability += rule.k() * dt_;
        if(probability <= threshold)
        {
            continue;
        }
        if(1.0 < probability)
        {
            std::cerr << "BDPropagator::attempt_single_reaction_2D: "
                         "probability exceeds 1.0" << std::endl;
        }
        switch(rule.products().size())
        {
            case 0:
            {
                world_.remove_particle(pid);
                last_reactions_.emplace_back(rule,
                        make_degradation_reaction_info(world_.t(), pid, p));
                return true;
            }
            case 1:
            {
                return attempt_1to1_reaction_2D(pid, p, fid, rule);
            }
            case 2:
            {
                return attempt_1to2_reaction_2D(pid, p, fid, rule);
            }
            default:
            {
                throw_exception<NotSupported>("ngfrd::BDPropagator attempts 1 "
                    "to ", rule.products().size(), " reaction which is not "
                    "supported. => ", rule.as_string());
            }
        }
    }
    // No reaction rule is triggered. Nothing happens.
    return false;
}

// return true if reaction is accepted.
bool BDPropagator::attempt_1to1_reaction_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const ReactionRule& rule)
{
    assert(rule.products().size() == 1);
    const auto species_new = rule.products().front();
    const auto molinfo     = container_.get_molecule_info(species_new);
    const Real radius_new  = molinfo.radius;
    const Real D_new       = molinfo.D;

    const auto pos_fid = std::make_pair(p.position(), fid);

    if(not is_inside_of_shells_2D(pos_fid, radius_new))
    {
        sim_.determine_positions(pos_fid, radius_new);
    }

    if(world_.has_overlapping_particle(pos_fid, radius_new, /*ignore = */ pid))
    {
        return false;
    }

    Particle particle_new(species_new, p.position(), radius_new, D_new);

    world_.update_particle(pid, particle_new, fid);

    last_reactions_.emplace_back(rule, make_unimolecular_reaction_info(
                world_.t(), pid, p, pid, particle_new));
    return true;
}

// return true if reaction is accepted.
bool BDPropagator::attempt_1to2_reaction_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const ReactionRule& rule)
{
    assert(rule.products().size() == 2);

    const auto sp1      = rule.products().at(0);
    const auto sp2      = rule.products().at(1);
    const auto molinfo1 = world_.get_molecule_info(sp1);
    const auto molinfo2 = world_.get_molecule_info(sp2);

    const Real D1  = molinfo1.D;
    const Real D2  = molinfo2.D;
    const Real r1  = molinfo1.radius;
    const Real r2  = molinfo2.radius;
    const Real D12 = D1 + D2;
    const Real r12 = r1 + r2;

    if(D1 == 0. && D2 == 0)
    {
        throw_exception<NotSupported>("ngfrd::BDPropagator attempts 1 to 2 "
            "reaction but both particle is immovable => ", rule.as_string());
    }

    const Real3 n = world_.polygon().triangle_at(fid).normal();

    std::pair<Real3, FaceID> pos1_new = std::make_pair(p.position(), fid);
    std::pair<Real3, FaceID> pos2_new = std::make_pair(p.position(), fid);

    const Polygon& polygon = world_.polygon();

    const Real separation_length = r12 * NGFRDSimulator::SAFETY; // (1 + 1e-5)
    std::size_t separation_count = 1 + max_retry_count_;
    while(separation_count != 0)
    {
        --separation_count;

        const Real3 ipv = draw_ipv_2D(separation_length, D12, n);
        Real3 disp1 = ipv * (D1 / D12);
        Real3 disp2 = disp1 - ipv; // disp1 + (-disp2) = ipv

        pos1_new = ecell4::polygon::travel(polygon, pos1_new, disp1);
        pos2_new = ecell4::polygon::travel(polygon, pos2_new, disp2);

        // check distance on polygon
        const Real dist = ecell4::polygon::distance(polygon, pos1_new, pos2_new);
        if(dist <= r12)
        {
            pos1_new = std::make_pair(p.position(), fid); // rollback
            pos2_new = std::make_pair(p.position(), fid);

            continue;
        }

        if(world_.has_overlapping_particle(pos1_new, r1, /*ignore = */ pid) ||
           world_.has_overlapping_particle(pos2_new, r2, /*ignore = */ pid))
        {
            continue;
        }
    }
    if(separation_count == 0)
    {
        return false; // could not find an appropreate configuration.
    }

    // ------------------------------------------------------------------------
    // check configurations

    if(not is_inside_of_shells_2D(pos1_new, r1))
    {
        sim_.determine_positions(pos1_new, r1);
    }
    if(not is_inside_of_shells_2D(pos2_new, r2))
    {
        sim_.determine_positions(pos2_new, r2);
    }

    if(world_.has_overlapping_particle(pos1_new, r1, /*ignore = */ pid) ||
       world_.has_overlapping_particle(pos2_new, r2, /*ignore = */ pid))
    {
        this->rejected_move_count_ += 1;
        return false; // overlaps with a particle at outside of the domain
    }

    // ------------------------------------------------------------------------
    // reaction is accepted! update particle positions.

    Particle p1_new(sp1, pos1_new.first, r1, D1);
    Particle p2_new(sp2, pos2_new.first, r2, D2);

    const auto result1 = world_.update_particle(pid, p1_new, pos1_new.second);
    const auto result2 = world_.new_particle   (     p2_new, pos2_new.second);

    assert(not result1); // should be already exist (it returns true if it's new)
    assert(result2.second); // should succeed

    const auto pid2 = result2.first.first;

    last_reactions_.emplace_back(rule, make_unbinding_reaction_info(world_.t(),
                pid, p, pid, p1_new, pid2, p2_new));

    // ------------------------------------------------------------------------
    // trial move.

    const auto try_to_move_particle =
        [&](const ParticleID& pid_, Particle p_, const FaceID fid_) {
            auto  pos  = std::make_pair(p_.position(), fid_);
            Real3 disp = draw_2D_displacement(p_, fid_);
            auto  new_pos = ecell4::polygon::travel(polygon, pos, disp);

            if(not is_inside_of_shells_2D(new_pos, p_.radius()))
            {
                sim_.determine_positions(new_pos, p_.radius());
            }
            if(not world_.has_overlapping_particle(new_pos, p_.radius(), pid_))
            {
                p_.position() = new_pos.first;
                world_.update_particle(pid_, p_, new_pos.second);
            }
            return ;
        };

    // move 1 or 2 particles, randomly. If 2 particles are ready to move,
    // the order is also randomized.
    switch(rng_.uniform_int(0, 3))
    {
        case 0:
        {
            try_to_move_particle(pid,  p1_new, pos1_new.second);
            break;
        }
        case 1:
        {
            try_to_move_particle(pid2, p2_new, pos2_new.second);
            break;
        }
        case 2:
        {
            try_to_move_particle(pid,  p1_new, pos1_new.second);
            try_to_move_particle(pid2, p2_new, pos2_new.second);
            break;
        }
        case 3:
        {
            try_to_move_particle(pid2, p2_new, pos2_new.second);
            try_to_move_particle(pid,  p1_new, pos1_new.second);
            break;
        }
    }
    return true;
}

// ============================================================================
// 2D pair reaction

bool BDPropagator::attempt_pair_reaction_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const std::vector<std::pair<std::pair<ParticleID, Particle>, Real>& overlapped)
{
    Real probability = 0.0;
    const Real threshold = rng_.uniform(0.0, 1.0);

    for(const auto& pidp : overlapped) // try all the particles in the reactive region
    {
        const auto& rules = model_.query_reaction_rules(p1.species(), p2.species());
        if(rules.empty())
        {
            continue;
        }

        const Real k_tot = std::accumurate(rules.begin(), rules.end(), Real(0.0),
                [](const Real k, const ReactionRule& rule) -> Real {
                    return k + rule.k();
                });

        probability += k_tot * calc_pair_acceptance_coef_2D(p, pidp.second);

        if(probability <= threshold)
        {
            continue;
        }

        if(1.0 < probability)
        {
            std::cerr << "BDPropagator::attempt_pair_reaction_2D: "
                         "probability exceeds 1.0" << std::endl;
        }

        const auto& rule = determine_reaction_rule(rules, k_tot);

        switch(rule.products().size())
        {
            case 0:
            {
                world_.remove_particle(pid);
                world_.remove_particle(pidp.first);

                const auto found = std::find(queue_.begin(), queue_.end(), pidp.first);
                if(found != queue_.end())
                {
                    queue_.erase(found);
                }

                last_reactions_.emplace_back(rule, make_degradation_reaction_info(
                            world_.t(), pid, p, pidp.first, pidp.second);
                return true;
            }
            case 1:
            {
                const auto fid2 = world_.on_which_face(pidp.frist);
                assert(fid2.has_value());

                return attempt_2to1_reaction_2D(p, pid, fid,
                        pidp.first, pidp.second, *fid2, rule);
            }
            default:
            {
                throw_exception<NotSupported>("ngfrd::BDPropagator: attempts "
                        "2 to ", rule.products().size(), " reaction which is "
                        "not allowed. => ", rule.as_string());
            }
        }
    }
    return false;
}

bool BDPropagator::attempt_2to1_reaction_2D(
        const ParticleID& pid1, const Particle& p1, const FaceID& fid1,
        const ParticleID& pid2, const Particle& p2, const FaceID& fid2,
        const ReactionRule& rule)
{
    const auto species_new = rule.products().front();
    const auto molinfo     = world_.get_molecule_info(species_new);
    const Real radius_new  = molinfo.radius;
    const Real D_new       = molinfo.D;

    const auto pos1 = std::make_pair(p1.position(), fid1);
    const auto pos2 = std::make_pair(p2.position(), fid2);
    const Real D1   = p1.D();
    const Real D2   = p2.D();
    const Real D12  = D1 + D2;

    std::pair<Real3, FaceID> new_pos;
    if(D1 == 0.0)
    {
        new_pos = std::make_pair(pos1, fid1);
    }
    else if(D2 == 0.0)
    {
        new_pos = std::make_pair(pos2, fid2);
    }
    else
    {
        const Real3 dp = ecell4::polygon::direction(
                world_.polygon(), pos1, pos2) * (D1 / D12);
        new_pos = ecell4::polygon::travel(world_.polygon(), pos1, dp);
    }

    if(not is_inside_of_shells_2D(new_pos, radius_new))
    {
        sim_.determine_positions(new_pos, radius_new);
    }
    if(world_.has_overlapping_particle(new_pos, radius_new, pid1, pid2))
    {
        return false;
    }

    Particle particle_new(species_new, new_pos.first, radius_new, D_new);

    world_.remove_particle(pid2);
    world_.update_particle(pid1, particle_new, new_pos.second);

    const auto found2 = std::find(queue_.begin(), queue_.end(), pid2);
    if(found2 != queue_.end())
    {
        queue_.erase(found2);
    }

    last_reactions_.emplace_back(rule, make_binding_reaction_info(world_.t(),
                pid1, p1, pid2, p2, pid1, particle_new));
    return true;
}

} // ngfrd
} // ecell4
