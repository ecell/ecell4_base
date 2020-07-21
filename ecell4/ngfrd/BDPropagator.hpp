#ifndef ECELL4_NGFRD_BD_PROPAGATOR_HPP
#define ECELL4_NGFRD_BD_PROPAGATOR_HPP
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/ngfrd/ReactionInfo.hpp>
#include <greens_functions/freeFunctions.hpp>
#include <boost/math/constants/constants.hpp>

namespace ecell4
{
namespace ngfrd
{

class BDPropagator
{
public:

    BDPropagator(Model& model, NGFRDWorld& world, NGFRDSimulator& sim,
        RandomNumberGenerator& rng, const Real& dt,
        const std::size_t max_retry_count,
        std::vector<ParticleID> particles, std::vector<ShellID> shells,
        std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions)
        : model_(model), world_(world), sim_(sim), rng_(rng), dt_(dt),
          last_reactions_(last_reactions), queue_(particles), shells_(shells),
          rejected_move_count_(0), max_retry_count_(max_retry_count)
    {
        shuffle(rng, queue_);
    }

    bool operator()()
    {
        if(queue_.empty())
        {
            return false;
        }
        const auto pid = queue_.back();
        queue_.pop_back();

        Particle p = world_.get_particle(pid);
        if(auto fid = world_.on_which_face(pid))
        {
            this->propagate_2D_particle(pid, std::move(p), std::move(*fid));
        }
        else // 3D particle
        {
            this->propagate_3D_particle(pid, std::move(p));
        }
        return true;
    }

    std::size_t get_rejected_move_count() const
    {
        return rejected_move_count_;
    }

private:

    void propagate_2D_particle(const ParticleID&, Particle, FaceID);
    void propagate_3D_particle(const ParticleID&, Particle);

    bool attempt_single_reaction_2D(const ParticleID&, const Particle&, const FaceID&);
    bool attempt_1to1_reaction_2D(
        const ParticleID&, const Particle&, const FaceID&, const ReactionRule&);
    bool attempt_1to2_reaction_2D(
        const ParticleID&, const Particle&, const FaceID&, const ReactionRule&);

    bool attempt_pair_reaction_2D(const ParticleID&, const Particle&, const FaceID&,
        const std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>& overlapped);
    bool attempt_2to1_reaction_2D(
        const ParticleID&, const Particle&, const FaceID&,
        const ParticleID&, const Particle&, const FaceID&,
        const ReactionRule&);

    bool attempt_single_reaction_3D(const ParticleID&, const Particle&);
    bool attempt_1to1_reaction_3D(
        const ParticleID&, const Particle&, const FaceID&, const ReactionRule&);
    bool attempt_1to2_reaction_3D(
        const ParticleID&, const Particle&, const FaceID&, const ReactionRule&);

    bool attempt_pair_reaction_3D(const ParticleID&, const Particle&,
                                  const ParticleID&, const Particle&);
    bool attempt_2to1_reaction_3D(const ParticleID&, const Particle&,
                                  const ParticleID&, const Particle&,
                                  const ReactionRule&);

    bool is_inside_of_shells_3D(const Real3&,                    const Real& radius) const;
    bool is_inside_of_shells_2D(const std::pair<Real3, FaceID>&, const Real& radius) const;

    Real3 draw_2D_displacement(const Particle&, const FaceID&);
    Real3 draw_3D_displacement(const Particle&);

    Real calc_pair_acceptance_coef_2D(
            const Particle& p1, const Particle& p2) const noexcept
    {
        const Real radius_sum     = p1.radius() + p2.radius();
        const Real reaction_range = radius_sum + reaction_length_;

        const Real reaction_area = boost::math::constants::pi<Real>() *
            (radius_sum * radius_sum - reaction_range * reaction_range);

        if(p1.D() == 0.0 || p2.D() == 0)
        {
            // immovable particles immediately return after attempting 1st order
            // reaction.
            // to attempt 2nd order reaction with them, we need to double the
            // acceptance coefficient.
            return this->dt_ / reaction_area;
        }
        else
        {
            // movable particles checks 2nd order reaction. If the both reactants
            // are movable, both particle attempts reaction. So here it halves
            // the acceptance coefficient to avoid double-counting.
            return 0.5 * this->dt_ / reaction_area;
        }
    }
    ReactionRule const&
    determine_reaction_rule(const std::vector<ReactionRule>& rules,
                            const Real k_tot) noexcept
    {
        assert(!rules.empty());
        if(rules.size() == 1)
        {
            return rules.front();
        }
        const Real rnd = this->rng_.uniform(0.0, 1.0) * k_tot;
        Real k_cumm = 0.0;
        for(const auto& rule : rules)
        {
            k_cumm += rule.k();
            if(rnd < k_cumm)
            {
                return rule;
            }
        }
        return rules.back();
    }

private:

    Model&                 model_;
    NGFRDWorld&            world_;
    NGFRDSimulator&        sim_;
    RandomNumberGenerator& rng_;
    Real                   dt_;
    std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions_;
    std::vector<ParticleID> queue_;
    std::vector<ShellID>    shells_;
    std::size_t rejected_move_count_;
    std::size_t max_retry_count_;
};

} // egfrd
} // ecell4
#endif /* BD_PROPAGATOR_HPP */
