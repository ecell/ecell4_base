#ifndef ECELL4_NGFRD_BD_PROPAGATOR_HPP
#define ECELL4_NGFRD_BD_PROPAGATOR_HPP
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/ngfrd/ReactionInfo.hpp>
#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/Shell.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>

#include <greens_functions/freeFunctions.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/roots.hpp>

namespace ecell4
{
namespace ngfrd
{

// determine the distance between particles that are dissociated from each other
namespace bd_math
{
Real drawR_gbd_3D(const Real sigma, const Real t, const Real D, const Real rnd) noexcept;
} // bd_math

class NGFRDSimulator; // forward

class BDPropagator
{
public:

    BDPropagator(const Model& model, NGFRDWorld& world, NGFRDSimulator& sim,
        RandomNumberGenerator& rng, const Real& dt,
        const std::size_t max_retry_count,
        std::vector<ParticleID> particles, std::vector<ShellID> shells,
        std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions)
        : model_(model), world_(world), sim_(sim), rng_(rng), dt_(dt),
          last_reactions_(last_reactions), particles_(particles),
          queue_(particles), shells_(shells),
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

        Particle p = world_.get_particle(pid).second;
        if(auto fid = world_.on_which_face(pid))
        {
//             std::cerr << "BDPropagator: moving 2D particle " << pid << std::endl;
            this->propagate_2D_particle(pid, std::move(p), std::move(*fid));
        }
        else // 3D particle
        {
//             std::cerr << "BDPropagator: moving 3D particle " << pid << std::endl;
            this->propagate_3D_particle(pid, std::move(p));
        }
        return true;
    }

    std::size_t get_rejected_move_count() const
    {
        return rejected_move_count_;
    }

    // to wrap resulting particles by domains
    std::vector<ParticleID> particles() const noexcept
    {
        return particles_;
    }

    bool particle_escaped() const
    {
        for(const auto& pid : this->particles_)
        {
            const auto& p = world_.get_particle(pid).second;
            if(auto fid = world_.on_which_face(pid))
            {
                if(!is_inside_of_shells_2D(std::make_pair(p.position(), *fid), p.radius()))
                {
                    return true;
                }
            }
            else
            {
                if(!is_inside_of_shells_3D(p.position(), p.radius()))
                {
                    return true;
                }
            }
        }
        return false;
    }

private:

    void remove_particle(const ParticleID& pid)
    {
        // remove from internal particle Id cache
        const auto fnd = std::remove(particles_.begin(), particles_.end(), pid);
        assert(std::distance(fnd, particles_.end()) == 1);
        this->particles_.erase(fnd, particles_.end());
        // remove from (outer-) world
        this->world_.remove_particle(pid);
        return;
    }

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
        const ParticleID&, const Particle&, const ReactionRule&);
    bool attempt_1to2_reaction_3D(
        const ParticleID&, const Particle&, const ReactionRule&);

    bool attempt_pair_reaction_3D(const ParticleID&, const Particle&,
                                  const ParticleID&, const Particle&);
    bool attempt_2to1_reaction_3D(const ParticleID&, const Particle&,
                                  const ParticleID&, const Particle&,
                                  const ReactionRule&);

    bool is_inside_of_shells_3D(const Real3&, const Real& radius) const;
    bool is_inside_of_shells_2D(
            const std::pair<Real3, FaceID>&, const Real& radius) const;

    Real3 draw_2D_displacement(const Particle&, const FaceID&);
    Real3 draw_3D_displacement(const Particle&);

    Real3 random_circular_uniform(const Real len, const Real3& normal)
    {
        assert(std::abs(length(normal) - 1.0) < 1e-12); // |n| == 1
        constexpr Real pi = boost::math::constants::pi<Real>();

        const Real  theta = this->rng_.uniform(0.0, 2.0 * pi);
        const Real3 rxy(len * std::cos(theta), len * std::sin(theta), 0.0);

        // dot(z, n), z:=(0,0,1)
        const Real  phi  = std::acos(boost::algorithm::clamp(normal[2], -1.0, 1.0));

        if     (std::abs(phi)      < 1e-12) {return rxy;}
        else if(std::abs(phi - pi) < 1e-12) {return rxy * (-1.0);}
        else if(std::abs(phi + pi) < 1e-12) {return rxy * (-1.0);}

        const Real3 axis(-normal[1], normal[0], 0); // cross(z, n)
        return rotate(phi, axis * (1.0 / length(axis)), rxy);
    }
    Real3 draw_ipv_2D(const Real r, const Real D, const Real3& normal)
    {
        const Real rl    = r + this->reaction_length_;
        const Real r_sq  = r * r;
        const Real rd    = rl * rl - r_sq;
        const Real ipvl  = std::sqrt(r_sq + this->rng_.uniform(0., 1.) * rd);
        return random_circular_uniform(ipvl, normal);
    }

    Real3 random_spherical_uniform(const Real len)
    {
        Real a(0), b(0), r2(1);
        while (r2 > 0.25)
        {
            a = rng_.uniform(0, 1) - 0.5;
            b = rng_.uniform(0, 1) - 0.5;
            r2 = a * a + b * b;
        }
        const Real scale(8 * len * std::sqrt(0.25 - r2));
        return Real3(a * scale, b * scale, len * (8 * r2 - 1));
    }

    Real3 draw_ipv_3D(const Real r12, const Real dt, const Real D12)
    {
        return random_spherical_uniform(
                bd_math::drawR_gbd_3D(r12, dt, D12, rng_.uniform(0.0, 1.0)));
    }

    Real calc_pair_acceptance_coef_2D(
            const Particle& p1, const Particle& p2) const noexcept;

    ReactionRule const& determine_reaction_rule(
            const std::vector<ReactionRule>& rules, const Real k_tot) noexcept;

private:

    const Model&           model_;
    NGFRDWorld&            world_;
    NGFRDSimulator&        sim_;
    RandomNumberGenerator& rng_;
    Real                   dt_;
    Real                   reaction_length_;
    std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions_;
    std::vector<ParticleID> particles_; // resulting particles
    std::vector<ParticleID> queue_;
    std::vector<ShellID>    shells_;
    std::size_t rejected_move_count_;
    std::size_t max_retry_count_;
};

} // ngfrd
} // ecell4
#endif /* BD_PROPAGATOR_HPP */
