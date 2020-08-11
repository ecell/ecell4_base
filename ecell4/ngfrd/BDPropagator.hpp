#ifndef ECELL4_NGFRD_BD_PROPAGATOR_HPP
#define ECELL4_NGFRD_BD_PROPAGATOR_HPP
#include <ecell4/core/exception.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/ngfrd/ReactionInfo.hpp>
#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/Shell.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>
#include <ecell4/ngfrd/NGFRDSimulator.hpp>

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
Real Igbd_3d(const Real sigma, const Real t, const Real D) noexcept
{
    constexpr Real sqrtPi = boost::math::constants::root_pi<Real>();

    const Real Dt      = D * t;
    const Real Dt2     = Dt + Dt;
    const Real sqrtDt  = std::sqrt(Dt);
    const Real sigmasq = sigma * sigma;

    constexpr Real term1 = 1 / (3 * sqrtPi);
    const     Real term2 = sigmasq - Dt2;
    const     Real term3 = Dt2 - 3 * sigmasq;
    const     Real term4 = sqrtPi * sigmasq * sigma *
                           boost::math::erfc<Real>(sigma / sqrtDt);

    return term1 * (-sqrtDt * (term2 * std::exp(-sigmasq / Dt) + term3) + term4);
}

Real Igbd_r_3d(const Real r, const Real sigma, const Real t, Real D) noexcept
{
    constexpr Real one_div_root_pi = boost::math::constants::one_div_root_pi<Real>();

    const Real Dt  = D * t;
    const Real Dt2 = Dt + Dt;
    const Real Dt4 = Dt2 + Dt2;
    const Real sqrtDt  = std::sqrt(Dt);
    const Real sqrtDt4 = 2 * sqrtDt;
    const Real sigmasq = sigma * sigma;
    const Real sigmacb = sigmasq * sigma;
    const Real rcb     = r * r * r;

    const Real rsigma =  r * sigma;
    const Real rps_sq = (r + sigma) * (r + sigma);
    const Real rms_sq = (r - sigma) * (r - sigma);

    const Real term1 = -2 * sqrtDt * one_div_root_pi;
    const Real term2 =  std::exp(-sigmasq / Dt) * (sigmasq - Dt2);
    const Real term3 = -std::exp(-rps_sq / Dt4) * (rms_sq + rsigma - Dt2);
    const Real term4 =  std::exp(-rms_sq / Dt4) * (rps_sq - rsigma - Dt2);
    const Real term5 = -sigmasq * 3 + Dt2;

    const Real term6 =  (sigmacb - rcb)     * boost::math::erf((r - sigma) / sqrtDt4);
    const Real term7 = -(sigmacb + sigmacb) * boost::math::erf( sigma      / sqrtDt );
    const Real term8 =  (sigmacb + rcb)     * boost::math::erf((r + sigma) / sqrtDt4);

    return (term1 * (term2 + term3 + term4 + term5) + term6 + term7 + term8) / 6;
}

Real drawR_gbd_3D(const Real sigma, const Real t, const Real D, const Real rnd) noexcept
{
    constexpr Real abs_tol = 1e-18;
    constexpr Real rel_tol = 1e-12;
    const Real target = Igbd_3d(sigma, t, D) * rnd; // rnd is in [0.0, 1.0)

    const Real low  = sigma;
    const Real high = sigma + 10 * std::sqrt(6 * D * t);

    const auto tol = [=](const Real a, const Real b) noexcept {
        return std::abs(a - b) < abs_tol || std::abs(a / b - 1.0) < rel_tol;
    };
    const auto Reqn = [=](const Real x) noexcept {
        return Igbd_r_3D(x, sigma, t, D) - target;
    };

    constexpr std::size_t max_iteration = 100;
    std::size_t iterated = max_iteration;

    const std::pair<Real, Real> t_range =
        boost::math::tools::toms748_solve(Reqn, low, high, tol, iterated);

    if(iterated == max_iteration)
    {
        throw_exception<std::runtime_error>(
            "ngfrd::BDPropagator::bd_math::drawR_gbd_3D: did not converge. "
            "(rnd=", rnd, ", r12=", sigma, ", dt=", t, ", D12=", D, ")")
    }
    return t_range.first;
}
} // bd_math

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

        Particle p = world_.get_particle(pid).second;
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
        cosnt Real3 rxy(len * std::cos(theta), len * std::sin(theta), 0.0);

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
        const Real scale(8 * r * std::sqrt(0.25 - r2));
        return Real3(a * scale, b * scale, r * (8 * r2 - 1));
    }

    Real3 draw_ipv_3D(const Real r12, const Real dt, const Real D12)
    {
        return random_spherical_uniform(
                bd_math::drawR_gbd_3D(r12, dt, D12, rng.uniform(0.0, 1.0)));
    }

    Real calc_pair_acceptance_coef_2D(
            const Particle& p1, const Particle& p2) const noexcept;

    ReactionRule const& determine_reaction_rule(
            const std::vector<ReactionRule>& rules, const Real k_tot) noexcept;

private:

    Model&                 model_;
    NGFRDWorld&            world_;
    NGFRDSimulator&        sim_;
    RandomNumberGenerator& rng_;
    Real                   dt_;
    Real                   reaction_length_;
    std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions_;
    std::vector<ParticleID> queue_;
    std::vector<ShellID>    shells_;
    std::size_t rejected_move_count_;
    std::size_t max_retry_count_;
};

} // egfrd
} // ecell4
#endif /* BD_PROPAGATOR_HPP */
