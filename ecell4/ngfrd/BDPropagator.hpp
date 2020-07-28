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

    // {distance, parameter `t` for the segment (p + td)}
    Real distance_segment_triangle(
            const Real3& p, const Real3& disp, const Triangle& t)
    {
        // Here we don't need to consider the periodic images. Those are
        // correctly handled by the core algorithm.

        // First, check if the segment intersects triangle.
        {
            const auto n = t.normal(); // always |n| = 1
            const auto dot_v0p_d = dot_product(p - t.vertices()[0], n);

            assert(dot_v0p_d != 0) // already collides

            Real3 q;       // where it collides
            Bacycentric b; // where it collides (parametric)
            if(0 < dot_v0p_d) // front side.
            {
                // intersect_ray_triangle checks only if 0 < n*d.
                // We need to FLIP front/back
                Triangle tri(t.vertices()[0], t.vertices()[2], t.vertices()[1]);
                if(intersect_ray_triangle(p, disp, tri, b, q))
                {
                    // collides within the displacement range?
                    if(length_sq(q - p) < length_sq(disp))
                    {
                        return 0.0;
                    }
                    // else, do nothing. go ahead.
                }
            }
            else // back side
            {
                Triangle tri(t.vertices()[0], t.vertices()[1], t.vertices()[2]);
                if(intersect_ray_triangle(p, disp, tri, b, q))
                {
                    if(length_sq(q - p) < length_sq(disp))
                    {
                        return 0.0;
                    }
                    // else, do nothing. go ahead.
                }
            }
        }
        // It does not collides with the triangle.
        // Next, check the minimum distance between segment and triangle.
        Real mindist_sq = std::numeric_limits<Real>::max();

        // - distance between Segment PQ / edge AB
        {
            Real  s,  t;
            Real3 c1, c2;
            const Real dist_sq = closest_point_segment_segment(p, p + disp,
                    t.vertices()[0], t.vertices()[1], s, t, c1, c2);
            mindist_sq = std::min(mindist_sq, dist_sq);
        }
        // - distance between Segment PQ / edge BC
        {
            Real  s,  t;
            Real3 c1, c2;
            const Real dist_sq = closest_point_segment_segment(p, p + disp,
                    t.vertices()[1], t.vertices()[2], s, t, c1, c2);
            mindist_sq = std::min(mindist_sq, dist_sq);
        }
        // - distance between Segment PQ / edge CA
        {
            Real  s,  t;
            Real3 c1, c2;
            const Real dist_sq = closest_point_segment_segment(p, p + disp,
                    t.vertices()[2], t.vertices()[0], s, t, c1, c2);
            mindist_sq = std::min(mindist_sq, dist_sq);
        }
        // - distance between point P / Triangle
        {
            mindist_sq = length_sq(p - closest_point_point_triangle(p, t));
        }
        // - distance between point Q / Triangle
        {
            mindist_sq = length_sq(p + disp -
                    closest_point_point_triangle(p + disp, t));
        }
        return std::sqrt(mindist_sq);
    }

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
