#ifndef __ECELL4_BD_BD_PROPAGATOR_2D_HPP
#define __ECELL4_BD_BD_PROPAGATOR_2D_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>

#include "BDPropagator.hpp"
#include "BDPolygon.hpp"
#include "functions2d.hpp"

namespace ecell4
{

namespace bd
{

class BDPropagator2D
{
public:

    typedef ReactionInfo reaction_info_type;
    typedef BDWorld::molecule_info_type molecule_info_type;
    typedef ecell4::Polygon::FaceID FaceID;

public:

    BDPropagator2D(
        Model& model, BDWorld& world, RandomNumberGenerator& rng,
        const Real dt, const Real rl,
        std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions)
    : model_(model), world_(world), poly_(world.container_2D().polygon()),
      rng_(rng), dt_(dt), last_reactions_(last_reactions),
      max_retry_count_(1), reaction_length_(rl)
    {
        queue_ = world.container_2D().list_particles();
        shuffle(rng_, queue_);
    }

    bool operator()();

    inline Real dt() const
    {
        return dt_;
    }

    inline RandomNumberGenerator& rng()
    {
        return rng_;
    }

    bool attempt_reaction(const ParticleID& pid, const Particle& particle,
                          const FaceID& fid);
    bool attempt_reaction(
        const ParticleID& pid1, const Particle& particle1, const FaceID& f1,
        const ParticleID& pid2, const Particle& particle2, const FaceID& f2);

    void remove_particle(const ParticleID& pid);

    Real3 draw_displacement(const Particle& particle, const Real3& normal)
    {
        return random_displacement_2d(rng_, this->dt(), particle.D(), normal);
    }

    Real3 draw_ipv(const Real r, const Real D)
    {
        return random_ipv_2d(rng_, r, D, reaction_length_);
    }

    Real3 draw_ipv(const Real r, const Real D, const Real3& normal)
    {
        return random_ipv_2d(rng_, r, D, reaction_length_, normal);
    }

private:

    class particle_finder
        : public std::unary_function<std::pair<ParticleID, Particle>, bool>
    {
    public:

        particle_finder(const ParticleID& pid)
            : pid_(pid)
        {
            ;
        }

        bool operator()(std::pair<ParticleID, Particle> pid_particle_pair)
        {
            return (pid_particle_pair.first == pid_);
        }

    protected:

        ParticleID pid_;
    };

    Real calc_reaction_area(const Real radius_sum) const
    {
        const Real rad_react(radius_sum + this->reaction_length_);
        return M_PI * (rad_react * rad_react - radius_sum * radius_sum);
    }

protected:

    Model& model_;
    BDWorld& world_;
    const BDPolygon& poly_; // XXX additional
    RandomNumberGenerator& rng_;
    Real dt_;
    std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions_;
    Integer max_retry_count_;
    Real reaction_length_;

    BDWorld::particle_container_type queue_;
};

} // bd

} // ecell4

#endif /* __ECELL4_BD_BD_PROPAGATOR_HPP */
