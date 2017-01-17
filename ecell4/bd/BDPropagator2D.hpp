#ifndef __ECELL4_BD_BD_PROPAGATOR_2D_HPP
#define __ECELL4_BD_BD_PROPAGATOR_2D_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>

#include "BDPropagator.hpp"
#include "BDPolygon.hpp"
#include "rotate_vector.hpp"

namespace ecell4
{

namespace bd
{

class BDPropagator2D
{
public:

    typedef ReactionInfo reaction_info_type;
    typedef BDPolygon polygon_type;
    typedef typename BDPolygon::face_type face_type;

public:

    BDPropagator2D(
        Model& model, BDWorld& world, RandomNumberGenerator& rng, const Real& dt,
        std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions)
        : model_(model), world_(world), poly_(world.polygon()), rng_(rng), dt_(dt),
        last_reactions_(last_reactions), max_retry_count_(1)
    {
        queue_ = world.list_2D_particles();
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

    bool attempt_reaction(const ParticleID& pid, const Particle& particle);
    bool attempt_reaction(
        const ParticleID& pid1, const Particle& particle1,
        const ParticleID& pid2, const Particle& particle2);

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

    void remove_particle(const ParticleID& pid);

    inline Real3 draw_displacement(const Particle& particle, const Real3& normal)
    {
//      return random_displacement_3d(rng(), dt(), particle.D());
        assert(length_sq(normal) - 1.0 < 1e-12);
        const Real sigma = std::sqrt(2 * particle.D() * dt());
        const Real x = rng_.gaussian(sigma);
        const Real y = rng_.gaussian(sigma);
        const Real leng  = std::sqrt(x * x + y * y);

        const Real theta = rng_.random() * 2. * M_PI;
        const Real3 dr = Real3(normal[1], -normal[0], 0) /
            std::sqrt(normal[1] * normal[1] + normal[0] * normal[0]);
        const Real3 disp = rotate(theta, normal, dr);

        if(dot_product(normal, disp) > 1e-10)
        {
            std::cerr << "large dot product" << dot_product(normal, disp) << std::endl;
            std::cerr << "theta  = " << theta << std::endl;
            std::cerr << "normal = " << normal[0] << ", " << normal[1] << ", " << normal[2] << std::endl;
            std::cerr << "disp   = " << disp[0] << ", " << disp[1] << ", " << disp[2] << std::endl;
        }

        return disp;
    }

    inline Real3 draw_ipv(const Real& sigma, const Real& t, const Real& D)
    {// XXX!
        return random_ipv_3d(rng(), sigma, t, D);
    }

protected:

    Model& model_;
    BDWorld& world_;
    BDPolygon& poly_; // XXX additional
    RandomNumberGenerator& rng_;
    Real dt_;
    std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions_;
    Integer max_retry_count_;

    BDWorld::particle_container_type queue_;
};

} // bd

} // ecell4

#endif /* __ECELL4_BD_BD_PROPAGATOR_HPP */
