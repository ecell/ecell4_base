#ifndef ECELL4_BD_BD_PROPAGATOR_HPP
#define ECELL4_BD_BD_PROPAGATOR_HPP

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>

#include "functions3d.hpp"
#include "BDWorld.hpp"


namespace ecell4
{

namespace bd
{

class ReactionInfo
{
public:

    typedef std::pair<ParticleID, Particle> particle_id_pair_type;
    typedef std::vector<particle_id_pair_type> container_type;

public:

    ReactionInfo(
        const Real t,
        const container_type& reactants,
        const container_type& products)
        : t_(t), reactants_(reactants), products_(products)
    {}

    ReactionInfo(const ReactionInfo& another)
        : t_(another.t()), reactants_(another.reactants()), products_(another.products())
    {}

    Real t() const
    {
        return t_;
    }

    const container_type& reactants() const
    {
        return reactants_;
    }

    void add_reactant(const particle_id_pair_type& pid_pair)
    {
        reactants_.push_back(pid_pair);
    }

    const container_type& products() const
    {
        return products_;
    }

    void add_product(const particle_id_pair_type& pid_pair)
    {
        products_.push_back(pid_pair);
    }

protected:

    Real t_;
    container_type reactants_, products_;
};

class BDPropagator
{
public:

    typedef ReactionInfo reaction_info_type;

public:

    BDPropagator(
        Model& model, BDWorld& world, RandomNumberGenerator& rng, const Real& dt,
        std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions)
        : model_(model), world_(world), rng_(rng), dt_(dt),
        last_reactions_(last_reactions), max_retry_count_(1)
    {
        queue_ = world_.container_3D().list_particles();
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

    inline Real3 draw_displacement(const Particle& particle)
    {
        return random_displacement_3d(rng(), dt(), particle.D());
    }

    inline Real3 draw_ipv(const Real& sigma, const Real& t, const Real& D)
    {
        return random_ipv_3d(rng(), sigma, t, D);
    }

protected:

    Model& model_;
    BDWorld& world_;
    RandomNumberGenerator& rng_;
    Real dt_;
    std::vector<std::pair<ReactionRule, reaction_info_type> >& last_reactions_;
    Integer max_retry_count_;

    BDWorld::particle_container_type queue_;
};

} // bd

} // ecell4

#endif /* ECELL4_BD_BD_PROPAGATOR_HPP */
