#ifndef ECELL4_NGFRD_REACTION_INFO_HPP
#define ECELL4_NGFRD_REACTION_INFO_HPP
#include <ecell4/core/Particle.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <boost/container/static_vector.hpp>

namespace ecell4
{
namespace ngfrd
{

class ReactionInfo
{
public:

    using reactant_container_type = boost::container::static_vector<
        std::pair<ParticleID, Particle>, 2>;

    using product_container_type  = boost::container::static_vector<
        std::pair<ParticleID, Particle>, 2>;

public:

    explicit ReactionInfo(const Real t)
        : t_(t)
    {}
    ReactionInfo(const Real t,
                 const std::vector<std::pair<ParticleID, Particle>>& reactants,
                 const std::vector<std::pair<ParticleID, Particle>>& products)
        : t_(t), reactants_(reactants.begin(), reactants.end()),
                 products_ (products.begin(), products.end())
    {}
    ~ReactionInfo() = default;
    ReactionInfo(const ReactionInfo&) = default;
    ReactionInfo(ReactionInfo&&)      = default;
    ReactionInfo& operator=(const ReactionInfo&) = default;
    ReactionInfo& operator=(ReactionInfo&&)      = default;

    Real t() const noexcept
    {
        return t_;
    }

    reactant_container_type const& reactants() const noexcept
    {
        return reactants_;
    }

    void add_reactant(const std::pair<ParticleID, Particle>& pidp)
    {
        reactants_.push_back(pidp);
    }

    product_container_type const& products() const noexcept
    {
        return products_;
    }

    void add_product(const std::pair<ParticleID, Particle>& pidp)
    {
        products_.push_back(pidp);
    }

protected:

    Real t_;
    reactant_container_type reactants_;
    product_container_type  products_;
};

inline ReactionInfo
make_degradation_reaction_info(const Real t,
        const ParticleID& pid, const Particle& p)
{
    return ReactionInfo(t, {std::make_pair(pid, p)}, {});
}

inline ReactionInfo
make_degradation_reaction_info(const Real t,
        const ParticleID& pid1, const Particle& p1,
        const ParticleID& pid2, const Particle& p2)
{
    return ReactionInfo(t,
            {std::make_pair(pid1, p1), std::make_pair(pid2, p2)}, {});
}

inline ReactionInfo
make_birth_reaction_info(
        const Real t, const ParticleID& pid, const Particle& p)
{
    return ReactionInfo(t, {}, {std::make_pair(pid, p)});
}

inline ReactionInfo
make_unimolecular_reaction_info(const Real t,
        const ParticleID& pid1, const Particle& p1,
        const ParticleID& pid2, const Particle& p2)
{
    return ReactionInfo(t, {std::make_pair(pid1, p1)},
                           {std::make_pair(pid2, p2)});
}

inline ReactionInfo
make_binding_reaction_info(const Real t,
        const ParticleID& pid1, const Particle& p1,
        const ParticleID& pid2, const Particle& p2,
        const ParticleID& pid3, const Particle& p3)
{
    return ReactionInfo(t, {std::make_pair(pid1, p1), std::make_pair(pid2, p2)},
                           {std::make_pair(pid3, p3)});
}

inline ReactionInfo
make_unbinding_reaction_info(const Real t,
        const ParticleID& pid1, const Particle& p1,
        const ParticleID& pid2, const Particle& p2,
        const ParticleID& pid3, const Particle& p3)
{
    return ReactionInfo(t, {std::make_pair(pid1, p1)},
            {std::make_pair(pid2, p2), std::make_pair(pid3, p3)});
}

} // ngfrd
} // ecell4
#endif// ECELL4_NGFRD_REACTION_AND_MOLECULE_INFO_HPP
