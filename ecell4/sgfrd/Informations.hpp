#ifndef ECELL4_SGFRD_REACTION_AND_MOLECULE_INFO
#define ECELL4_SGFRD_REACTION_AND_MOLECULE_INFO
#include <ecell4/core/Particle.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/ReactionRule.hpp>

namespace ecell4
{
namespace sgfrd
{

struct MoleculeInfo
{
    MoleculeInfo(const Real r, const Real d): radius(r), D(d){}
    MoleculeInfo(MoleculeInfo const& rhs): radius(rhs.radius), D(rhs.D){}
    const Real radius;
    const Real D;
};

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

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_REACTION_AND_MOLECULE_INFO
