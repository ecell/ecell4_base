#ifndef ECELL4_SPATIOCYTE_REACTIONS_HPP
#define ECELL4_SPATIOCYTE_REACTIONS_HPP

#include <boost/shared_ptr.hpp>
#include <ecell4/core/VoxelPool.hpp>
#include <ecell4/core/ReactionRule.hpp>

namespace ecell4
{

namespace spatiocyte
{

class ReactionInfo
{
public:

    typedef std::pair<ParticleID, Voxel> particle_id_pair_type;
    typedef std::vector<particle_id_pair_type> container_type;

public:

    ReactionInfo() : t_(0), reactants_(), products_() {}

    ReactionInfo(const Real t) : t_(t), reactants_(), products_() {}

    ReactionInfo(
        const Real t,
        const container_type& reactants,
        const container_type& products)
        : t_(t), reactants_(reactants), products_(products) {}

    ReactionInfo(const ReactionInfo& another)
        : t_(another.t()), reactants_(another.reactants()), products_(another.products()) {}

    Real t() const
    {
        return t_;
    }

    bool has_occurred() const
    {
        return reactants_.size() > 0 || products_.size() > 0;
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

// Application of reactions

class SpatiocyteWorld;

ReactionInfo apply_a2b(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p,
        const Species& product_species);

ReactionInfo apply_a2bc(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p,
        const Species& product_species0,
        const Species& product_species1);

ReactionInfo apply_second_order_reaction(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionRule& reaction_rule,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1);

ReactionInfo apply_vanishment(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1);

ReactionInfo apply_ab2c(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1,
        const Species& product_species);

// ReactionInfo apply_ab2cd_in_order(
//         boost::shared_ptr<SpatiocyteWorld> world,
//         const ReactionInfo::particle_id_pair_type& p0,
//         const ReactionInfo::particle_id_pair_type& p1,
//         const Species& product_species0,
//         const Species& product_species1,
//         const SpatiocyteWorld::coordinate_type coord0,
//         const SpatiocyteWorld::coordinate_type coord1);

ReactionInfo apply_ab2cd(
        boost::shared_ptr<SpatiocyteWorld> world,
        const ReactionInfo::particle_id_pair_type& p0,
        const ReactionInfo::particle_id_pair_type& p1,
        const Species& product_species0,
        const Species& product_species1);

} // spatiocyte

} // ecell4

#endif /* ECELL4_SPATIOCYTE_REACTIONS_HPP */
