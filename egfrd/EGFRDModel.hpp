#ifndef __EGFRD_MODEL_HPP
#define __EGFRD_MODEL_HPP

// #include "get_mapper_mf.hpp"

#include <ecell4/core/types.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "ParticleModel.hpp" // epdp

namespace ecell4
{

namespace egfrd
{

class EGFRDModel
    : public Model
{
public:

    typedef NetworkModel::species_container_type species_container_type;
    typedef NetworkModel::reaction_rule_container_type
    reaction_rule_container_type;

protected:

    typedef utils::get_mapper_mf<Species::serial_type, ::SpeciesTypeID>::type
    species_type_id_map_type;

public:

    EGFRDModel()
        : network_model_(), particle_model_()
    {
        ;
    }

    ::ParticleModel& particle_model()
    {
        return particle_model_;
    }

    species_container_type const& species() const
    {
        return network_model_.species();
    }

    reaction_rule_container_type const& reaction_rules() const
    {
        return network_model_.reaction_rules();
    }

    std::vector<ReactionRule> query_reaction_rules(Species const& sp) const;
    std::vector<ReactionRule> query_reaction_rules(
        Species const& sp1, Species const& sp2) const;

    void add_species(Species const& sp);
    void remove_species(Species const& sp);
    bool has_species(Species const& sp) const;

    void add_reaction_rule(ReactionRule const& rr);
    void remove_reaction_rule(ReactionRule const& rr);
    bool has_reaction_rule(ReactionRule const& rr) const;

    ::SpeciesTypeID get_species_type_id(Species const& sp)
    {
        species_type_id_map_type::iterator i(sid_map_.find(sp.serial()));
        if (i == sid_map_.end())
        {
            throw NotFound("species not found");
        }
        return (*i).second;
    }

protected:

    NetworkModel network_model_;
    ::ParticleModel particle_model_;
    species_type_id_map_type sid_map_;
};

} // egfrd

} // ecell4

#endif /* __EGFRD_MODEL_HPP */
