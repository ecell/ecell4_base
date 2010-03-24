#ifndef SINGLE_HPP
#define SINGLE_HPP

#include "Domain.hpp"

template<typename Ttraits_, typename Tshell_>
class Single: public Domain<Ttraits_>
{
public:
    typedef Domain<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::world_type::surface_id_type surface_id_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef Tshell_ shell_type;
    typedef std::pair<shell_id_type, shell_type> shell_id_pair;
    typedef typename network_rules_type::reaction_rule_vector reaction_rule_vector;
    typedef typename traits_type::rate_type rate_type;

public:
    virtual ~Single() {}

    Single(surface_id_type const& surface_id,
           shell_id_pair const& shell,
           particle_id_pair const& particle,
           reaction_rule_vector const& reactions)
        : base_type(surface_id), shell_(shell),
          particle_(particle), reactions_(reactions),
          k_tot_(calculate_k_tot(reactions)) {}

    particle_id_pair const& particle() const
    {
        return particle_;
    }

    shell_id_pair const& shell() const
    {
        return shell_;
    }

    reaction_rule_vector const& reactions()
    {
        return reactions_;
    }

    rate_type const& k_tot()
    {
        return k_tot_;
    }

private:
    static rate_type calculate_k_tot(reaction_rule_vector const& reactions)
    {
        rate_type k_tot(0);
        for (typename reaction_rule_vector::const_iterator
             i(reactions.begin()), e(reactions.end());
             i != e; ++i)
        {
            k_tot += (*i).k();
        }
        return k_tot;
    }

protected:
    const shell_id_pair shell_;
    const particle_id_pair particle_;
    reaction_rule_vector const& reactions_;
    const rate_type k_tot_;
};

#endif /* SINGLE_HPP */
