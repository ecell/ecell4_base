#ifndef DOMAIN_FACTORY_HPP
#define DOMAIN_FACTORY_HPP

template<typename Ttraits_>
struct DomainFactory
{
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::shell_id_pair shell_id_pair;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::single_type single_type;
    typedef typename traits_type::pair_type pair_type;
    typedef typename traits_type::domain_id_type domain_id_type;
    typedef typename traits_type::network_rules_type network_rules_type;

    virtual single_type*
    create_single(domain_id_type const& domain_id,
                  particle_id_pair const& particle_id_pair,
                  shell_id_pair const& shell_id_pair,
                  typename network_rules_type::reaction_rule_vector const& reactions) const = 0;

    virtual pair_type*
    create_pair(domain_id_type const& domain_id,
                position_type const& com,
                particle_id_pair const& single1,
                particle_id_pair const& single2,
                shell_id_pair const& shell_id_pair,
                typename network_rules_type::reaction_rule_vector const& reactions) const = 0;
};

#endif /* DOMAIN_FACTORY_HPP */
