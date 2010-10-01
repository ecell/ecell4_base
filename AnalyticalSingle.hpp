#ifndef ANALYTICAL_SINGLE_HPP
#define ANALYTICAL_SINGLE_HPP

#include "Single.hpp"

template<typename Ttraits_, typename Tshell_>
class AnalyticalSingle: public Single<Ttraits_>
{
public:
    typedef Single<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::domain_id_type identifier_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef Tshell_ shell_type;
    typedef std::pair<shell_id_type, shell_type> shell_id_pair;
    typedef typename network_rules_type::reaction_rule_vector reaction_rule_vector;
    typedef typename traits_type::rate_type rate_type;

public:
    virtual ~AnalyticalSingle() {}

    AnalyticalSingle(identifier_type const& id,
                     particle_id_pair const& particle,
                     shell_id_pair const& shell)
        : base_type(id, particle), shell_(shell) {}

    shell_id_pair const& shell() const
    {
        return shell_;
    }

    length_type mobility_radius() const
    {
        return shape_size(shape(shell_.second)) - base_type::particle().second.radius();
    }

    virtual char const* type_name() const
    {
        return retrieve_domain_type_name(*this);
    }

    virtual position_type const& position() const
    {
        return shape_position(shape(shell_.second));
    }

    virtual length_type const& size() const
    {
        return shape_size(shape(shell_.second));
    }

    virtual typename Domain<traits_type>::size_type num_shells() const
    {
        return 1;
    }

    virtual typename Domain<traits_type>::size_type multiplicity() const
    {
        return 1;
    }

protected:
    const shell_id_pair shell_;
};

#endif /* ANALYTICAL_SINGLE_HPP */
