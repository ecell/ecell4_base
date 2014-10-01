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
    typedef std::pair<const shell_id_type, shell_type> shell_id_pair;
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

    shell_id_pair& shell()
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

    virtual position_type& position()
    {
        return shape_position(shape(shell_.second));
    }

    virtual length_type const& size() const
    {
        return shape_size(shape(shell_.second));
    }

    virtual length_type& size()
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

    virtual void accept(ImmutativeDomainVisitor<traits_type> const& visitor) const
    {
        visitor(*this);
    }

    virtual void accept(MutativeDomainVisitor<traits_type> const& visitor)
    {
        visitor(*this);
    }

    virtual std::string as_string() const
    {
        return (boost::format(
            "%s(id=%s, event=%s, last_time=%.16g, dt=%.16g, particle=(%s:%s), shell=(%d:%s))") %
            type_name() %
            boost::lexical_cast<std::string>(base_type::id_) %
            boost::lexical_cast<std::string>(base_type::event_.first) %
            base_type::last_time_ % base_type::dt_ %
            boost::lexical_cast<std::string>(base_type::particle().first) %
            boost::lexical_cast<std::string>(base_type::particle().second) %
            boost::lexical_cast<std::string>(shell_.first) %
            boost::lexical_cast<std::string>(shell_.second)).str();
    }

protected:
    shell_id_pair shell_;
};

#endif /* ANALYTICAL_SINGLE_HPP */
