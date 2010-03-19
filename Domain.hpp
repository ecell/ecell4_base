#ifndef DOMAIN_HPP
#define DOMAIN_HPP

template<typename Ttraits_>
class Domain
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::world_type::surface_id_type surface_id_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::domain_id_type domain_id_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::event_id_type event_id_type;
    typedef typename traits_type::time_type time_type;

public:
    virtual ~Domain() {}

    Domain(domain_id_type const& id, surface_id_type const& surface_id)
        : id_(id), surface_id_(surface_id) {}

    domain_id_type const& id() const
    {
        return id_;
    }

    surface_id_type const& surface_id() const
    {
        return surface_id_;
    }

    event_id_type const& event_id() const
    {
        return event_id_;
    }

    event_id_type& event_id()
    {
        return event_id_;
    }

    time_type const& last_time() const
    {
        return last_time_;
    }

    time_type& last_time()
    {
        return last_time_;
    }

    time_type const& dt() const
    {
        return dt_;
    }

    time_type& dt()
    {
        return dt_;
    }

    event_kind_type const& event_kind() const
    {
        return event_kind_;
    }

    event_kind_type& event_kind()
    {
        return event_kind_;
    }

protected:
    domain_id_type id_;
    surface_id_type surface_id_;
    event_id_type event_id_;
    time_type last_time_;
    time_type dt_;
    event_kind_type event_kind_;
};

#endif /* DOMAIN_HPP */
