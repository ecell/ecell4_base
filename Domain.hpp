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
    typedef typename traits_type::world_type::structure_id_type structure_id_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::event_id_type event_id_type;
    typedef typename traits_type::time_type time_type;

public:
    virtual ~Domain() {}

    Domain(structure_id_type const& structure_id)
        : structure_id_(structure_id) {}

    structure_id_type const& structure_id() const
    {
        return structure_id_;
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

protected:
    structure_id_type structure_id_;
    event_id_type event_id_;
    time_type last_time_;
    time_type dt_;
};

#endif /* DOMAIN_HPP */
