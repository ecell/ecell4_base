#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <string>
#include <cstddef>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

template<typename Ttraits_>
class ImmutativeDomainVisitor;

template<typename Ttraits_>
class MutativeDomainVisitor;

template<typename Ttraits_>
class Domain
{
public:
    typedef Ttraits_ traits_type;
    typedef std::size_t size_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::domain_id_type identifier_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::event_id_pair_type event_id_pair_type;
    typedef typename traits_type::time_type time_type;

public:
    virtual ~Domain() {}

    Domain(identifier_type const& id)
        : id_(id), last_time_(0.), dt_(0.) {}

    identifier_type const& id() const
    {
        return id_;
    }

    event_id_pair_type const& event() const
    {
        return event_;
    }

    event_id_pair_type& event()
    {
        return event_;
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

    virtual size_type num_shells() const = 0;

    virtual size_type multiplicity() const = 0;

    virtual char const* type_name() const = 0;

    virtual std::string as_string() const
    {
        return (boost::format(
            "%s(id=%s, event=%s, last_time=%.16g, dt=%.16g)") %
            type_name() %
            boost::lexical_cast<std::string>(id_).c_str() %
            boost::lexical_cast<std::string>(event_.first).c_str() %
            last_time_ % dt_).str();
    }

    virtual void accept(ImmutativeDomainVisitor<traits_type> const&) const = 0;

    virtual void accept(MutativeDomainVisitor<traits_type> const&) = 0;

protected:
    identifier_type id_;
    event_id_pair_type event_;
    time_type last_time_;
    time_type dt_;
};

template<typename Tstrm, typename Ttraits, typename TdomTraits>
inline std::basic_ostream<Tstrm, Ttraits>&
operator<<(std::basic_ostream<Tstrm, Ttraits>& lhs,
           Domain<TdomTraits> const& rhs)
{
    lhs << rhs.as_string();
    return lhs;
}

template<typename Ttraits>
inline char const* retrieve_domain_type_name(Domain<Ttraits> const&)
{
    return "Domain";
}

#endif /* DOMAIN_HPP */
