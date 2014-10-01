#ifndef ANALYTICAL_PAIR_HPP
#define ANALYTICAL_PAIR_HPP

#include <cmath>
#include <boost/array.hpp>
#include "Pair.hpp"
#include "AnalyticalSingle.hpp"

template<typename Ttraits_, typename Tshell_>
class AnalyticalPair: public Pair<Ttraits_>
{
public:
    typedef Pair<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::world_type::traits_type::D_type D_type;
    typedef typename traits_type::domain_id_type identifier_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename network_rules_type::reaction_rule_type reaction_rule_type;
    typedef typename network_rules_type::reaction_rule_vector reaction_rule_vector;
    typedef Tshell_ shell_type;
    typedef std::pair<const shell_id_type, shell_type> shell_id_pair;

public:
    virtual ~AnalyticalPair() {}

    AnalyticalPair(identifier_type const& id,
                   particle_id_pair const& p0, particle_id_pair const& p1,
                   shell_id_pair const& shell,
                   position_type const& iv,
                   reaction_rule_vector const& reactions)
        : base_type(id, p0, p1), shell_(shell), iv_(iv), reactions_(reactions)
    {
        // determine a_r and a_R
        {
            D_type D0(base_type::particles_[0].second.D());
            D_type D1(base_type::particles_[1].second.D());
            length_type R0(base_type::particles_[0].second.radius());
            length_type R1(base_type::particles_[1].second.radius());
            const length_type sigma(R0 + R1);
            const length_type D_tot(D0 + D1);
            const length_type D_geom(std::sqrt(D0 * D1));
            const length_type shell_size(shape(shell.second).radius() / traits_type::SAFETY);
            const length_type r0(this->r0());
            BOOST_ASSERT(r0 >= sigma);
            if (((D_geom - D0) * r0) / D_tot + shell_size
                + std::sqrt(D0 / D1) * (R1 - shell_size) - R0 < 0)
            {
                std::swap(D0, D1);
                std::swap(R0, R1);
            }
            a_R_ = D_geom * (D0 * (shell_size - R1)
                              + D1 * (shell_size - r0 - R1)) /
                   (D1 * D1 + D1 * D0 + D_geom * D_tot);
            a_r_ = (D_geom * r0 + D_tot * (shell_size - R1)) / (D1 + D_geom);
            BOOST_ASSERT(a_r_ > 0);
            BOOST_ASSERT(a_r_ > r0);
            BOOST_ASSERT(a_R_ > 0 || (a_R_ == 0. && (D1 == 0. || D0 == 0.)));
            BOOST_ASSERT(a_R_ + a_r_ * D1 / D_tot + R1 >=
                         a_R_ + a_r_ * D0 / D_tot + R0);
            BOOST_ASSERT(std::abs(a_R_ + a_r_ * D1 / D_tot + R1 - shell_size) <
                         1e-12 * shell_size);
        }
    }

    shell_id_pair const& shell() const
    {
        return shell_;
    }

    shell_id_pair& shell()
    {
        return shell_;
    }

    position_type const& iv() const
    {
        return iv_;
    }

    length_type r0() const
    {
        return length(iv_);
    }

    length_type const& a_R() const
    {
        return a_R_;
    }

    length_type const& a_r() const
    {
        return a_r_;
    }

    length_type sigma() const
    {
        return base_type::particles_[0].second.radius()
               + base_type::particles_[1].second.radius();
    }

    D_type D_tot() const
    {
        return base_type::particles_[0].second.D()
               + base_type::particles_[1].second.D();
    }

    D_type D_geom() const
    {
        return std::sqrt(
            base_type::particles_[0].second.D() *
            base_type::particles_[1].second.D());
    }

    D_type D_R() const
    {
        return base_type::particles_[0].second.D() *
               base_type::particles_[1].second.D() / D_tot();
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

    virtual char const* type_name() const
    {
        return retrieve_domain_type_name(*this);
    }

    reaction_rule_vector const& reactions() const
    {
        return reactions_;
    }

    virtual typename Domain<traits_type>::size_type num_shells() const
    {
        return 1;
    }

    virtual typename Domain<traits_type>::size_type multiplicity() const
    {
        return 2;
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
            "%s(id=%s, event=%s, last_time=%.16g, dt=%.16g, particles=[(%s:%s), (%s:%s)], iv=%s, shell=(%s:%s))") %
            type_name() %
            boost::lexical_cast<std::string>(base_type::id_) %
            boost::lexical_cast<std::string>(base_type::event_.first) %
            base_type::last_time_ % base_type::dt_ %
            boost::lexical_cast<std::string>(base_type::particles()[0].first) %
            boost::lexical_cast<std::string>(base_type::particles()[0].second) %
            boost::lexical_cast<std::string>(base_type::particles()[1].first) %
            boost::lexical_cast<std::string>(base_type::particles()[1].second) %
            boost::lexical_cast<std::string>(iv_) %
            boost::lexical_cast<std::string>(shell_.first) %
            boost::lexical_cast<std::string>(shell_.second)).str();
    }
protected:
    shell_id_pair shell_;
    position_type const iv_;
    reaction_rule_vector const& reactions_;
    mutable length_type a_R_;
    mutable length_type a_r_;
};

#endif /* ANALYTICAL_PAIR_HPP */
