#ifndef PAIR_HPP
#define PAIR_HPP

#include <cmath>
#include <boost/array.hpp>
#include "Domain.hpp"

template<typename Ttraits_, typename Tshell_>
class Pair: public Domain<Ttraits_>
{
public:
    typedef Domain<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::world_type::surface_id_type surface_id_type;
    typedef typename traits_type::world_type::traits_type::D_type D_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef Tshell_ shell_type;
    typedef std::pair<shell_id_type, shell_type> shell_id_pair;
    typedef typename network_rules_type::reaction_rule_vector reaction_rule_vector;
    typedef boost::array<particle_id_pair, 2> particle_array_type;

public:
    virtual ~Pair() {}

    Pair(surface_id_type const& surface_id,
         shell_id_pair const& shell, position_type const& com,
         particle_id_pair const& p0, particle_id_pair const& p1,
         length_type const& r0, length_type const& rt)
        : base_type(surface_id), shell_(shell),
          r0_(r0), rt_(rt)
    {
        if (p0.second.D() < p1.second.D())
        {
            new(&particles_[0]) particle_id_pair(p0);
            new(&particles_[1]) particle_id_pair(p1);
        }
        else
        {
            new(&particles_[0]) particle_id_pair(p1);
            new(&particles_[1]) particle_id_pair(p0);
        }
        // determine a_r and a_R
        {
            D_type D0(particles_[0].second.D());
            D_type D1(particles_[1].second.D());
            length_type R0(particles_[0].second.radius());
            length_type R1(particles_[1].second.radius());
            const length_type sigma(R0 + R1);
            const length_type D_tot(D0 + D1);
            const length_type D_geom(std::sqrt(D0 * D1));
            const length_type shell_size(shape(shell.second).radius() / SAFETY);
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

    particle_array_type const& particles() const
    {
        return particles_;
    }

    shell_id_pair const& shell() const
    {
        return shell_;
    }

    length_type const& r0() const
    {
        return r0_;
    }

    length_type const& rt() const
    {
        return rt_;
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
        return particles_[0].second.radius() + particles_[1].second.radius();
    }

    D_type D_tot() const
    {
        return particles_[0].second.D() + particles_[1].second.D();
    }

    D_type D_geom() const
    {
        return std::sqrt(particles_[0].second.D() * particles_[1].second.D());
    }

    D_type D_R() const
    {
        return particles_[0].second.D() * particles_[1].second.D() / D_tot();
    }

protected:
    const shell_id_pair shell_;
    mutable particle_array_type particles_;
    const length_type r0_;
    const length_type rt_;
    mutable length_type a_R_;
    mutable length_type a_r_;
    static const double SAFETY = 1.0 + 1e-5;
};

#endif /* PAIR_HPP */
