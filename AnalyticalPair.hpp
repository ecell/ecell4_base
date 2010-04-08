#ifndef ANALYTICAL_PAIR_HPP
#define ANALYTICAL_PAIR_HPP

#include <cmath>
#include <boost/array.hpp>
#include "Pair.hpp"

template<typename Ttraits_, typename Tshell_>
class AnalyticalPair: public Pair<Ttraits_>
{
public:
    typedef Pair<Ttraits_> base_type;
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

public:
    virtual ~AnalyticalPair() {}

    AnalyticalPair(surface_id_type const& surface_id,
                   particle_id_pair const& p0, particle_id_pair const& p1,
                   shell_id_pair const& shell,
                   length_type const& r0, length_type const& rt)
        : base_type(surface_id, p0, p1), shell_(shell),
          r0_(r0), rt_(rt)
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

protected:
    const shell_id_pair shell_;
    const length_type r0_;
    const length_type rt_;
    mutable length_type a_R_;
    mutable length_type a_r_;
    static const double SAFETY = 1.0 + 1e-5;
};

#endif /* ANALYTICAL_PAIR_HPP */
