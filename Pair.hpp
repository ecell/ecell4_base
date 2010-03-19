#ifndef PAIR_HPP
#define PAIR_HPP

#include <cmath>
#include "Domain.hpp"
#include "twofold_container.hpp"

template<typename Ttraits_>
class Pair: public Domain<Ttraits_>
{
public:
    typedef Domain<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::world_type::surface_id_type surface_id_type;
    typedef typename traits_type::world_type::D_type D_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::domain_id_type domain_id_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef Tshell_ shell_type;
    typedef std::pair<shell_id_type, shell_type> shell_id_pair;
    typedef typename network_rules_type::reaction_rule_vector reaction_rule_vector;
    typedef boost::array<particle_id_pair, 2> particle_array_type;

public:
    virtual ~Pair() {}

    Pair(domain_id_type const& id, surface_id_type const& surface_id,
         shell_id_pair const& shell, position_type const& com,
         particle_id_pair const& p0, particle_id_pair const& p1,
         length_type const& r0, length_type const& rt)
        : base_type(id, surface_id), shell_(shell),
          r0_(r0), rt_(rt)
    {
        if (p0.D() < p1.D())
        {
            new(&particles_[0]) particle_id_pair(p0);
            new(&particles_[1]) particle_id_pair(p1);
        }
        else
        {
            new(&particles_[0]) particle_id_pair(p1);
            new(&particles_[1]) particle_id_pair(p0);
        }

        D_tot_ = D0 + D1;
        D_geom_ = std::sqrt(D0 * D1);
        sigma_ = R0 + R1;

        // determine a_r and a_R
        {
            D_type D0(particles_[0].D);
            D_type D1(particles_[1].D);
            length_type R0(particles_[0].radius);
            length_type R1(particles_[1].radius);
            const length_type shell_size(shape(shell).radius() / SAFETY)
            BOOST_ASSERT(r0 >= sigma_);
            if (((D_geom_ - D1) * r0) / D_tot_ + shell_size
                + std::sqrt(D1 / D0) * (R0 - shell_size) - R1 < 0)
            {
                std::swap(D0, D1);
                std::swap(R0, R1);
            }
            a_R_ = D_geom_ * (D1 * (shell_size - R0)
                              + D0 * (shell_size - r0 - R0)) /
                   (D0 * D0 + D0 * D1 + D_geom_ * D_tot_);
            a_r_ = (D_geom_ * r0 + D_tot * (shell_size - R0)) /
                   (D0 + D_geom_);
            BOOST_ASSERT(a_r_ > 0);
            BOOST_ASSERT(a_r_ > r0);
            BOOST_ASSERT(a_R_ > 0 || (a_R_ == 0. && (D0 == 0. || D1 == 0.)));
            BOOST_ASSERT(a_R_ + a_r_ * D0 / D_tot + R0 >=
                         a_R_ + a_r_ * D1 / D_tot + R1);
            BOOST_ASSERT(std::abs(a_R_ + a_r_ * D0 / D_tot + R0 - shell_size) <
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
        return a_R;
    }

    length_type const& a_r() const
    {
        return a_r;
    }

    length_type const& sigma() const
    {
        return sigma_;
    }

    D_type const& D_tot() const
    {
        return D_tot_;
    }

    D_type const& D_geom() const
    {
        return D_geom_;
    }

    D_type const& D_R() const
    {
        return particles_[0].D() * particles_[1].D() / D_tot_;
    }

protected:
    const particle_array_type particles_;
    const shell_id_pair shell_;
    const length_type r0_;
    const length_type rt_;
    const length_type a_R_;
    const length_type a_r_;
    const length_type sigma_;
    const D_type D_tot_;
    const D_type D_geom_;
    static const double SAFETY = 1.0 + 1e-5;
};

#endif /* PAIR_HPP */
