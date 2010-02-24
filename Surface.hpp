#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "Vector3.hpp"
#include "Single.hpp"
#include "Pair.hpp"
#include "Multi.hpp"

template<typename Tsimulator_>
class Surface
{
public:
    typedef Tsimulator_ simulator_type;
    typedef simulator_type::traits_type traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::D_type D_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::world_type::shell_id_pair shell_id_pair;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::domain_type domain_type;
    typedef typename traits_type::single_type single_type;
    typedef typename traits_type::pair_type pair_type;
    typedef typename traits_type::network_rules_type network_rules_type;

public:

    virtual ~Surface() {}

    std::string const& name() const
    {
        return name_;
    }

    std::string& name()
    {
        return name_;
    }

    simulator_type* simulator()
    {
        return simulator_;
    }

    virtual length_type distance(position_type const& pos) const = 0;

    virtual position_type draw_bd_displacement(time_type const& dt, D_type const& D) = 0;

    virtual single_type*
    create_single(domain_id_type const& domain_id,
                  particle_id_pair const& particle_id_pair,
                  shell_id_pair const& shell_id_pair,
                  typename network_rules_type::reaction_rule_vector const& reactions) const = 0;

    virtual pair_type* create_pair(domain_id_type const& domain_id,
                                   position_type const& com,
                                   particle_id_pair const& single1,
                                   particle_id_pair const& single2,
                                   shell_id_pair const& shell_id_pair,
                                   typename network_rules_type::reaction_rule_vector const& reactions) const = 0;

    Surface(simulator_type* sim, std::string const& name)
        : sim_(sim), name_(name) {}

protected:
    simulator_type* sim_;
    std::string name_;
};

template<typename Ttraits_, Tpos_>
inline typename Ttraits_::length_type distance(Surface<Ttraits_> const& p1, Tpos_ const& p2)
{
    p1.distance(p2);
}

#endif /* SURFACE_HPP */
