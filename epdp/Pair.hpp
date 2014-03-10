#ifndef PAIR_HPP
#define PAIR_HPP

#include <cmath>
#include <utility>
#include <boost/array.hpp>
#include "ShapedDomain.hpp"

template<typename Ttraits_>
class Pair: public ShapedDomain<Ttraits_>
{
public:
    typedef ShapedDomain<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::particle_id_type particle_id_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::domain_id_type identifier_type;
    typedef boost::array<particle_id_pair, 2> particle_array_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::traits_type::position_type position_type;

public:
    virtual ~Pair() {}

    Pair(identifier_type const& id,
         particle_id_pair const& p0, particle_id_pair const& p1)
        : base_type(id)
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
    }

    particle_array_type const& particles() const
    {
        return particles_;
    }

    particle_array_type& particles()
    {
        return particles_;
    }

protected:
    particle_array_type particles_;
};

#endif /* PAIR_HPP */
