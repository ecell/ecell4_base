#ifndef PAIR_HPP
#define PAIR_HPP

#include <cmath>
#include <boost/array.hpp>
#include "Domain.hpp"

template<typename Ttraits_>
class Pair: public Domain<Ttraits_>
{
public:
    typedef Domain<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::world_type::surface_id_type surface_id_type;
    typedef boost::array<particle_id_pair, 2> particle_array_type;

public:
    virtual ~Pair() {}

    Pair(surface_id_type const& surface_id,
         particle_id_pair const& p0, particle_id_pair const& p1)
        : base_type(surface_id)
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

protected:
    mutable particle_array_type particles_;
};

#endif /* PAIR_HPP */
