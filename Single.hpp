#ifndef SINGLE_HPP
#define SINGLE_HPP

#include "Domain.hpp"

template<typename Ttraits_>
class Single: public Domain<Ttraits_>
{
public:
    typedef Domain<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::world_type::surface_id_type surface_id_type;

public:
    virtual ~Single() {}

    Single(surface_id_type const& surface_id,
           particle_id_pair const& particle)
        : base_type(surface_id), particle_(particle) {}

    particle_id_pair const& particle() const
    {
        return particle_;
    }

protected:
    const particle_id_pair particle_;
};

#endif /* SINGLE_HPP */
