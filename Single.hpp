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
    typedef typename traits_type::world_type::structure_id_type structure_id_type;

public:
    virtual ~Single() {}

    Single(structure_id_type const& structure_id,
           particle_id_pair const& particle)
        : base_type(structure_id), particle_(particle) {}

    particle_id_pair const& particle() const
    {
        return particle_;
    }

protected:
    const particle_id_pair particle_;
};

#endif /* SINGLE_HPP */
