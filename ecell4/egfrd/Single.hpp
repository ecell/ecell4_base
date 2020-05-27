#ifndef ECELL4_EGFRD_SINGLE_HPP
#define ECELL4_EGFRD_SINGLE_HPP

#include <utility>
#include "ShapedDomain.hpp"

namespace ecell4
{
namespace egfrd
{
template<typename Ttraits_>
class Single: public ShapedDomain<Ttraits_>
{
public:
    typedef ShapedDomain<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::traits_type::position_type position_type;
    typedef typename traits_type::world_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::domain_id_type identifier_type;
    typedef typename traits_type::world_type::traits_type::D_type D_type;

public:
    virtual ~Single() {}

    Single(identifier_type const& id,
           particle_id_pair const& particle)
        : base_type(id), particle_(particle) {}

    particle_id_pair const& particle() const
    {
        return particle_;
    }

    particle_id_pair& particle()
    {
        return particle_;
    }

    D_type const& D() const
    {
        return particle_.second.D();
    }

protected:
    particle_id_pair particle_;
};

} // egfrd
} // ecell4
#endif /* SINGLE_HPP */
