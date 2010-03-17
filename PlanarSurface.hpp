#ifndef PLANAR_SURFACE_HPP
#define PLANAR_SURFACE_HPP

#include <string>
#include "Shape.hpp"
#include "Surface.hpp"
#include "DomainFactory.hpp"

template<typename Tid_, typename Tbox_>
class PlanarSurface: public Surface<Tid_,
                                    typename shape_position_type<Tbox_>::type,
                                    typename shape_length_type<Tbox_>::type>
{
public:
    typedef Surface<Tid_,
                    typename shape_position_type<Tbox_>::type,
                    typename shape_length_type<Tbox_>::type> base_type;
    typedef Tbox_ shape_type;
    typedef Tid_ identifier_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::length_type length_type;

public:
    PlanarSurface(identifier_type const& id, position_type const& pos,
                  position_type const& vector_x, position_type const& vector_y)
        : base_type(id), shape_()
    {
        // BOOST_ASSERT(dot_product(vector_x, vector_y) == 0.);
        shape_position(shape_) = pos;
        shape_.unit_x() = normalize(vector_x);
        shape_.unit_y() = normalize(vector_y);
        shape_.unit_z() = cross_product(vector_x, vector_y);
    }

    shape_type const& shape() const
    {
        return shape_;
    }

    shape_type& shape()
    {
        return shape_;
    }

protected:
    shape_type shape_;
};

template<typename Tsurface_>
struct PlanarSurfaceDomainFactory: public DomainFactory
{
    virtual single_type*
    create_single(domain_id_type const& domain_id,
                  particle_id_pair const& particle_id_pair,
                  shell_id_pair const& shell_id_pair,
                  typename network_rules_type::reaction_rule_vector const& reaction_types) const
    {
        // TBD
        return 0;
    }

    virtual pair_type* create_pair(domain_id_type const& domain_id,
                                   position_type const& com,
                                   particle_id_pair const& single1,
                                   particle_id_pair const& single2,
                                   shell_id_pair const& shell_id_pair,
                                   typename network_rules_type::reaction_rule_vector const& reactions) const
    {
        // TBD
        return 0;
    }
};

#endif /* PLANAR_SURFACE_HPP */
