#ifndef PLANAR_SURFACE_HPP
#define PLANAR_SURFACE_HPP

#include <string>
#include "Surface.hpp"

template<typename Tsimulator_>
class PlanarSurface: public Surface<Tsimulator_>
{
public:
    typedef Surface<Tsimulator_> base_type;
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;
    typedef typename base_type::D_type D_type;
    typedef typename base_type::particle_id_pair particle_id_pair;
    typedef typename base_type::shell_id_pair shell_id_pair;
    typedef typename base_type::time_type time_type;
    typedef typename base_type::domain_type domain_type;
    typedef typename base_type::single_type single_type;
    typedef typename base_type::pair_type pair_type;
    typedef typename base_type::network_rules_type network_rules_type;
    typedef Box<length_type> shape_type;

public:
    PlanarSurface(simulator_type* sim, std::string const& name,
                  shape_type const& shape)
        : base_type(sim, name)
    {
        BOOST_ASSERT(dot_product(shape.vector_x(), shape.vector_y()) == 0.);
        shape_ = shape;
        shape_.vector_z() = cross_product(shape.vector_x(), shape.vector_y());

    }

    shape_type const& shape() const
    {
        return shape_;
    }

    shape_type& shape()
    {
        return shape_;
    }

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

    virtual length_type distance(position_type const& pos) const
    {
        return distance(shape_, pos);
    }

    virtual position_type draw_bd_displacement(time_type const& dt, D_type const& D)
    {
        std::sqrt(2. * D * dt);
        sim_->get_rng().normal(0., r)
    }

protected:
    shape_type shape_;
}

#endif /* PLANAR_SURFACE_HPP */
