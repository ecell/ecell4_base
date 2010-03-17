#ifndef STRUCTURE_UTILS_HPP
#define STRUCTURE_UTILS_HPP

#include <string>
#include <typeinfo>
#include "linear_algebra.hpp"
#include "geometry.hpp"
#include "exceptions.hpp"
#include "Surface.hpp"
#include "Region.hpp"

template<typename Ttraits_>
struct StructureUtils
{
    typedef Ttraits_ traits_type;
    typedef typename traits_type::world_type::position_type position_type;
    typedef typename traits_type::world_type::length_type length_type;
    typedef typename traits_type::world_type::surface_id_type surface_id_type;
    typedef typename traits_type::world_type::surface_type surface_type;
    typedef typename traits_type::sphere_type sphere_type;
    typedef typename traits_type::cylinder_type cylinder_type;
    typedef typename traits_type::box_type box_type;
    typedef typename traits_type::spherical_surface_type spherical_surface_type;
    typedef typename traits_type::cylindrical_surface_type cylindrical_surface_type;
    typedef typename traits_type::planar_surface_type planar_surface_type;
    typedef typename traits_type::cuboidal_region_type cuboidal_region_type;
    typedef typename traits_type::rng_type rng_type;

    static planar_surface_type* create_planar_surface(
            surface_id_type const& id,
            position_type const& pos,
            position_type const& vector_x,
            position_type const& vector_y,
            length_type const& lx,
            length_type const& ly,
            length_type const& lz = 0)
    {
        BOOST_ASSERT(dot_product(vector_x, vector_y) == 0.);
        const position_type unit_x(normalize(vector_x));
        const position_type unit_y(normalize(vector_y));
        const position_type unit_z(cross_product(unit_x, unit_y));

        return new planar_surface_type(id,
                                       box_type(pos, unit_x, unit_y, unit_z,
                                                lx, ly, lz));
    }

    static spherical_surface_type* create_spherical_surface(
            surface_id_type const& id,
            position_type const& pos,
            length_type const& radius)
    {
        return new spherical_surface_type(id, sphere_type(pos, radius));
    }

    static cylindrical_surface_type* create_cylindrical_surface(
            surface_id_type const& id,
            position_type const& pos,
            length_type const& radius,
            position_type const& unit_z,
            length_type const& size)
    {
        return new cylindrical_surface_type(id,
                cylinder_type(pos, radius, unit_z, size));
    }


    static cuboidal_region_type* create_cuboidal_region(
            surface_id_type const& id,
            position_type const& corner,
            boost::array<length_type, 3> const& extent)
    {
        const boost::array<length_type, 3> half_of_extent(divide(extent, 2));
        return new cuboidal_region_type(id,
                box_type(add(corner, half_of_extent),
                         create_vector<position_type>(1, 0, 0),
                         create_vector<position_type>(0, 1, 0),
                         create_vector<position_type>(0, 0, 1),
                         half_of_extent));
    }

    static position_type draw_bd_displacement(planar_surface_type const& surface, length_type const& r, rng_type& rng)
    {
        length_type const x(rng.normal(0., r)), y(rng.normal(0., r));
        return add(
            multiply(surface.shape().unit_x(), x),
            multiply(surface.shape().unit_y(), y));
    }

    static position_type draw_bd_displacement(cylindrical_surface_type const& surface, length_type const& r, rng_type& rng)
    {
        return multiply(surface.shape().unit_z(), rng.normal(0., r));
    }

    static position_type draw_bd_displacement(cuboidal_region_type const& region, length_type const& r, rng_type& rng)
    {
        return create_vector<position_type>(
            rng.normal(0., r),
            rng.normal(0., r),
            rng.normal(0., r));
    }

    static position_type draw_bd_displacement(surface_type const& surface, length_type const& r, rng_type& rng)
    {
        {
            planar_surface_type const* _surface(dynamic_cast<planar_surface_type const*>(&surface));
            if (_surface)
            {
                return draw_bd_displacement(*_surface, r, rng);
            }
        }
        {
            cylindrical_surface_type const* _surface(dynamic_cast<cylindrical_surface_type const*>(&surface));
            if (_surface)
            {
                return draw_bd_displacement(*_surface, r, rng);
            }
        }
        {
            cuboidal_region_type const* _surface(dynamic_cast<cuboidal_region_type const*>(&surface));
            if (_surface)
            {
                return draw_bd_displacement(*_surface, r, rng);
            }
        }

        throw unsupported(std::string("Unsupported structure type: ") + typeid(surface).name());
    }
};

#endif /* STRUCTURE_UTILS_HPP */
