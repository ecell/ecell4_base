#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include "binding_common.hpp"

namespace binding {

static Position 
calculate_pair_CoM(Position const& p1, 
                   Position const& p2, 
                   element_type_of< Position >::type const& D1,
                   element_type_of< Position >::type const& D2,
                   element_type_of< Position >::type const& world_size)
{
    typedef element_type_of<Position>::type element_type;   

    Position retval;

    const Position p2t(cyclic_transpose<Position>(p2, p1, world_size));

    return modulo(
        divide(
            add(multiply(p1, D2), multiply(p2t, D1)),
            add(D1, D2)),
        world_size);
}

Position draw_bd_displacement(Structure const& surface, Length const& length, GSLRandomNumberGenerator& rng)
{
    return StructureUtils::draw_bd_displacement(surface, length, rng);
}

void register_module_functions()
{
    using namespace boost::python;
    def( "length_sq", &length_sq<Position> );
    def( "length", &length<Position> );
    def( "distance", (WorldTraits::length_type(*)(Position const&, Position const&))&distance<Position> );
    def( "distance_cyclic", &distance_cyclic<Position, Position> );
    def( "apply_boundary", &apply_boundary<Position, WorldTraits::length_type> );
    def( "calculate_pair_CoM", &calculate_pair_CoM );

    def( "normalize", (Position(*)(Position const&))&normalize<Position> );
    def( "normalize", (Position(*)(Position const&, WorldTraits::length_type const&))&normalize<Position> );
    def( "cyclic_transpose", &cyclic_transpose<Position, element_type_of<Position>::type> );
    def("create_planar_surface", &StructureUtils::create_planar_surface,
            return_value_policy<manage_new_object>());
    def("create_spherical_surface", &StructureUtils::create_spherical_surface,
            return_value_policy<manage_new_object>());
    def("create_cylindrical_surface", &StructureUtils::create_cylindrical_surface,
            return_value_policy<manage_new_object>());
    def("create_cuboidal_region", &StructureUtils::create_cuboidal_region,
            return_value_policy<manage_new_object>());
    def("draw_bd_displacement", &draw_bd_displacement);
}

} // namespace binding
