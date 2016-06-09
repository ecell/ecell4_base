
#include "ShapeContainer.hpp"
#include "collision.hpp"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

namespace ecell4
{

PlanarSurfaceContainer::PlanarSurfaceContainer(void)
{
    ;
}

Real3 
PlanarSurfaceContainer::apply_reflection(const Real3 &from, const Real3 &displacement) const
{
    std::size_t bound_surface = std::numeric_limits<std::size_t>::infinity();
    std::vector<signed int> save_isinside(this->num_surfaces()); // NEVER TO USE vector<bool>
    surface_container_type const surfaces = this->list_surfaces();

    Real3 remaining_displacement = displacement;
    Real3 current_pos = from;

    for(std::size_t i = 0; i < this->num_surfaces(); i++) {
        // To check after the apply_reflection.
        save_isinside[i] = collision::sgn(surfaces[i].second.second.is_inside(from));
    }
    do {
        bool reflection_happen = false;
        boost::tuple<bool, Real3, Real3> closest;
        std::size_t closest_surface = -1;
        for(std::size_t i = 0; i < this->num_surfaces(); i++) {
            if (i != bound_surface) {
                boost::tuple<bool, Real3, Real3> t = collision::reflect_PlanarSurface(
                        surfaces[i].second.second, from, remaining_displacement);
                if (t.get<0>() == true) {
                    if (false == reflection_happen) {
                        reflection_happen = true;
                        closest = t;
                        closest_surface = i;
                    } else if (length(from - t.get<1>() ) < length(from - closest.get<1>())) {
                        // update
                        closest = t;
                        closest_surface = i;
                    }
                }
            }
        }
        if (reflection_happen == false) {
            break;
        } else {
            current_pos = closest.get<1>();
            remaining_displacement = closest.get<2>();
            bound_surface = closest_surface;
            //std::cout << "bound_surface (" << bound_surface << ") at " << from << std::endl;
        }
    } while(true);
    Real3 final_pos = current_pos + remaining_displacement;
    for(std::size_t i = 0; i < this->num_surfaces(); i++) {
        if( save_isinside[i] != collision::sgn(surfaces[i].second.second.is_inside(final_pos)) ){
            throw IllegalState("Particle moved to beyond the PlanarSurface");
        }
    }
    return final_pos;
}

}


