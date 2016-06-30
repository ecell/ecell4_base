
#include "ShapeContainer.hpp"
#include "collision.hpp"
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>


#include <boost/format.hpp>

namespace ecell4
{

class PlanarSurfaceReflectionLogger
{
public:
    typedef boost::tuple<bool, Real3, Real3> reflection_discriptor_type;
    typedef std::vector< std::pair<std::size_t, reflection_discriptor_type> > log_container_type;
public:
    PlanarSurfaceReflectionLogger(const Real3 &from, const Real3 &displacement)
        :from_(from), displacement_(displacement)
    {;}
    void add(const std::size_t surface_index, const reflection_discriptor_type discriptor)
    {
        std::pair<std::size_t, reflection_discriptor_type> a(
                std::make_pair(surface_index, discriptor));
        this->log_container_.push_back(a);
    }
    void to_str() const
    {
        std::cerr << "From: " << this->from_ << std::endl;
        std::cerr << "Move: " << this->displacement_ << std::endl;
        for(log_container_type::const_iterator it = this->log_container_.begin(); it != this->log_container_.end(); it++)
        {
            std::cerr << "surface: " << it->first << "  pos: " << it->second.get<1>() << std::endl;
        }
        std::cerr << "Bound log output done" << std::endl;
    }
    void clear()
    {   this->log_container_.clear();   }
private:
    log_container_type log_container_;
    Real3 const &from_;
    Real3 const &displacement_;
};

PlanarSurfaceContainer::PlanarSurfaceContainer(void)
{
    ;
}

Real3 
PlanarSurfaceContainer::apply_reflection(const Real3 &from, const Real3 &displacement) const
{
    std::vector<Real> save_isinside(this->num_surfaces()); // NEVER TO USE vector<bool>
    surface_container_type const surfaces = this->list_surfaces();

    Real3 current_pos = from;
    Real3 remaining_displacement = displacement;
    std::size_t bound_surface = this->num_surfaces() + 100;

    for(std::size_t i = 0; i < this->num_surfaces(); i++) {
        // To check after the apply_reflection.
        save_isinside[i] = surfaces[i].second.second.is_inside(from);
    }
    do {
        bool reflection_happen = false;
        boost::tuple<bool, Real3, Real3> closest;
        std::size_t closest_surface = -1;
        for(std::size_t i = 0; i < this->num_surfaces(); i++) {
            if (i != bound_surface) {
                boost::tuple<bool, Real3, Real3> t = collision::reflect_PlanarSurface(
                        surfaces[i].second.second, current_pos, remaining_displacement);
                if (t.get<0>() == true) {
                    if (false == reflection_happen) {
                        reflection_happen = true;
                        closest = t;
                        closest_surface = i;
                    } else {
                        if (length(current_pos - t.get<1>() ) < length(current_pos - closest.get<1>())) {
                        // update
                            closest = t;
                            closest_surface = i;
                        }
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
        }
    } while(true);
    Real3 final_pos = current_pos + remaining_displacement;
    for(std::size_t i = 0; i < this->num_surfaces(); i++) {
        if( save_isinside[i] * surfaces[i].second.second.is_inside(final_pos) < 0 ){
            std::cerr << "surface " << i << std::endl;
            throw IllegalState("Particle moved to beyond the PlanarSurface");
        }
    }
    return final_pos;
}

}


