#ifndef PARTICLE_CONTAINER_HPP
#define PARTICLE_CONTAINER_HPP

#include <utility>
#include <boost/shared_ptr.hpp>
#include "generator.hpp"
#include "utils/get_default_impl.hpp"
#include "utils/unassignable_adapter.hpp"

#include <ecell4/core/types.hpp>
#include <ecell4/core/Space.hpp>
#include <ecell4/core/ShapeContainer.hpp>

template<typename Ttraits_>
class Transaction;

template<typename Ttraits_>
class ParticleContainer
    : public ecell4::Space
{
public:

    typedef Ttraits_ traits_type;

    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename traits_type::molecule_info_type molecule_info_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::particle_id_pair_generator
        particle_id_pair_generator;
    typedef typename traits_type::particle_id_pair_and_distance
        particle_id_pair_and_distance;
    typedef typename traits_type::particle_id_pair_and_distance_list
        particle_id_pair_and_distance_list;
    typedef Transaction<traits_type> transaction_type;

    // Surface-related type definitions
    typedef ecell4::PlanarSurfaceContainer surface_container_type;
    typedef surface_container_type::surface_type surface_type;
    typedef surface_container_type::surface_id_type surface_id_type;

public:

    virtual ~ParticleContainer() {};

    virtual ecell4::Integer num_particles() const = 0;
    // virtual size_type num_particles() const = 0;

    virtual molecule_info_type const& get_molecule_info(species_id_type const& id) = 0;
    virtual molecule_info_type const& find_molecule_info(species_id_type const& id) const = 0;

    virtual boost::shared_ptr<structure_type> get_structure(
        structure_id_type const& id) const = 0;

    virtual particle_id_pair new_particle(
        species_id_type const& sid, position_type const& pos) = 0;

    virtual bool update_particle(particle_id_pair const& pi_pair) = 0;

    virtual bool remove_particle(particle_id_type const& id) = 0;

    virtual particle_id_pair get_particle(particle_id_type const& id) const = 0;

    virtual bool has_particle(particle_id_type const& id) const = 0;

    virtual particle_id_pair_and_distance_list* check_overlap(
        particle_shape_type const& s) const = 0;

    virtual particle_id_pair_and_distance_list* check_overlap(
        particle_shape_type const& s, particle_id_type const& ignore) const = 0;

    virtual particle_id_pair_and_distance_list* check_overlap(
        particle_shape_type const& s, particle_id_type const& ignore1,
        particle_id_type const& ignore2) const = 0;

    virtual bool no_overlap(particle_shape_type const& s) const
    {
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            check_overlap(s));
        // boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
        //     check_overlap<particle_shape_type>(s));
        return (!overlapped || ::size(*overlapped) == 0);
    }

    virtual bool no_overlap(particle_shape_type const& s,
        particle_id_type const& ignore) const
    {
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            check_overlap(s, ignore));
        return (!overlapped || ::size(*overlapped) == 0);
    }

    virtual bool no_overlap(particle_shape_type const& s,
        particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            check_overlap(s, ignore1, ignore2));
        return (!overlapped || ::size(*overlapped) == 0);
    }

    virtual particle_id_pair_generator* get_particles() const = 0;

    virtual transaction_type* create_transaction() = 0;

    virtual length_type distance(
        position_type const& lhs, position_type const& rhs) const = 0;

    virtual position_type apply_boundary(position_type const& v) const = 0;

    // virtual length_type apply_boundary(length_type const& v) const = 0;

    virtual position_type cyclic_transpose(
        position_type const& p0, position_type const& p1) const = 0;

    // virtual length_type cyclic_transpose(
    //     length_type const& p0, length_type const& p1) const = 0;

    virtual void save(const std::string& filename) const
    {
        throw ecell4::NotSupported(
            "save(const std::string) is not supported by this space class");
    }

    std::pair<std::pair<surface_id_type, surface_type>, bool>
    //new_surface(const Species &sp, const surface_type &surface)
    new_surface(const ecell4::Species &sp, const surface_type &surface)
    {
        ecell4::PlanarSurfaceID psid(psidgen_());
        surfaces_.update_surface(psid, sp, surface);
        return std::make_pair( std::make_pair(psid, surface), true);
    }
    Integer num_surfaces(void) const
    {
        return surfaces_.num_surfaces();
    }

    position_type apply_reflection(position_type const &from, position_type const &displacement) const 
    {
        return this->surfaces_.apply_reflection(from, displacement);
    }
    
private:
    surface_container_type surfaces_;
    ecell4::SerialIDGenerator<ecell4::PlanarSurfaceID> psidgen_;
};


#endif /* PARTICLE_CONTAINER_HPP */
