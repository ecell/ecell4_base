#ifndef PARTICLE_CONTAINER_HPP
#define PARTICLE_CONTAINER_HPP

#include <utility>
#include <boost/shared_ptr.hpp>
#include "generator.hpp"
#include "utils/get_default_impl.hpp"
#include "utils/unassignable_adapter.hpp"

#include <ecell4/core/types.hpp>
#include <ecell4/core/Space.hpp>

template<typename Ttraits_>
class Transaction;

template<typename Ttraits_>
class ParticleContainer
    : ecell4::Space
{
public:

    typedef Ttraits_ traits_type;

    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename traits_type::species_type species_type;
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

public:

    virtual ~ParticleContainer() {};

    virtual ecell4::Integer num_particles() const = 0;
    // virtual size_type num_particles() const = 0;

    virtual length_type world_size() const = 0;

    virtual species_type const& get_species(species_id_type const& id) const = 0;

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

    virtual particle_id_pair_generator* get_particles() const = 0;

    virtual transaction_type* create_transaction() = 0;

    virtual length_type distance(
        position_type const& lhs, position_type const& rhs) const = 0;

    virtual position_type apply_boundary(position_type const& v) const = 0;

    virtual length_type apply_boundary(length_type const& v) const = 0;

    virtual position_type cyclic_transpose(
        position_type const& p0, position_type const& p1) const = 0;

    virtual length_type cyclic_transpose(
        length_type const& p0, length_type const& p1) const = 0;
};


#endif /* PARTICLE_CONTAINER_HPP */
