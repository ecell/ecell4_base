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
    : public ecell4::Space
{
public:

    typedef Ttraits_ traits_type;
    typedef Transaction<traits_type> transaction_type;

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

public:

    virtual ~ParticleContainer() {};

    virtual molecule_info_type get_molecule_info(species_id_type const& id) const = 0;
    // virtual molecule_info_type const& get_molecule_info(species_id_type const& id) = 0;
    // virtual molecule_info_type const& find_molecule_info(species_id_type const& id) const = 0;

    virtual boost::shared_ptr<structure_type> get_structure(
        structure_id_type const& id) const = 0;

    virtual std::pair<particle_id_pair, bool> new_particle(
        species_id_type const& sid, position_type const& pos) = 0;

    virtual particle_id_pair_and_distance_list check_overlap(
        particle_shape_type const& s) const = 0;

    virtual particle_id_pair_and_distance_list check_overlap(
        particle_shape_type const& s, particle_id_type const& ignore) const = 0;

    virtual particle_id_pair_and_distance_list check_overlap(
        particle_shape_type const& s, particle_id_type const& ignore1,
        particle_id_type const& ignore2) const = 0;

    virtual bool no_overlap(particle_shape_type const& s) const
    {
        const particle_id_pair_and_distance_list overlapped(
            check_overlap(s));
        return (overlapped.size() == 0);
    }

    virtual bool no_overlap(particle_shape_type const& s,
        particle_id_type const& ignore) const
    {
        const particle_id_pair_and_distance_list overlapped(
            check_overlap(s, ignore));
        return (overlapped.size() == 0);
    }

    virtual bool no_overlap(particle_shape_type const& s,
        particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        const particle_id_pair_and_distance_list overlapped(
            check_overlap(s, ignore1, ignore2));
        return (overlapped.size() == 0);
    }

    virtual transaction_type* create_transaction() = 0;

    /**
     * ParticleSpace has functions similar, but a bit different.
     */

    virtual bool update_particle(const particle_id_type& pid, const particle_type& p) = 0;

    virtual position_type periodic_transpose(
        position_type const& p0, position_type const& p1) const = 0;

    /**
     * Space
     */

    virtual void save(const std::string& filename) const
    {
        throw ecell4::NotSupported(
            "save(const std::string) is not supported by this space class");
    }

    /**
     * ParticleSpace
     */

    virtual void remove_particle(particle_id_type const& id) = 0;

    virtual length_type distance(
        position_type const& lhs, position_type const& rhs) const = 0;

    virtual position_type apply_boundary(position_type const& v) const = 0;

    virtual void add_surface(const boost::array<position_type, 3>& vertices)
    {
        return;
    }

    virtual position_type apply_reflection(const position_type& pos, const position_type& disp)
    {
        return pos + disp;
    }

    virtual position_type apply_structure(const position_type& pos, const position_type& disp)
    {
        return pos + disp;
    }
};


#endif /* PARTICLE_CONTAINER_HPP */
