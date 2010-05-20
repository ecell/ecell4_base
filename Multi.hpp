#ifndef MULTI_HPP
#define MULTI_HPP

#include "exceptions.hpp"
#include "Domain.hpp"
#include "ParticleContainer.hpp"
#include "ParticleContainerBase.hpp"
#include "Sphere.hpp"
#include "utils/array_helper.hpp"
#include "utils/range.hpp"

template<typename Ttraits_>
class MultiParticleContainer
    : public Ttraits_::world_type::particle_container_type
{
public:
    typedef ParticleContainerUtils<typename Ttraits_::world_type::traits_type> utils;
    typedef typename Ttraits_::world_type world_type;
    typedef typename world_type::traits_type traits_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename particle_type::shape_type particle_shape_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef Transaction<traits_type> transaction_type;
    typedef abstract_limited_generator<particle_id_pair> particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef unassignable_adapter<particle_id_pair_and_distance, get_default_impl::std::vector> particle_id_pair_and_distance_list;
    typedef std::map<particle_id_type, particle_type> particle_map;

    virtual ~MultiParticleContainer() {}

    virtual size_type num_particles() const
    {
        return particles_.size();
    }

    virtual length_type world_size() const
    {
        return world_.world_size();
    }

    virtual species_type const& get_species(species_id_type const& id) const
    {
        return world_.get_species(id);
    }

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return world_.get_structure(id);
    }

    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        particle_id_pair const retval(world_.new_particle(sid, pos));
        particles_.insert(retval);
        return retval;
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        bool const retval(world_.update_particle(pi_pair));
        particles_[pi_pair.first] = pi_pair.second;
        return retval;
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        world_.remove_particle(id);
        return particles_.erase(id);
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        typename particle_map::const_iterator i(particles_.find(id));
        if (particles_.end() == i)
        {
            throw not_found(std::string("No such particle: id=")
                    + boost::lexical_cast<std::string>(id));
        }
        return *i;
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return check_overlap(s, array_gen<particle_id_type>());
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return check_overlap(s, array_gen(ignore));
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return check_overlap(s, array_gen(ignore1, ignore2));
    }

    template<typename Tsph_, typename Tset_>
    particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s, Tset_ const& ignore) const
    {
        typename utils::template overlap_checker<Tset_> checker(ignore);
        for (typename particle_map::const_iterator i(particles_.begin()),
                                                   e(particles_.end());
             i != e; ++i)
        {
            length_type const dist(world_.distance(shape((*i).second), s.position()));
            if (dist < s.radius())
            {
                checker(i, dist);
            }
        }
        return checker.result();
    }

    virtual particle_id_pair_generator* get_particles() const
    {
        return make_range_generator<particle_id_pair>(particles_);
    }

    virtual transaction_type* create_transaction()
    {
        return new TransactionImpl<MultiParticleContainer>(*this);
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return world_.distance(lhs, rhs);
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return world_.apply_boundary(v);
    }

    virtual length_type apply_boundary(length_type const& v) const
    {
        return world_.apply_boundary(v);
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return world_.cyclic_transpose(p0, p1);
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return world_.cyclic_transpose(p0, p1);
    }

    MultiParticleContainer(world_type& world): world_(world) {}

private:
    world_type& world_;
    particle_map particles_;
};

template<typename Tsim_>
class Multi: public Domain<typename Tsim_::traits_type>
{
public:
    typedef Tsim_ simulator_type;
    typedef typename simulator_type::traits_type traits_type;
    typedef Domain<traits_type> base_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename particle_type::shape_type particle_shape_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::spherical_shell_type spherical_shell_type;
    typedef typename traits_type::spherical_shell_id_pair spherical_shell_id_pair;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef Transaction<traits_type> transaction_type;
    typedef abstract_limited_generator<particle_id_pair> particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef unassignable_adapter<particle_id_pair_and_distance, get_default_impl::std::vector> particle_id_pair_and_distance_list;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::domain_id_type domain_id_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::sphere_type sphere_type;
    typedef MatrixSpace<particle_type, particle_id_type, get_mapper_mf> particle_matrix_type;
    typedef std::map<shell_id_type, particle_id_type> shell_particle_id_map;
    typedef sized_iterator_range<typename shell_particle_id_map::const_iterator> shell_particle_id_pair_range;

public:
    virtual ~Multi() {}

    Multi(domain_id_type const& id, structure_id_type const& structure_id,
          simulator_type& main)
        : base_type(id, structure_id),
          main_(main), shell_ids_(), escaped_(false) {}

    spherical_shell_id_pair
    add_particle_and_shell(domain_id_type const& id,
                           particle_id_pair const& pp,
                           length_type const& shell_size)
    {
        BOOST_ASSERT(main_.world().update_particle(pp));
        spherical_shell_id_pair ssp(
            main_.new_spherical_shell(
                id, sphere_type(pp.second.position(), shell_size)));
        shell_ids_[ssp.first] = pp.first;
        return ssp;
    }

    shell_particle_id_pair_range get_shell_particle_pairs() const
    {
        return shell_particle_id_pair_range(shell_ids_);
    }

    template<typename Tset_>
    void clear_volume(particle_shape_type const& sphere, Tset_ const& ignore)
    {
        // check if the particle has escaped
        if (!within_shell(sphere))
        {
            escaped_ = true;
            clear_outer_volume(sphere, ignore);
        }
    }

protected:
    bool within_shell(particle_shape_type const& sphere) const
    {
        for (typename shell_particle_id_map::const_iterator
                i(shell_ids_.begin()), e(shell_ids_.end());
             i != e; ++i)
        {
            spherical_shell_type shell(main_.get_spherical_shell((*i).first));
            position_type ppos(main_.world().cyclic_transpose(sphere.position(), shell.position()));
            if (distance(ppos, shell.position()) < shell.radius())
            {
                return true;
            }
        }
        return false;
    }

    template<typename Tset_>
    void clear_outer_volume(particle_shape_type const& sphere, Tset_ const& ignore)
    {
        main_.clear_volume(sphere, ignore);
        if (std::auto_ptr<particle_id_pair_and_distance_list>(main_.world().check_overlap(sphere, ignore)).get())
        {
            throw no_space(__FUNCTION__);
        }
    }

protected:
    simulator_type& main_;
    shell_particle_id_map shell_ids_;
    bool escaped_;
};

#endif /* MULTI_HPP */
