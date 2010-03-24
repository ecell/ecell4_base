#ifndef MULTI_HPP
#define MULTI_HPP

#include "exceptions.hpp"
#include "Domain.hpp"
#include "ParticleContainer.hpp"
#include "Sphere.hpp"
#include "utils/range.hpp"

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
    typedef typename traits_type::surface_id_type surface_id_type;
    typedef typename traits_type::surface_type surface_type;
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

    Multi(domain_id_type const& id, surface_id_type const& surface_id,
          simulator_type& main)
        : base_type(id, surface_id),
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
