#ifndef TRANSACTION_HPP
#define TRANSACTION_HPP

#include <vector>
#include <map>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include "utils.hpp"
#include "exceptions.hpp"
#include "ParticleContainer.hpp"
#include "sorted_list.hpp"
#include "generator.hpp"
#include "utils/unassignable_adapter.hpp"
#include "utils/stringizer.hpp"

template<typename Ttraits_>
class Transaction: public ParticleContainer<Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef typename traits_type::particle_id_pair_generator
        particle_id_pair_generator;

    virtual ~Transaction() {}

    virtual particle_id_pair_generator* get_added_particles() const = 0;

    virtual particle_id_pair_generator* get_removed_particles() const = 0;

    virtual particle_id_pair_generator* get_modified_particles() const = 0;

    virtual void rollback() = 0;
};

template<typename Tpc_>
class TransactionImpl: public Transaction<typename Tpc_::traits_type>
{
public:
    typedef Tpc_ particle_container_type;
    typedef typename particle_container_type::traits_type traits_type;
    typedef Transaction<traits_type> base_type;

    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename traits_type::molecule_info_type molecule_info_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::particle_id_pair_and_distance
        particle_id_pair_and_distance;
    typedef typename traits_type::particle_id_pair_and_distance_list
        particle_id_pair_and_distance_list;
    typedef typename traits_type::particle_id_pair_generator
        particle_id_pair_generator;

    typedef typename base_type::time_type time_type;

private:
    typedef std::map<typename particle_id_pair::first_type,
            typename particle_id_pair::second_type> particle_id_pair_set_type;
    typedef sorted_list<std::vector<particle_id_type> > particle_id_list_type;

public:
    virtual std::pair<particle_id_pair, bool> new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        std::pair<particle_id_pair, bool> retval(pc_.new_particle(sid, pos));
        const bool result(added_particles_.push_no_duplicate(retval.first.first));
        BOOST_ASSERT(result);
        return retval;
    }

    virtual bool update_particle(const particle_id_type& pid, const particle_type& p)
    {
        BOOST_ASSERT(removed_particles_.end() ==
                removed_particles_.find(pid));
        std::pair<typename particle_id_pair_set_type::iterator, bool> r(
                orig_particles_.insert(particle_id_pair(
                    pid, particle_type())));
        if (r.second &&
            added_particles_.end() == added_particles_.find(pid))
        {
            modified_particles_.push_no_duplicate(pid);
            particle_type _v(pc_.get_particle(pid).second);
            std::swap((*r.first).second, _v);
        }
        return pc_.update_particle(pid, p);
    }

    virtual void remove_particle(particle_id_type const& id)
    {
        std::pair<typename particle_id_pair_set_type::iterator, bool> r(
                orig_particles_.insert(particle_id_pair(
                    id, particle_type())));
        if (r.second)
        {
            particle_type _v(pc_.get_particle(id).second);
            std::swap((*r.first).second, _v);
        }

        if (added_particles_.erase(id) == 0)
        {
            modified_particles_.erase(id);
            const bool result(removed_particles_.push_no_duplicate(id));
            BOOST_ASSERT(result);
        }
        else
        {
            orig_particles_.erase(id);
        }

        pc_.remove_particle(id);
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        return pc_.get_particle(id);
    }

    virtual bool has_particle(particle_id_type const& id) const
    {
        return pc_.has_particle(id);
    }

    virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s) const
    {
        return pc_.check_overlap(s);
    }

    virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return pc_.check_overlap(s, ignore);
    }

    virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return pc_.check_overlap(s, ignore1, ignore2);
    }

    virtual Transaction<traits_type>* create_transaction()
    {
        return new TransactionImpl<particle_container_type>(*this);
    }

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return pc_.get_structure(id);
    }

    // virtual molecule_info_type const& find_molecule_info(species_id_type const& id) const
    // {
    //     return pc_.find_molecule_info(id);
    // }

    // virtual molecule_info_type const& get_molecule_info(species_id_type const& id)
    // {
    //     return pc_.get_molecule_info(id);
    // }

    virtual molecule_info_type get_molecule_info(species_id_type const& id) const
    {
        return pc_.get_molecule_info(id);
    }

    virtual ecell4::Integer num_particles() const
    {
        return pc_.num_particles();
    }

    // virtual size_type num_particles() const
    // {
    //     return pc_.num_particles();
    // }

    virtual const position_type& edge_lengths() const
    {
        return pc_.edge_lengths();
    }

    virtual particle_id_pair_generator* get_added_particles() const
    {
        return make_range_generator<true>(
            make_transform_iterator_range(added_particles_,
                boost::bind(&TransactionImpl::get_particle, this, _1)));
    }

    virtual particle_id_pair_generator* get_removed_particles() const
    {
        return make_range_generator<true>(
            make_transform_iterator_range(removed_particles_,
                boost::bind(&TransactionImpl::get_original_particle, this, _1)));
    }

    virtual particle_id_pair_generator* get_modified_particles() const
    {
        return make_range_generator<true>(
            make_transform_iterator_range(modified_particles_,
                boost::bind(&TransactionImpl::get_particle, this, _1)));
    }

    virtual void rollback()
    {
        for (typename particle_id_pair_set_type::iterator
                i(orig_particles_.begin()), e(orig_particles_.end());
                i != e; ++i)
        {
            pc_.update_particle((*i).first, (*i).second);
        }

        for (typename particle_id_list_type::iterator
                i(added_particles_.begin()), e(added_particles_.end());
                i != e; ++i)
        {
            pc_.remove_particle(*i);
        }
        added_particles_.clear();
        modified_particles_.clear();
        removed_particles_.clear();
        orig_particles_.clear();
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return pc_.distance(lhs, rhs);
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return pc_.apply_boundary(v);
    }

    // virtual length_type apply_boundary(length_type const& v) const
    // {
    //     return pc_.apply_boundary(v);
    // }

    virtual position_type periodic_transpose(position_type const& p0, position_type const& p1) const
    {
        return pc_.periodic_transpose(p0, p1);
    }

    // virtual length_type periodic_transpose(length_type const& p0, length_type const& p1) const
    // {
    //     return pc_.periodic_transpose(p0, p1);
    // }

    virtual ~TransactionImpl() {}

    TransactionImpl(particle_container_type& pc): pc_(pc) {}

    /** ecell4::Space
     */
    // virtual const time_type& t() const
    virtual const time_type t() const
    {
        return pc_.t();
    }

    virtual void set_t(const time_type& t)
    {
        pc_.set_t(t);
    }

private:
    particle_id_pair get_original_particle(particle_id_type const& id) const
    {
        typename particle_id_pair_set_type::const_iterator i(orig_particles_.find(id));
        if (orig_particles_.end() == i)
        {
            throw not_found(std::string("No such particle: id=")
                    + boost::lexical_cast<std::string>(id));
        }
        return *i;
    }

private:
    particle_container_type& pc_;
    particle_id_list_type added_particles_;
    particle_id_list_type modified_particles_;
    particle_id_pair_set_type orig_particles_;
    particle_id_list_type removed_particles_;
};

#endif /* TRANSACTION_HPP */
