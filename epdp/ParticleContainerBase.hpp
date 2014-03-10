#ifndef PARTICLE_CONTAINER_BASE_HPP
#define PARTICLE_CONTAINER_BASE_HPP
#include "utils/range.hpp"
#include "utils/get_mapper_mf.hpp"
#include "utils/unassignable_adapter.hpp"
#include "MatrixSpace.hpp"
#include "abstract_set.hpp"
#include "generator.hpp"
#include "exceptions.hpp"
#include "ParticleContainer.hpp"
#include "Transaction.hpp"

template<typename Ttraits_>
struct ParticleContainerUtils
{
    typedef Ttraits_ traits_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef unassignable_adapter<particle_id_pair_and_distance, get_default_impl::std::vector> particle_id_pair_and_distance_list;

    struct distance_comparator:
            public std::binary_function<
                typename particle_id_pair_and_distance_list::placeholder,
                typename particle_id_pair_and_distance_list::placeholder,
                bool>
    {
        typedef typename particle_id_pair_and_distance_list::placeholder
                first_argument_type;
        typedef typename particle_id_pair_and_distance_list::const_caster const_caster;
        bool operator()(first_argument_type const& lhs,
                        first_argument_type const& rhs) const
        {
            return c_(lhs).second < c_(rhs).second;
        }

        const_caster c_;
    };

    template<typename Tset_>
    struct overlap_checker
    {
        overlap_checker(Tset_ const& ignore = Tset_()): ignore_(ignore), result_(0) {}

        template<typename Titer_>
        void operator()(Titer_ const& i, length_type const& dist)
        {
            if (!contains(ignore_, (*i).first))
            {
                if (!result_)
                {
                    result_ = new particle_id_pair_and_distance_list();
                }
                result_->push_back(std::make_pair(*i, dist));
            }
        }

        particle_id_pair_and_distance_list* result() const
        {
            if (result_)
            {
                std::sort(result_->pbegin(), result_->pend(), compare_);
            }
            return result_;
        }

    private:
        Tset_ const& ignore_;
        particle_id_pair_and_distance_list* result_;
        distance_comparator compare_;
    };
};

template<typename Tderived_, typename Ttraits_ = typename Tderived_::traits_type>
class ParticleContainerBase
    : public ParticleContainer<Ttraits_>
{
public:
    typedef ParticleContainerUtils<Ttraits_> utils;
    typedef ParticleContainer<Ttraits_> base_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator particle_id_generator;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::particle_type::shape_type particle_shape_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef Transaction<traits_type> transaction_type;

    typedef MatrixSpace<particle_type, particle_id_type, get_mapper_mf> particle_matrix_type;
    typedef abstract_limited_generator<particle_id_pair> particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    typedef sized_iterator_range<typename particle_matrix_type::const_iterator> particle_id_pair_range;

    typedef unassignable_adapter<particle_id_pair_and_distance, get_default_impl::std::vector> particle_id_pair_and_distance_list;

protected:
public:
    ParticleContainerBase(length_type world_size, size_type size)
        : pmat_(world_size, size) {}

    virtual size_type num_particles() const
    {
        return pmat_.size();
    }

    virtual length_type world_size() const
    {
        return pmat_.world_size();
    }

    length_type cell_size() const
    {
        return pmat_.cell_size();
    }

    size_type matrix_size() const
    {
        return pmat_.matrix_size();
    }

    template<typename T_>
    length_type distance(T_ const& lhs, position_type const& rhs) const
    {
        return traits_type::distance(lhs, rhs, world_size());
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return traits_type::distance(lhs, rhs, world_size());
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return traits_type::apply_boundary(v, world_size());
    }

    virtual length_type apply_boundary(length_type const& v) const
    {
        return traits_type::apply_boundary(v, world_size());
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return traits_type::cyclic_transpose(p0, p1, world_size());
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return traits_type::cyclic_transpose(p0, p1, world_size());
    }

    template<typename T1_>
    T1_ calculate_pair_CoM(
        T1_ const& p1, T1_ const& p2, 
        typename element_type_of<T1_>::type const& D1,
        typename element_type_of<T1_>::type const& D2)
    {
        typedef typename element_type_of< T1_ >::type element_type;   

        T1_ retval;

        const T1_ p2t(cyclic_transpose(p2, p1));

        return modulo(
            divide(
                add(multiply(p1, D2), multiply(p2t, D1)),
                add(D1, D2)),
            world_size());
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return check_overlap<particle_shape_type>(s);
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
    particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s, Tset_ const& ignore,
        typename boost::disable_if<boost::is_same<Tsph_, particle_id_pair> >::type* =0) const
    {
        typename utils::template overlap_checker<Tset_> oc(ignore);
        traits_type::take_neighbor(pmat_, oc, s);
        return oc.result();
    }

    template<typename Tsph_>
    particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s,
        typename boost::disable_if<boost::is_same<Tsph_, particle_id_pair> >::type* =0) const
    {
        typename utils::template overlap_checker<boost::array<particle_id_type, 0> > oc;
        traits_type::take_neighbor(pmat_, oc, s);
        return oc.result();
    }

    particle_id_pair get_particle(particle_id_type const& id, bool& found) const
    {
        typename particle_matrix_type::const_iterator i(pmat_.find(id));
        if (pmat_.end() == i) {
            found = false;
            return particle_id_pair();
        }
        found = true;
        return *i;
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        typename particle_matrix_type::const_iterator i(pmat_.find(id));
        if (pmat_.end() == i) {
            throw not_found(std::string("No such particle: id=")
                    + boost::lexical_cast<std::string>(id));
        }
        return *i;
    }

    virtual bool has_particle(particle_id_type const& id) const
    {
        return pmat_.end() != pmat_.find(id);
    }

    virtual transaction_type* create_transaction();

    virtual particle_id_pair_generator* get_particles() const
    {
        return make_range_generator<particle_id_pair>(pmat_);
    }

    particle_id_pair_range get_particles_range() const
    {
        return particle_id_pair_range(pmat_.begin(), pmat_.end(), pmat_.size());
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        return pmat_.update(pi_pair).second;
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        return pmat_.erase(id);
    }

protected:
    particle_matrix_type pmat_;
};

template<typename Tderived_, typename Ttraits_>
inline Transaction<Ttraits_>*
ParticleContainerBase<Tderived_, Ttraits_>::create_transaction()
{
    return new TransactionImpl<ParticleContainerBase>(*this);
}

#endif /* PARTICLE_CONTAINER_BASE_HPP */
