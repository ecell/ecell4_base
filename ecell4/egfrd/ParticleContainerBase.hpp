#ifndef PARTICLE_CONTAINER_BASE_HPP
#define PARTICLE_CONTAINER_BASE_HPP

#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/AABBSurface.hpp>
#include "utils/range.hpp"
#include "utils/unassignable_adapter.hpp"
#include "MatrixSpace.hpp"
#include "abstract_set.hpp"
#include "generator.hpp"
#include "exceptions.hpp"
#include "ParticleContainer.hpp"
#include "Transaction.hpp"
#include "Polygon.hpp"

template<typename Ttraits_>
struct ParticleContainerUtils
{
    typedef Ttraits_ traits_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::particle_id_pair_and_distance
        particle_id_pair_and_distance;
    typedef typename traits_type::particle_id_pair_and_distance_list
        particle_id_pair_and_distance_list;

    // struct distance_comparator:
    //         public std::binary_function<
    //             typename particle_id_pair_and_distance_list::placeholder,
    //             typename particle_id_pair_and_distance_list::placeholder,
    //             bool>
    // {
    //     typedef typename particle_id_pair_and_distance_list::placeholder
    //             first_argument_type;
    //     typedef typename particle_id_pair_and_distance_list::const_caster const_caster;
    //     bool operator()(first_argument_type const& lhs,
    //                     first_argument_type const& rhs) const
    //     {
    //         return c_(lhs).second < c_(rhs).second;
    //     }
    //     const_caster c_;
    // };

    struct distance_comparator
        : public std::binary_function<
            typename particle_id_pair_and_distance_list::value_type,
            typename particle_id_pair_and_distance_list::value_type,
            bool>
    {
        typedef typename particle_id_pair_and_distance_list::value_type
                first_argument_type;

        bool operator()(first_argument_type const& lhs,
                        first_argument_type const& rhs) const
        {
            return lhs.second < rhs.second;
        }
    };

    template<typename Tset_>
    struct overlap_checker
    {
        overlap_checker(Tset_ const& ignore = Tset_()): ignore_(ignore), result_(0) {}

        template<typename Titer_>
        void operator()(Titer_ const& i, length_type const& dist)
        {
            if (!collection_contains(ignore_, (*i).first))
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
                // std::sort(result_->pbegin(), result_->pend(), compare_);
                std::sort(result_->begin(), result_->end(), compare_);
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
    typedef typename traits_type::molecule_info_type molecule_info_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator particle_id_generator;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::particle_id_pair_generator
        particle_id_pair_generator;
    typedef typename traits_type::particle_id_pair_and_distance
        particle_id_pair_and_distance;
    typedef typename traits_type::particle_id_pair_and_distance_list
        particle_id_pair_and_distance_list;

    typedef typename base_type::transaction_type transaction_type;
    typedef typename base_type::time_type time_type;

    typedef MatrixSpace<particle_type, particle_id_type, ecell4::utils::get_mapper_mf> particle_matrix_type;
    typedef sized_iterator_range<typename particle_matrix_type::const_iterator> particle_id_pair_range;
    typedef typename particle_matrix_type::matrix_sizes_type matrix_sizes_type;

public:
    ParticleContainerBase(const position_type& edge_lengths, const matrix_sizes_type& sizes)
        : pmat_(new particle_matrix_type(edge_lengths, sizes)), t_(0.0) {}

    virtual ecell4::Integer num_particles() const
    {
        return (*pmat_).size();
    }

    // virtual size_type num_particles() const
    // {
    //     return (*pmat_).size();
    // }

    virtual const position_type& edge_lengths() const
    {
        return (*pmat_).edge_lengths();
    }

    virtual void reset(const position_type& lengths)
    {
        const matrix_sizes_type sizes((*pmat_).matrix_sizes());
        this->reset(lengths, sizes);
    }

    virtual void reset(const position_type& lengths, const matrix_sizes_type& sizes)
    {
        (*pmat_).clear();
        boost::scoped_ptr<particle_matrix_type>
            newpmat(new particle_matrix_type(lengths, sizes));
        pmat_.swap(newpmat);

        ; // newpmat will be released here
    }

    position_type cell_sizes() const
    {
        return (*pmat_).cell_sizes();
    }

    matrix_sizes_type matrix_sizes() const
    {
        return (*pmat_).matrix_sizes();
    }

    template<typename T_>
    length_type distance(T_ const& lhs, position_type const& rhs) const
    {
        return traits_type::distance(lhs, rhs, edge_lengths());
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return traits_type::distance(lhs, rhs, edge_lengths());
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return traits_type::apply_boundary(v, edge_lengths());
    }

    // virtual length_type apply_boundary(length_type const& v) const
    // {
    //     return traits_type::apply_boundary(v, edge_lengths());
    // }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return traits_type::cyclic_transpose(p0, p1, edge_lengths());
    }

    // virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    // {
    //     return traits_type::cyclic_transpose(p0, p1, world_size());
    // }

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
            edge_lengths());
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
        traits_type::take_neighbor(*pmat_, oc, s);
        return oc.result();
    }

    template<typename Tsph_>
    particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s,
        typename boost::disable_if<boost::is_same<Tsph_, particle_id_pair> >::type* =0) const
    {
        typename utils::template overlap_checker<boost::array<particle_id_type, 0> > oc;
        traits_type::take_neighbor(*pmat_, oc, s);
        return oc.result();
    }

    particle_id_pair get_particle(particle_id_type const& id, bool& found) const
    {
        typename particle_matrix_type::const_iterator i((*pmat_).find(id));
        if ((*pmat_).end() == i) {
            found = false;
            return particle_id_pair();
        }
        found = true;
        return *i;
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        typename particle_matrix_type::const_iterator i((*pmat_).find(id));
        if ((*pmat_).end() == i) {
            throw not_found(std::string("No such particle: id=")
                    + boost::lexical_cast<std::string>(id));
        }
        return *i;
    }

    virtual bool has_particle(particle_id_type const& id) const
    {
        return (*pmat_).end() != (*pmat_).find(id);
    }

    virtual transaction_type* create_transaction();

    virtual particle_id_pair_generator* get_particles() const
    {
        return make_range_generator<particle_id_pair>(*pmat_);
    }

    particle_id_pair_range get_particles_range() const
    {
        return particle_id_pair_range((*pmat_).begin(), (*pmat_).end(), (*pmat_).size());
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        return (*pmat_).update(pi_pair).second;
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        return (*pmat_).erase(id);
    }

    /** ecell4::Space
     */
    // virtual const time_type& t() const
    virtual const time_type t() const
    {
        return t_;
    }

    virtual void set_t(const time_type& t)
    {
        t_ = t;
    }

    // for polygon
    virtual void add_surface(const boost::array<position_type, 3>& vertices)
    {
        polygon_.emplace(vertices);
    }

//     virtual position_type
//     apply_reflection(const position_type& pos, const position_type& disp)
//     {
//         return polygon_.apply_reflection(pos, disp,
//                 (polygon_.get_faces_within_radius(pos, length(disp))).first,
//                 this->edge_lengths());
//     }

    virtual position_type
    apply_structure(const position_type& pos, const position_type& disp)
    {
        return this->apply_structure_rec(pos, disp, Polygon<position_type>::make_nonsence_id());
    }

  protected:

    position_type
    apply_structure_rec(const position_type& pos, const position_type& disp,
            const typename Polygon<position_type>::face_id_type ignore)
    {
        typedef typename Polygon<position_type>::face_id_type face_id_t;

        const ecell4::AABBSurface unitcell(
                position_type(0., 0., 0.), this->edge_lengths());
        const std::pair<bool, length_type> test_unitcell =
                unitcell.intersect_ray(pos, disp);
        const length_type dist_to_unit_cell =
                length(disp) * test_unitcell.second;

        const std::pair<bool, std::pair<length_type, face_id_t> > test_polygon = 
                this->polygon_.intersect_ray(pos, disp, ignore);

        if(not test_unitcell.first && not test_polygon.first)
            return pos + disp;

        if(test_polygon.first && test_polygon.second.first < dist_to_unit_cell)
        {
            const std::pair<std::pair<position_type, position_type>, face_id_t>
                    reflected = this->polygon_.apply_reflection(
                            pos, disp, test_polygon.second.second);
            return this->apply_structure_rec(reflected.first.first,
                    reflected.first.second - reflected.first.first, reflected.second);
        }
        else if(test_unitcell.first)
        {
            if(test_unitcell.second <= 0.0 || 1.0 < test_unitcell.second)
            {
                std::cerr << "aabb.is_inside(begin) = " << unitcell._is_inside(pos) << std::endl;
                std::cerr << "begin = " << pos << std::endl;
                std::cerr << "edge_length = " << this->edge_lengths() << std::endl;
                std::cerr << "test_unitcell.first = " << test_unitcell.first << std::endl;
                std::cerr << "test_unitcell.second = " <<  test_unitcell.second  << std::endl;
                std::cerr << "test_polygon.first = " << test_polygon.first << std::endl;
                std::cerr << "test_polygon.second.first = "  << test_polygon.second.first << std::endl;
                std::cerr << "test_polygon.second.second = " << test_polygon.second.second << std::endl;
                assert(0);
            }
            const std::pair<position_type, position_type> next_segment =
                apply_periodic_only_once(pos, disp, test_unitcell.second, unitcell);
            return this->apply_structure_rec(
                    next_segment.first, next_segment.second - next_segment.first,
                    Polygon<position_type>::make_nonsence_id());
        }
    }

    std::pair<position_type, position_type>
    apply_periodic_only_once(const position_type& pos, const position_type& disp,
                             const length_type tmin, const ecell4::AABBSurface& aabb)
    {
        //XXX: this function assumes the conditions described below is satisfied.
        // - aabb.lower = (0, 0, 0)
        // - periodic boundary is applied
        assert(0. < tmin && tmin <= 1.0);
        position_type next_begin = pos + disp * tmin;
        position_type next_end   = pos + disp;
        position_type pullback;
             if(std::abs(next_begin[0] - aabb.upper()[0]) < 1e-12)
        {
            next_begin[0] = aabb.lower()[0];
            next_end[0] -= (aabb.upper()[0] - aabb.lower()[0]);
        }
        else if(std::abs(next_begin[0] - aabb.lower()[0]) < 1e-12)
        {
            next_begin[0] = aabb.upper()[0];
            next_end[0] += (aabb.upper()[0] - aabb.lower()[0]);
        }
        else if(std::abs(next_begin[1] - aabb.upper()[1]) < 1e-12)
        {
            next_begin[1] = aabb.lower()[1];
            next_end[1] -= (aabb.upper()[1] - aabb.lower()[1]);
        }
        else if(std::abs(next_begin[1] - aabb.lower()[1]) < 1e-12)
        {
            next_begin[1] = aabb.upper()[1];
            next_end[1] += (aabb.upper()[1] - aabb.lower()[1]);
        }
        else if(std::abs(next_begin[2] - aabb.upper()[2]) < 1e-12)
        {
            next_begin[2] = aabb.lower()[2];
            next_end[2] -= (aabb.upper()[2] - aabb.lower()[2]);
        }
        else if(std::abs(next_begin[2] - aabb.lower()[2]) < 1e-12)
        {
            next_begin[2] = aabb.upper()[2];
            next_end[2] += (aabb.upper()[2] - aabb.lower()[2]);
        }
        else
        {
            throw std::logic_error("never reach here");
        }
        assert(aabb._is_inside(next_begin));

        return std::make_pair(next_begin, next_end);
    }

protected:

    std::pair<bool, typename particle_matrix_type::iterator>
    __has_particle(const particle_id_type& pid)
    {
        typename particle_matrix_type::iterator i((*pmat_).find(pid));
        return std::make_pair(i != (*pmat_).end(), i);
    }

    typename particle_matrix_type::iterator __update_particle(
        const typename particle_matrix_type::iterator& position,
        const particle_id_pair& pi_pair)
    {
        return (*pmat_).update(position, pi_pair);
    }

    virtual void clear()
    {
        (*pmat_).clear();
        t_ = 0.0;
    }

protected:
    boost::scoped_ptr<particle_matrix_type> pmat_;
    Polygon<position_type> polygon_;

    time_type t_;
};

template<typename Tderived_, typename Ttraits_>
inline typename ParticleContainerBase<Tderived_, Ttraits_>::transaction_type*
ParticleContainerBase<Tderived_, Ttraits_>::create_transaction()
{
    return new TransactionImpl<ParticleContainerBase>(*this);
}

#endif /* PARTICLE_CONTAINER_BASE_HPP */
