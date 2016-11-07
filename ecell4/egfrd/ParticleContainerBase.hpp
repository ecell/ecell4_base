#ifndef PARTICLE_CONTAINER_BASE_HPP
#define PARTICLE_CONTAINER_BASE_HPP

// #include <ecell4/core/get_mapper_mf.hpp>
// #include <ecell4/core/ParticleSpaceCellListImpl.hpp>
// 
// #include "utils/range.hpp"
// #include "utils/unassignable_adapter.hpp"
// #include "MatrixSpace.hpp"
// #include "abstract_set.hpp"
// #include "generator.hpp"
// #include "exceptions.hpp"
// #include "ParticleContainer.hpp"
// #include "Transaction.hpp"
// 
// template<typename Ttraits_>
// struct ParticleContainerUtils
// {
//     typedef Ttraits_ traits_type;
//     typedef typename traits_type::length_type length_type;
//     typedef typename traits_type::particle_type particle_type;
//     typedef typename traits_type::particle_id_type particle_id_type;
//     typedef typename traits_type::particle_id_pair particle_id_pair;
//     typedef typename traits_type::particle_id_pair_and_distance
//         particle_id_pair_and_distance;
//     typedef typename traits_type::particle_id_pair_and_distance_list
//         particle_id_pair_and_distance_list;
// 
//     // struct distance_comparator:
//     //         public std::binary_function<
//     //             typename particle_id_pair_and_distance_list::placeholder,
//     //             typename particle_id_pair_and_distance_list::placeholder,
//     //             bool>
//     // {
//     //     typedef typename particle_id_pair_and_distance_list::placeholder
//     //             first_argument_type;
//     //     typedef typename particle_id_pair_and_distance_list::const_caster const_caster;
//     //     bool operator()(first_argument_type const& lhs,
//     //                     first_argument_type const& rhs) const
//     //     {
//     //         return c_(lhs).second < c_(rhs).second;
//     //     }
//     //     const_caster c_;
//     // };
// 
//     struct distance_comparator
//         : public std::binary_function<
//             typename particle_id_pair_and_distance_list::value_type,
//             typename particle_id_pair_and_distance_list::value_type,
//             bool>
//     {
//         typedef typename particle_id_pair_and_distance_list::value_type
//                 first_argument_type;
// 
//         bool operator()(first_argument_type const& lhs,
//                         first_argument_type const& rhs) const
//         {
//             return lhs.second < rhs.second;
//         }
//     };
// 
//     template<typename Tset_>
//     struct overlap_checker
//     {
//         overlap_checker(Tset_ const& ignore = Tset_()): ignore_(ignore), result_(0) {}
// 
//         template<typename Titer_>
//         void operator()(Titer_ const& i, length_type const& dist)
//         {
//             if (!collection_contains(ignore_, (*i).first))
//             {
//                 if (!result_)
//                 {
//                     result_ = new particle_id_pair_and_distance_list();
//                 }
//                 result_->push_back(std::make_pair(*i, dist));
//             }
//         }
// 
//         particle_id_pair_and_distance_list* result() const
//         {
//             if (result_)
//             {
//                 // std::sort(result_->pbegin(), result_->pend(), compare_);
//                 std::sort(result_->begin(), result_->end(), compare_);
//             }
//             return result_;
//         }
// 
//     private:
//         Tset_ const& ignore_;
//         particle_id_pair_and_distance_list* result_;
//         distance_comparator compare_;
//     };
// };
// 
// template<typename Tderived_, typename Ttraits_ = typename Tderived_::traits_type>
// class ParticleContainerBase
//     : public ParticleContainer<Ttraits_>
// {
// public:
// 
//     typedef ParticleContainer<Ttraits_> base_type;
// 
//     typedef ParticleContainerUtils<Ttraits_> utils;
//     typedef Ttraits_ traits_type;
//     typedef typename traits_type::length_type length_type;
//     typedef typename traits_type::molecule_info_type molecule_info_type;
//     typedef typename traits_type::position_type position_type;
//     typedef typename traits_type::particle_type particle_type;
//     typedef typename traits_type::particle_id_type particle_id_type;
//     typedef typename traits_type::particle_id_generator particle_id_generator;
//     typedef typename traits_type::species_id_type species_id_type;
//     typedef typename traits_type::particle_shape_type particle_shape_type;
//     typedef typename traits_type::size_type size_type;
//     typedef typename traits_type::structure_id_type structure_id_type;
//     typedef typename traits_type::structure_type structure_type;
//     typedef typename traits_type::particle_id_pair particle_id_pair;
//     typedef typename traits_type::particle_id_pair_generator
//         particle_id_pair_generator;
//     typedef typename traits_type::particle_id_pair_and_distance
//         particle_id_pair_and_distance;
//     typedef typename traits_type::particle_id_pair_and_distance_list
//         particle_id_pair_and_distance_list;
// 
//     typedef typename base_type::transaction_type transaction_type;
//     typedef typename base_type::time_type time_type;
// 
//     typedef MatrixSpace<particle_type, particle_id_type, ecell4::utils::get_mapper_mf> particle_matrix_type;
//     typedef sized_iterator_range<typename particle_matrix_type::const_iterator> particle_id_pair_range;
//     typedef typename particle_matrix_type::matrix_sizes_type matrix_sizes_type;
// 
//     typedef ecell4::ParticleSpaceCellListImpl particle_space_type;
// 
// public:
// 
//     ParticleContainerBase(
//         const position_type& edge_lengths, const matrix_sizes_type& matrix_sizes)
//         : ps_(new particle_space_type(edge_lengths, matrix_sizes))
//     {
//         ;
//     }
// 
//     /**
//      * redirects
//      */
// 
//     virtual ecell4::Integer num_particles() const
//     {
//         return (*ps_).num_particles();
//     }
// 
//     virtual const position_type& edge_lengths() const
//     {
//         return (*ps_).edge_lengths();
//     }
// 
//     virtual void reset(const position_type& lengths)
//     {
//         (*ps_).reset(lengths);
//     }
// 
//     virtual void reset(const position_type& lengths, const matrix_sizes_type& sizes)
//     {
//         boost::scoped_ptr<particle_space_type>
//             newps(new particle_space_type(lengths, sizes));
//         ps_.swap(newps);
// 
//         ; // newpmat will be released here
//     }
// 
//     position_type cell_sizes() const
//     {
//         return (*ps_).cell_sizes();
//     }
// 
//     matrix_sizes_type matrix_sizes() const
//     {
//         return (*ps_).matrix_sizes();
//     }
// 
//     virtual bool has_particle(particle_id_type const& id) const
//     {
//         return (*ps_).has_particle(id);
//     }
// 
//     virtual particle_id_pair get_particle(particle_id_type const& id) const
//     {
//         return (*ps_).get_particle(id);
//     }
// 
//     virtual length_type distance(
//         position_type const& lhs, position_type const& rhs) const
//     {
//         return (*ps_).distance(lhs, rhs);
//     }
// 
//     virtual position_type apply_boundary(position_type const& v) const
//     {
//         return (*ps_).apply_boundary(v);
//     }
// 
//     virtual const time_type t() const
//     {
//         return (*ps_).t();
//     }
// 
//     virtual void set_t(const time_type& t)
//     {
//         (*ps_).set_t(t);
//     }
// 
//     /**
//      * wrappers
//      */
// 
//     virtual position_type periodic_transpose(
//         position_type const& p0, position_type const& p1) const
//     {
//         return (*ps_).periodic_transpose(p0, p1);
//     }
// 
//     template<typename T1_>
//     T1_ calculate_pair_CoM(
//         T1_ const& p1, T1_ const& p2,
//         typename element_type_of<T1_>::type const& D1,
//         typename element_type_of<T1_>::type const& D2)
//     {
//         typedef typename element_type_of<T1_>::type element_type;
// 
//         const T1_ p2_trans(periodic_transpose(p2, p1));
//         const element_type D12(add(D1, D2));
//         const element_type s(divide(D1, D12)), t(divide(D2, D12));
//         const T1_ com(add(multiply(p1, t), multiply(p2_trans, s)));
//         return apply_boundary(com);
//     }
// 
//     virtual bool remove_particle(particle_id_type const& id)
//     {
//         if (!has_particle(id))
//         {
//             return false;
//         }
//         (*ps_).remove_particle(id);
//         return true;
//     }
// 
//     virtual bool update_particle(particle_id_pair const& pid_particle_pair)
//     {
//         return (*ps_).update_particle(pid_particle_pair.first, pid_particle_pair.second);
//     }
// 
//     particle_id_pair get_particle(particle_id_type const& id, bool& found) const
//     {
//         found = (*ps_).has_particle(id);
//         if (!found)
//         {
//             return particle_id_pair();
//         }
//         return get_particle(id);
//     }
// 
//     virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s) const
//     {
//         return (*ps_).list_particles_within_radius(s.position(), s.radius());
//     }
// 
//     virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
//     {
//         return (*ps_).list_particles_within_radius(s.position(), s.radius(), ignore);
//     }
// 
//     virtual particle_id_pair_and_distance_list check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
//     {
//         return (*ps_).list_particles_within_radius(s.position(), s.radius(), ignore1, ignore2);
//     }
// 
//     // template<typename Tsph_, typename Tset_>
//     // particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s, Tset_ const& ignore,
//     //     typename boost::disable_if<boost::is_same<Tsph_, particle_id_pair> >::type* = 0) const
//     // {
//     //     typename utils::template overlap_checker<Tset_> oc(ignore);
//     //     traits_type::take_neighbor(*pmat_, oc, s);
//     //     return oc.result();
//     // }
// 
//     // template<typename Tsph_>
//     // particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s,
//     //     typename boost::disable_if<boost::is_same<Tsph_, particle_id_pair> >::type* = 0) const
//     // {
//     //     typename utils::template overlap_checker<boost::array<particle_id_type, 0> > oc;
//     //     traits_type::take_neighbor(*pmat_, oc, s);
//     //     return oc.result();
//     // }
// 
//     particle_id_pair_range get_particles_range() const
//     {
//         const particle_space_type::particle_container_type& particles((*ps_).particles());
//         return particle_id_pair_range(particles.begin(), particles.end(), particles.size());
//     }
// 
//     /**
//      *
//      */
// 
//     virtual transaction_type* create_transaction();
// 
//     template<typename T_>
//     length_type distance(T_ const& lhs, position_type const& rhs) const
//     {
//         // return (*ps_).distance(lhs, rhs);
//         return traits_type::distance(lhs, rhs, edge_lengths());
//     }
// 
//     void clear()
//     {
//         (*ps_).reset((*ps_).edge_lengths());
//     }
// 
// protected:
// 
//     boost::scoped_ptr<particle_space_type> ps_;
// };
// 
// template<typename Tderived_, typename Ttraits_>
// inline typename ParticleContainerBase<Tderived_, Ttraits_>::transaction_type*
// ParticleContainerBase<Tderived_, Ttraits_>::create_transaction()
// {
//     return new TransactionImpl<ParticleContainerBase>(*this);
// }

#endif /* PARTICLE_CONTAINER_BASE_HPP */
