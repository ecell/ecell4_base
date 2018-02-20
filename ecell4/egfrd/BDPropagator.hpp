#ifndef BD_PROPAGATOR_HPP
#define BD_PROPAGATOR_HPP

#include <algorithm>
#include <limits>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/scoped_ptr.hpp>
#include "Defs.hpp"
#include "generator.hpp"
#include "exceptions.hpp"
#include "utils.hpp"
#include "utils/random.hpp"
#include "utils/get_default_impl.hpp"
#include "Logger.hpp"

#include <ecell4/core/get_mapper_mf.hpp>
#include "PotentialField.hpp"
#include <greens_functions/freeFunctions.hpp>
// using namespace greens_functions;

template<typename Ttraits_>
class BDPropagator
{
public:
    typedef Ttraits_ traits_type;
    typedef typename Ttraits_::world_type::particle_container_type particle_container_type;
    typedef typename particle_container_type::species_id_type species_id_type;
    typedef typename particle_container_type::position_type position_type;
    typedef typename particle_container_type::particle_shape_type particle_shape_type;
    typedef typename particle_container_type::molecule_info_type molecule_info_type;
    typedef typename particle_container_type::length_type length_type;
    typedef typename particle_container_type::particle_id_type particle_id_type;
    typedef typename particle_container_type::particle_type particle_type;
    typedef typename particle_container_type::particle_id_pair particle_id_pair;
    typedef std::vector<particle_id_type> particle_id_vector_type;
    typedef typename particle_container_type::particle_id_pair_generator particle_id_pair_generator;
    typedef typename particle_container_type::particle_id_pair_and_distance particle_id_pair_and_distance;
    typedef typename particle_container_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;
    typedef typename particle_container_type::structure_type structure_type;
    typedef typename traits_type::world_type::traits_type::rng_type rng_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename network_rules_type::reaction_rules reaction_rules;
    typedef typename network_rules_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef typename traits_type::volume_clearer_type volume_clearer_type;

    typedef typename ecell4::utils::get_mapper_mf<particle_id_type, position_type>::type particle_id_position_map_type;

    typedef ecell4::PotentialField<particle_container_type> potential_field_type;
    typedef typename ecell4::utils::get_mapper_mf<species_id_type, boost::shared_ptr<potential_field_type> >::type potential_field_map_type;

public:
    template<typename Trange_>
    BDPropagator(
        particle_container_type& tx, network_rules_type const& rules,
        rng_type& rng, time_type dt, int max_retry_count,
        reaction_recorder_type* rrec, volume_clearer_type* vc,
        Trange_ const& particles,
        potential_field_map_type const& potentials = potential_field_map_type())
        : tx_(tx), rules_(rules), rng_(rng), dt_(dt),
          max_retry_count_(max_retry_count), rrec_(rrec), vc_(vc),
          queue_(), rejected_move_count_(0),
          potentials_(potentials)
    {
        call_with_size_if_randomly_accessible(
            boost::bind(&particle_id_vector_type::reserve, &queue_, _1),
            particles);
        for (typename boost::range_const_iterator<Trange_>::type
                i(boost::begin(particles)),
                e(boost::end(particles)); i != e; ++i)
        {
            queue_.push_back(*i);
        }
        shuffle(rng, queue_);
    }

    bool operator()()
    {
        if (queue_.empty())
            return false;

        particle_id_type pid(queue_.back());
        queue_.pop_back();
        particle_id_pair pp(tx_.get_particle(pid));

        LOG_DEBUG(("propagating particle %s", boost::lexical_cast<std::string>(pp.first).c_str()));

        try
        {
            if (attempt_reaction(pp))
                return true;
        }
        catch (propagation_error const& reason)
        {
            log_.info("first-order reaction rejected (reason: %s)", reason.what());
            ++rejected_move_count_;
            return true;
        }

        const species_id_type& species_id(pp.second.species());
        const molecule_info_type species(tx_.get_molecule_info(species_id));
        if (species.D == 0.)
            return true;

        position_type const displacement = drawR_free(species);
        position_type const new_pos    = tx_.apply_structure(pp.second.position(), displacement);
//         position_type const new_pos      = tx_.apply_boundary(reflected);

        // typename potential_field_map_type::const_iterator it = potentials_.find(species_id);
        // if (it != potentials_.end())
        // {
        //     if (!(*it).second->try_move(rng_, pp, new_pos, tx_))
        //     {
        //         return true;
        //     }
        // }

        particle_id_pair particle_to_update(
                pp.first, particle_type(species_id,
                    new_pos, species.radius,
                    species.D));
        particle_id_pair_and_distance_list overlapped(
            tx_.check_overlap(shape(particle_to_update.second),
                              particle_to_update.first));
        switch (overlapped.size())
        {
        case 0:
            break;

        case 1:
            {
                particle_id_pair_and_distance const& closest(overlapped.at(0));
                try
                {
                    if (!attempt_reaction(pp, closest.first))
                    {
                        LOG_DEBUG(("collision with a nonreactive particle %s. move rejected", boost::lexical_cast<std::string>(closest.first.first).c_str()));
                        ++rejected_move_count_;
                    }
                }
                catch (propagation_error const& reason)
                {
                    log_.info("second-order reaction rejected (reason: %s)", reason.what());
                    ++rejected_move_count_;
                }
            }
            /* reject the move even if the reaction has not occurred */
            return true;

        default:
            log_.info("collision involving two or more particles; move rejected");
            ++rejected_move_count_;
            return true;
        }
        if (vc_)
        {
            if (!(*vc_)(shape(particle_to_update.second), 
                        particle_to_update.first))
            {
                log_.info("propagation move rejected.");
                return true;
            }
        }
        tx_.update_particle(particle_to_update.first, particle_to_update.second);
        return true;
    }

    std::size_t get_rejected_move_count() const
    {
        return rejected_move_count_;
    }

private:
    position_type drawR_free(molecule_info_type const& species)
    {
        return tx_.get_structure(species.structure_id)->bd_displacement(std::sqrt(2.0 * species.D * dt_), rng_);
    }

    bool attempt_reaction(particle_id_pair const& pp)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.species()));
        if (::size(rules) == 0)
        {
            return false;
        }

        const Real rnd(rng_.random() / dt_);
        Real prob = 0.;

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);
            prob += r.k();
            if (prob > rnd)
            {
                typename reaction_rule_type::species_id_range products(
                        r.get_products());
                switch (::size(products))
                {
                case 0:
                    remove_particle(pp.first);
                    break;

                case 1:
                    {
                        const molecule_info_type s0(tx_.get_molecule_info(products[0]));
                        const particle_id_pair new_p(
                            pp.first, particle_type(products[0],
                                pp.second.position(), s0.radius, s0.D));
                        if (!tx_.no_overlap(shape(new_p.second), new_p.first))
                        {
                            throw propagation_error("no space");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(shape(new_p.second), pp.first))
                            {
                                throw propagation_error("no space");
                            }
                        }

                        tx_.update_particle(new_p.first, new_p.second);

                        if (rrec_)
                        {
                            // (*rrec_)(
                            //     reaction_record_type(
                            //         r.id(), array_gen(new_p.first), pp.first));
                            (*rrec_)(
                                reaction_record_type(r.id(), array_gen(new_p), pp));
                        }
                    }
                    break;

                case 2:
                    {
                        const species_id_type& product_id0(products[0]),
                            product_id1(products[1]);
                        const molecule_info_type s0(tx_.get_molecule_info(product_id0)),
                                s1(tx_.get_molecule_info(product_id1));
                        const Real D01(s0.D + s1.D);
                        const length_type r01(s0.radius + s1.radius);
                        int i = max_retry_count_;
                        position_type np0, np1;

                        for (;;)
                        {
                            if (--i < 0)
                            {
                                throw propagation_error("no space");
                            }

                            const Real rnd(rng_.random());
                            length_type pair_distance(
                                greens_functions::drawR_gbd_3D(rnd, r01, dt_, D01));
                            const position_type m(random_unit_vector() * pair_distance);
                            np0 = tx_.apply_boundary(pp.second.position()
                                    + m * (s0.D / D01));
                            np1 = tx_.apply_boundary(pp.second.position()
                                    - m * (s1.D / D01));

                            const particle_shape_type sphere1(np0, s0.radius);
                            const particle_shape_type sphere2(np1, s1.radius);
                            if (tx_.no_overlap(sphere1, pp.first)
                                && tx_.no_overlap(sphere2, pp.first))
                            {
                                break;
                            }
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(particle_shape_type(np0, s0.radius), pp.first) || !(*vc_)(particle_shape_type(np1, s1.radius), pp.first))
                            {
                                throw propagation_error("no space");
                            }
                        }

                        tx_.remove_particle(pp.first);
                        const particle_id_pair
                            npp0(tx_.new_particle(product_id0, np0).first),
                            npp1(tx_.new_particle(product_id1, np1).first);

                        if (rrec_)
                        {
                            // (*rrec_)(
                            //     reaction_record_type(
                            //         r.id(),
                            //         array_gen(npp0.first, npp1.first),
                            //         pp.first));
                            (*rrec_)(reaction_record_type(r.id(), array_gen(npp0, npp1), pp));
                        }
                    }
                    break;
                default:
                    throw not_implemented("monomolecular reactions that produce more than two products are not supported");
                }
                return true;
            }
        }
        return false;
    }

    bool attempt_reaction(particle_id_pair const& pp0, particle_id_pair const& pp1)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp0.second.species(), pp1.second.species()));
        if (::size(rules) == 0)
        {
            return false;
        }

        const molecule_info_type s0(tx_.get_molecule_info(pp0.second.species())),
                s1(tx_.get_molecule_info(pp1.second.species()));
        const length_type r01(s0.radius + s1.radius);

        const Real rnd(rng_.random());
        Real prob = 0;

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);
            const Real p(r.k() * dt_ / ((greens_functions::I_bd_3D(r01, dt_, s0.D) + greens_functions::I_bd_3D(r01, dt_, s1.D)) * 4.0 * M_PI));
            BOOST_ASSERT(p >= 0.);
            prob += p;
            if (prob >= 1.)
            {
                throw propagation_error(
                    "invalid acceptance ratio ("
                    + boost::lexical_cast<std::string>(p)
                    + ") for reaction rate "
                    + boost::lexical_cast<std::string>(r.k())
                    + ".");
            }
            if (prob > rnd)
            {
                LOG_DEBUG(("fire reaction"));
                const typename reaction_rule_type::species_id_range products(
                    r.get_products());

                switch (::size(products))
                {
                case 1:
                    {
                        const species_id_type product(products[0]);
                        const molecule_info_type sp(tx_.get_molecule_info(product));

                        const position_type new_pos(
                            tx_.apply_boundary(
                                divide(
                                    add(multiply(pp0.second.position(), s1.D),
                                        multiply(tx_.periodic_transpose(
                                            pp1.second.position(),
                                            pp0.second.position()), s0.D)),
                                    (s0.D + s1.D))));
                        if (!tx_.no_overlap(
                            particle_shape_type(new_pos, sp.radius),
                            pp0.first, pp1.first))
                        {
                            throw propagation_error("no space");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(
                                    particle_shape_type(new_pos, sp.radius), 
                                    pp0.first, pp1.first))
                            {
                                throw propagation_error("no space");
                            }
                        }

                        remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        particle_id_pair npp(tx_.new_particle(product, new_pos).first);
                        if (rrec_)
                        {
                            // (*rrec_)(
                            //     reaction_record_type(
                            //         r.id(), array_gen(npp.first), pp0.first, pp1.first));
                            (*rrec_)(reaction_record_type(r.id(), array_gen(npp), pp0, pp1));
                        }
                        break;
                    }
                case 0:
                    remove_particle(pp0.first);
                    remove_particle(pp1.first);
                    break;
                default:
                    throw not_implemented("bimolecular reactions that produce more than one product are not supported");
                }

                return true;
            }
        }
        return false;
    }

    void remove_particle(particle_id_type const& pid)
    {
        LOG_DEBUG(("remove particle %s", boost::lexical_cast<std::string>(pid).c_str()));
        tx_.remove_particle(pid);
        typename particle_id_vector_type::iterator i(
            std::find(queue_.begin(), queue_.end(), pid));
        if (queue_.end() != i)
            queue_.erase(i);
    }

private:
    position_type random_unit_vector()
    {
        position_type v(rng_.random() - 0.5, rng_.random() - 0.5, rng_.random() - 0.5);
        return v / length(v);
    }

private:
    particle_container_type& tx_;
    network_rules_type const& rules_;
    rng_type& rng_;
    Real const dt_;
    int const max_retry_count_;
    reaction_recorder_type* const rrec_;
    volume_clearer_type* const vc_;
    particle_id_vector_type queue_;
    int rejected_move_count_;
    potential_field_map_type const& potentials_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& BDPropagator<Ttraits_>::log_(Logger::get_logger("ecell.BDPropagator"));

#endif /* BD_PROPAGATOR_HPP */

