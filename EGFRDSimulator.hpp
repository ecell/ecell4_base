#ifndef EGFRDSIMULATOR_HPP
#define EGFRDSIMULATOR_HPP

#include <boost/bind.hpp>
#include <boost/array.hpp>
#include <boost/format.hpp>
#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/sequence/intrinsic/at_key.hpp>
#include <boost/fusion/sequence/intrinsic/value_at_key.hpp>
#include <boost/fusion/include/at_key.hpp>
#include "utils/array_helper.hpp"
#include "utils/get_mapper_mf.hpp"
#include "utils/fun_composition.hpp"
#include "utils/fun_wrappers.hpp"
#include "utils/pointer_as_ref.hpp"
#include "utils/pair.hpp"
#include "utils/math.hpp"
#include "ShellID.hpp"
#include "DomainID.hpp"
#include "Shell.hpp"
#include "EventScheduler.hpp"
#include "PairGreensFunction.hpp"
#include "ParticleSimulator.hpp"
#include "MatrixSpace.hpp"
#include "AnalyticalSingle.hpp"
#include "AnalyticalPair.hpp"
#include "Multi.hpp"
#include "GreensFunction3DRadAbs.hpp"
#include "GreensFunction3DRadInf.hpp"
#include "GreensFunction3DAbsSym.hpp"
#include "GreensFunction3DAbs.hpp"
#include "GreensFunction3D.hpp"

template<typename Tworld_>
struct EGFRDSimulatorTraitsBase: public ParticleSimulatorTraitsBase<Tworld_>
{
    typedef ParticleSimulatorTraitsBase<Tworld_> base_type;
    typedef Tworld_ world_type;
    typedef ShellID shell_id_type;
    typedef DomainID domain_id_type;
    typedef typename ParticleSimulatorTraitsBase<Tworld_>::sphere_type sphere_type;
    typedef typename ParticleSimulatorTraitsBase<Tworld_>::cylinder_type cylinder_type;
    typedef SerialIDGenerator<shell_id_type> shell_id_generator;
    typedef SerialIDGenerator<domain_id_type> domain_id_generator;
    typedef Domain<EGFRDSimulatorTraitsBase> domain_type;
    typedef std::pair<const domain_id_type, boost::shared_ptr<domain_type> > domain_id_pair;
    typedef int event_id_type;
    typedef EventScheduler<typename base_type::time_type> event_scheduler_type;
    typedef typename event_scheduler_type::Event event_type;
    typedef typename event_scheduler_type::value_type event_id_pair_type;

    template<typename Tshape_>
    struct shell_generator
    {
        typedef Shell<Tshape_, domain_id_type> type;
    };

    static const Real SAFETY = 1. + 1e-5;
    static const Real SINGLE_SHELL_FACTOR = .1;
    static const Real DEFAULT_DT_FACTOR = 1e-5;
    static const Real CUTOFF_FACTOR = 5.6;
};

namespace detail {

template<typename T_>
struct get_greens_function {};

template<typename T_>
struct get_greens_function<Sphere<T_> >
{
    typedef GreensFunction3DAbsSym type;
};

template<typename T_>
struct get_greens_function<Cylinder<T_> >
{
    typedef GreensFunction3DAbsSym type;
};

template<typename T_>
struct get_pair_greens_function {};

template<typename T_>
struct get_pair_greens_function<Sphere<T_> >
{
    typedef GreensFunction3DRadAbs iv_type;
    typedef GreensFunction3DAbsSym com_type;
};

template<typename T_>
struct get_pair_greens_function<Cylinder<T_> >
{
    typedef GreensFunction3DRadAbs iv_type;
    typedef GreensFunction3DAbsSym com_type;
};

} // namespace detail

template<typename Ttraits_>
class EGFRDSimulator: public ParticleSimulator<Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef ParticleSimulator<Ttraits_> base_type;
    typedef typename traits_type::world_type world_type;
    typedef typename traits_type::domain_id_type domain_id_type;
    typedef typename traits_type::shell_id_type shell_id_type;
    typedef typename traits_type::sphere_type sphere_type;
    typedef typename traits_type::cylinder_type cylinder_type;
    typedef typename traits_type::template shell_generator<sphere_type>::type spherical_shell_type;
    typedef typename traits_type::template shell_generator<cylinder_type>::type cylindrical_shell_type;
    typedef std::pair<const shell_id_type, spherical_shell_type> spherical_shell_id_pair;
    typedef std::pair<const shell_id_type, cylindrical_shell_type> cylindrical_shell_id_pair;
    typedef typename traits_type::shell_id_generator shell_id_generator;
    typedef typename traits_type::domain_id_generator domain_id_generator;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename world_type::traits_type::length_type length_type;
    typedef typename world_type::traits_type::position_type position_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename world_type::traits_type::particle_type particle_type;
    typedef typename world_type::traits_type::D_type D_type;
    typedef typename world_type::traits_type::species_type species_type;
    typedef typename world_type::traits_type::species_id_type species_id_type;
    typedef typename world_type::traits_type::structure_type structure_type;
    typedef typename world_type::particle_shape_type particle_shape_type;
    typedef typename world_type::traits_type::particle_id_type particle_id_type;
    typedef typename world_type::particle_id_pair particle_id_pair;
    typedef typename world_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;

    typedef typename traits_type::domain_type domain_type;
    typedef typename traits_type::domain_id_pair domain_id_pair;
    typedef typename traits_type::time_type time_type;

    typedef Single<traits_type> single_type;
    typedef Pair<traits_type> pair_type;
    typedef Multi<EGFRDSimulator> multi_type;
    typedef AnalyticalSingle<traits_type, spherical_shell_type> spherical_single_type;
    typedef AnalyticalSingle<traits_type, cylindrical_shell_type> cylindrical_single_type;
    typedef AnalyticalPair<traits_type, spherical_shell_type> spherical_pair_type;
    typedef AnalyticalPair<traits_type, cylindrical_shell_type> cylindrical_pair_type;

    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef typename traits_type::cylindrical_surface_type cylindrical_surface_type;
    typedef typename traits_type::planar_surface_type planar_surface_type;
    typedef typename traits_type::cuboidal_region_type cuboidal_region_type;
    typedef typename traits_type::event_scheduler_type event_scheduler_type;
    typedef typename traits_type::event_type event_type;
    typedef typename traits_type::event_id_type event_id_type;
    typedef typename traits_type::event_id_pair_type event_id_pair_type;

protected:
    typedef boost::fusion::map<
        boost::fusion::pair<spherical_shell_type, 
                            MatrixSpace<spherical_shell_type,
                                        shell_id_type, get_mapper_mf>&>,
        boost::fusion::pair<cylindrical_shell_type, MatrixSpace<cylindrical_shell_type,
                                        shell_id_type, get_mapper_mf>&> >
            shell_matrix_map_type;
    typedef typename boost::remove_reference<
        typename boost::fusion::result_of::value_at_key<
            shell_matrix_map_type,
            spherical_shell_type>::type>::type
                spherical_shell_matrix_type;
    typedef typename boost::remove_reference<
        typename boost::fusion::result_of::value_at_key<
            shell_matrix_map_type,
            cylindrical_shell_type>::type>::type
                cylindrical_shell_matrix_type;
    typedef typename get_mapper_mf<domain_id_type, boost::shared_ptr<domain_type> >::type domain_map;
    typedef typename network_rules_type::reaction_rules reaction_rules;
    typedef typename network_rules_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::rate_type rate_type;

    struct domain_event_base: public event_type
    {
        domain_event_base(time_type const& time): event_type(time) {}

        virtual domain_type& domain() const = 0;
    };

    template<typename Tdomain_, typename TeventKind_>
    class domain_event: public domain_event_base
    {
    public:
        typedef domain_event_base base_type;
        typedef Tdomain_ domain_type;
        typedef TeventKind_ event_kind_type;

    public:
        virtual domain_type& domain() const { return domain_; }

        virtual ~domain_event() {}

        event_kind_type kind() const { return kind_; }

        domain_event(time_type const& time,
                     domain_type& domain,
                     event_kind_type kind)
            : base_type(time), domain_(domain), kind_(kind) {}

    private:
        domain_type& domain_;
        event_kind_type kind_;
    };

    enum single_event_kind
    {
        SINGLE_EVENT_REACTION,
        SINGLE_EVENT_ESCAPE
    };

    enum pair_event_kind
    {
        PAIR_EVENT_SINGLE_REACTION_0,
        PAIR_EVENT_SINGLE_REACTION_1,
        PAIR_EVENT_COM_ESCAPE,
        PAIR_EVENT_IV
    };

    typedef domain_event<single_type, single_event_kind> single_event;
    typedef domain_event<pair_type, pair_event_kind> pair_event;

    class multi_event: public domain_event_base
    {
    public:
        typedef domain_event_base base_type;
        typedef multi_type domain_type;

    public:
        virtual domain_type& domain() const { return domain_; }

        virtual ~multi_event() {}

        multi_event(time_type const& time,
                     domain_type& domain)
            : base_type(time), domain_(domain) {}

    private:
        domain_type& domain_;
    };

    template<typename TdistFn_>
    struct intruder_collector
    {
        typedef TdistFn_ distance_function;

        intruder_collector(distance_function const& dist_fn,
                           domain_id_type const& ignore)
            : dist_fn(dist_fn), ignore(ignore),
              closest(domain_id_type(),
                      std::numeric_limits<length_type>::infinity()) {}

        template<typename Titer>
        void operator()(Titer const& i, position_type const& off)
        {
            domain_id_type const& did((*i).second.did());
            if (did == ignore)
                return;

            length_type const distance(dist_fn(shape(offset((*i).second, off))));
            if (distance < closest.second)
            {
                closest.first = did;
                closest.second = distance;
            }
            else
            {
                if (!intruders.container())
                {
                    intruders.container().set(new std::vector<domain_id_type>());
                }
                intruders.push_no_duplicate(did);
            }
        }

        distance_function const& dist_fn;
        domain_id_type ignore;
        std::pair<domain_id_type, length_type> closest;
        sorted_list<std::vector<domain_id_type>,
                    std::less<domain_id_type>,
                    pointer_as_ref<std::vector<domain_id_type> > > intruders;
    };

    struct no_filter
    {
        bool operator()(domain_id_type const&) const { return true; }
    };

    struct one_id_filter
    {
        bool operator()(domain_id_type const& did) const { return did != ignore; }
        one_id_filter(domain_id_type const& did): ignore(did) {}

        domain_id_type const ignore;
    };

    template<typename TfilterFn_>
    struct domain_collector 
    {
        typedef TfilterFn_ filter_function;

        domain_collector(world_type const& world,
                        particle_shape_type const& cmp,
                        filter_function const& filter)
            : world(world), cmp(cmp), filter(filter) {}

        template<typename Titer>
        void operator()(Titer const& i, position_type const& off)

        {
            domain_id_type const& did((*i).second.did());
            if (!filter(did))
                return;

            length_type const distance(world.distance(shape(offset((*i).second, off)), cmp.position()));
            if (distance < cmp.radius())
            {
                if (!neighbors.container())
                {
                    neighbors.container().set(new std::vector<domain_id_type>());
                }
                neighbors.push_no_duplicate(did);
            }
        }

        world_type const& world;
        particle_shape_type cmp;
        filter_function const& filter;

        std::pair<domain_id_type, length_type> closest;
        sorted_list<std::vector<domain_id_type>,
                    std::less<domain_id_type>,
                    pointer_as_ref<std::vector<domain_id_type> > > neighbors;
    };

    template<typename TdistFn_, typename TdidSet_>
    struct closest_object_finder
    {
        typedef TdistFn_ distance_function;
        typedef TdidSet_ domain_id_set;

        closest_object_finder(distance_function const& dist_fn,
                              domain_id_set const& ignore)
            : dist_fn(dist_fn), ignore(ignore),
              closest(domain_id_type(),
                      std::numeric_limits<length_type>::infinity()) {}

        template<typename Titer>
        void operator()(Titer const& i, position_type const& off)
        {
            domain_id_type const& did((*i).second.did());
            if (contains(ignore, did))
                return;

            length_type const distance(dist_fn(shape(offset((*i).second, off))));
            if (distance < closest.second)
            {
                closest.first = did;
                closest.second = distance;
            }
        }

        distance_function const& dist_fn;
        domain_id_set const& ignore;
        std::pair<domain_id_type, length_type> closest;
        sorted_list<std::vector<domain_id_type>,
                    std::less<domain_id_type>,
                    pointer_as_ref<std::vector<domain_id_type> > > intruders;
    };

    template<typename Tcol_>
    struct shell_collector_applier
    {
        typedef Tcol_ collector_type;

        shell_collector_applier(collector_type& col,
                                   position_type const& pos)
            : col_(col), pos_(pos) {}

        template<typename T>
        void operator()(T const& smat) const
        {
            world_type::traits_type::each_neighbor(smat.second, col_, pos_);
        }

    private:
        collector_type& col_;
        position_type pos_;
    };

    struct distance_calculator
    {
        distance_calculator(world_type const& world,
                            position_type const& pos)
            : world_(world), pos_(pos) {}

        template<typename Tshape_>
        length_type operator()(Tshape_ const& shape) const
        {
            return world_.distance(shape, pos_);
        }

    private:
        world_type const& world_;
        position_type pos_;
    };

    template<typename Tmap_>
    struct domain_shell_map_builder
    {
        typedef Tmap_ domain_shell_association;
        domain_shell_map_builder(world_type const& world,
                                 domain_shell_association& did_map)
            : world_(world), did_map_(did_map) {}

        template<typename T>
        void operator()(T const& smat) const
        {
            BOOST_ASSERT(world_.world_size() == smat.second.world_size());
            BOOST_FOREACH (typename boost::remove_reference<typename T::second_type>::type::value_type pair, smat.second)
            {
                did_map_[pair.second.did()].push_back(pair.first);
            }
        }

    private:
        world_type const& world_;
        domain_shell_association& did_map_;
    };

    struct shell_population_summer
    {
        shell_population_summer(std::size_t& population)
            : population_(population) {}

        template<typename T>
        void operator()(T const& smat) const
        {
            population_ += smat.second.size();
        }

    private:
        std::size_t& population_;
    };


    struct draw_on_com_escape
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt) const
        {
            return add(
                domain.shell().second.position(),
                normalize(
                    create_vector<position_type>(
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.)),
                    domain.a_R()));

        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv) const
        {
            boost::scoped_ptr<PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(draw_r(rng_, *gf, dt, domain.a_R(), domain.sigma()));
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    array_gen(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt) const
        {
            boost::shared_ptr<structure_type> const _structure(
                world_.get_structure(
                    world_.get_species(
                        domain.particles()[0].second.sid())
                    .structure_id()));
            
            cylindrical_surface_type const* const structure(
                dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            return add(
                domain.shell().second.position(),
                multiply(structure->shape().unit_z(), domain.a_R()));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv) const
        {
            BOOST_ASSERT(::size(domain.reactions()) == 1);
            length_type const r(
                draw_r(rng_, GreensFunction3DRadAbs(domain.D_tot(),
                    domain.reactions()[0].k(), domain.r0(),
                    domain.sigma(), domain.a_r()),
                   dt, domain.a_r(), domain.sigma()));
            BOOST_ASSERT(r > domain.sigma() && r <= domain.a_r());
            return multiply(normalize(old_iv), r);
        }

        draw_on_com_escape(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };

    struct draw_on_single_reaction
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt) const
        {
            return add(
                domain.shell().second.position(),
                draw_r(rng_,
                        GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                        dt, domain.a_R()));

        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv) const
        {
            boost::scoped_ptr<PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(draw_r(
                rng_, *gf, dt, domain.a_R(), domain.sigma()));
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    array_gen(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt) const
        {
            boost::shared_ptr<structure_type> const _structure(
                world_.get_species(
                    domain.particles()[0].second.sid())
                .structure_id());
            
            cylindrical_surface_type const* const structure(
                dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            BOOST_ASSERT(structure);

            return add(
                domain.shell().second.position(),
                multiply(structure->shape().unit_z(), domain.a_R()));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv) const
        {
            BOOST_ASSERT(::size(domain.reactions()) == 1);
            length_type const r(
                draw_r(rng_, GreensFunction3DRadAbs(domain.D_tot(),
                    domain.reactions()[0].k(), domain.r0(),
                    domain.sigma(), domain().a_r()),
                   dt, domain.a_r(), domain.sigma()));
            BOOST_ASSERT(r > domain.sigma() && r <= domain.a_r());
            return multiply(normalize(old_iv), r);
        }

        draw_on_single_reaction(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };

    struct draw_on_iv_escape
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt)
        {
            return add(
                domain.shell().second.position(),
                normalize(
                    create_vector<position_type>(
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.)),
                    draw_r(
                        rng_,
                        GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                        dt, domain.a_R())));
        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            boost::scoped_ptr<PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(domain.a_R());
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    array_gen(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt)
        {
            boost::shared_ptr<structure_type> const _structure(
                world_.get_structure(
                    world_.get_species(
                        domain.particles()[0].second.sid())
                    .structure_id()));
            
            cylindrical_surface_type const* const structure(
                dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            BOOST_ASSERT(structure);

            length_type const r_R(draw_r(
                rng_,
                GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                dt, domain.a_R()));
            return add(
                domain.shell().second.position(),
                multiply(structure->shape().unit_z(), r_R));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            return multiply(normalize(old_iv), domain.a_r());
        }

        draw_on_iv_escape(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };

    struct draw_on_iv_reaction
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt)
        {
            return add(
                domain.shell().second.position(),
                normalize(
                    create_vector<position_type>(
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.)),
                    draw_r(
                        rng_,
                        GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                        dt, domain.a_R())));
        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            boost::scoped_ptr<PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(domain.sigma());
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    array_gen(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt)
        {
            boost::shared_ptr<structure_type> const _structure(
                world_.get_structure(
                    world_.get_species(
                        domain.particles()[0].second.sid()).structure_id()));
            
            cylindrical_surface_type const* const structure(
                dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            BOOST_ASSERT(structure);

            length_type const r_R(draw_r(
                rng_,
                GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                dt, domain.a_R()));
            return add(
                domain.shell().second.position(),
                multiply(structure->shape().unit_z(), r_R));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            return multiply(domain.sigma(), normalize(old_iv));
        }

        draw_on_iv_reaction(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };

    struct draw_on_burst
    {
        position_type draw_com(spherical_pair_type const& domain,
                               time_type dt)
        {
            return add(
                domain.shell().second.position(),
                normalize(
                    create_vector<position_type>(
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.),
                        rng_.uniform(-1., 1.)),
                    draw_r(
                        rng_,
                        GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                        dt, domain.a_R())));
        }

        position_type draw_iv(spherical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            boost::scoped_ptr<PairGreensFunction> const gf(
                choose_pair_greens_function(domain, dt));
            length_type const r(draw_r(rng_, *gf, dt, domain.a_R(), domain.sigma()));
            length_type const theta(draw_theta(rng_, *gf, dt, r));
            return adjust_iv_with_old_iv(
                spherical_to_cartesian(
                    array_gen(r, theta, rng_.uniform(0., 1.) * 2 * M_PI)),
                old_iv);
        }

        position_type draw_com(cylindrical_pair_type const& domain,
                               time_type dt)
        {
            boost::shared_ptr<structure_type> const _structure(
                world_.get_structure(
                    world_.get_species(
                        domain.particles()[0].second.sid())
                    .structure_id()));

            cylindrical_surface_type const* const structure(
                dynamic_cast<cylindrical_surface_type*>(_structure.get()));

            BOOST_ASSERT(structure);

            length_type const r_R(draw_r(
                rng_,
                GreensFunction3DAbsSym(domain.D_R(), domain.a_R()),
                dt, domain.a_R()));
            return add(
                domain.shell().second.position(),
                multiply(structure->shape().unit_z(), r_R));
        }

        position_type draw_iv(cylindrical_pair_type const& domain,
                              time_type dt, position_type const& old_iv)
        {
            BOOST_ASSERT(::size(domain.reactions()) == 1);
            length_type const r(
                draw_r(rng_,
                    GreensFunction3DRadAbs(
                        domain.D_tot(),
                        domain.reactions()[0].k(), domain.r0(),
                        domain.sigma(), domain.a_r()),
                    dt, domain.a_r(), domain.sigma()));
            BOOST_ASSERT(r > domain.sigma() && r <= domain.a_r());
            return multiply(normalize(old_iv), r);
        }

        draw_on_burst(rng_type& rng, world_type const& world)
            : rng_(rng), world_(world) {}

        rng_type& rng_;
        world_type const& world_;
    };
public:
    typedef abstract_limited_generator<domain_id_pair> domain_id_pair_generator;

public:
    virtual ~EGFRDSimulator()
    {
        //std::for_each(domains_.begin(), domains_.end(),
        //    compose_unary(delete_ptr<domain_type>(),
        //                  select_second<typename domain_map::value_type>()));
    }

    EGFRDSimulator(boost::shared_ptr<world_type> world,
                   boost::shared_ptr<network_rules_type const> network_rules,
                   rng_type& rng, reaction_recorder_type& rrec,
                   int dissociation_retry_moves,
                   length_type const& user_max_shell_size =
                    std::numeric_limits<length_type>::infinity())
        : base_type(world, network_rules, rng, rrec),
          num_retries_(dissociation_retry_moves),
          user_max_shell_size_(user_max_shell_size),
          ssmat_((*world).world_size(), (*world).matrix_size()),
          csmat_((*world).world_size(), (*world).matrix_size()),
          smatm_(boost::fusion::pair<spherical_shell_type,
                                     spherical_shell_matrix_type&>(ssmat_),
                 boost::fusion::pair<cylindrical_shell_type,
                                     cylindrical_shell_matrix_type&>(csmat_)),
          single_shell_factor_(.1),
          multi_shell_factor_(.05),
          rejected_moves_(0), zero_step_count_(0) {}

    length_type const& user_max_shell_size() const
    {
        return user_max_shell_size_;
    }

    length_type max_shell_size() const
    {
        return std::min((*base_type::world_).matrix_size() * .5 /
                        traits_type::SAFETY,
                   user_max_shell_size_);
    }

    template<typename Tshape>
    std::pair<const shell_id_type,
              typename traits_type::template shell_generator<Tshape>::type>
    new_shell(domain_id_type const& did, Tshape const& shape)
    {
        typedef typename traits_type::template shell_generator<Tshape>::type shell_type;
        std::pair<const shell_id_type, shell_type> const retval(shidgen_(), shell_type(did, shape));
        boost::fusion::at_key<shell_type>(smatm_).update(retval);
        return retval;
    }

    template<typename T>
    std::pair<const shell_id_type, T> const& get_shell(shell_id_type const& id) const
    {
        return boost::fusion::at_key<T>(smatm_)[id];
    }

    boost::shared_ptr<domain_type> get_domain(domain_id_type const& id) const
    {
        typename domain_map::const_iterator i(domains_.find(id));

        if (i == domains_.end())
        {
            throw not_found(
                (boost::format("domain id #%s not found") % boost::lexical_cast<std::string>(id)).str());
        }

        return (*i).second;
    }

    domain_id_pair_generator* get_domains() const
    {
        return make_range_generator<domain_id_pair>(domains_);
    }

    typename domain_map::size_type num_domains() const
    {
        return domains_.size();
    }

    int num_domains_per_type(std::type_info const& type) const
    {
        return domain_count_per_type_[&type];
    }

    std::vector<domain_id_type>*
    get_neighbor_domains(particle_shape_type const& p)
    {
        typedef domain_collector<no_filter> collector_type;
        no_filter f;
        collector_type col((*base_type::world_), p, f);
        boost::fusion::for_each(smatm_, shell_collector_applier<collector_type>(col, p.position()));
        return col.neighbors.container().get();
    }

    std::vector<domain_id_type>*
    get_neighbor_domains(particle_shape_type const& p, domain_id_type const& ignore)
    {
        typedef domain_collector<one_id_filter> collector_type;
        one_id_filter f(ignore);
        collector_type col((*base_type::world_), p, f);
        boost::fusion::for_each(smatm_, shell_collector_applier<collector_type>(col, p.position()));
        return col.neighbors.container().get();
    }

    virtual void clear_volume(particle_shape_type const& p)
    {
        boost::scoped_ptr<std::vector<domain_id_type> > domains(
            get_neighbor_domains(p));
        if (domains)
        {
            burst_domains(*domains);
        }
    }

    virtual void clear_volume(particle_shape_type const& p,
                              domain_id_type const& ignore)
    {
        boost::scoped_ptr<std::vector<domain_id_type> > domains(
            get_neighbor_domains(p, ignore));

        if (domains)
        {
            burst_domains(*domains);
        }
    }

    virtual void initialize()
    {
        domains_.clear();
        ssmat_.clear();
        csmat_.clear();
    }

    virtual void step()
    {
        step(base_type::dt_);
    }

    virtual bool step(time_type const& upto)
    {
        time_type const lt(upto - base_type::t_);
        if (lt <= 0.)
            return false;
        if (base_type::dt_ < lt)
        {
            _step(base_type::dt_);
        }
        else
        {
            _step(lt);
            base_type::t_ = upto;
        }
        return true;
    }

protected:
    template<typename Tshell>
    void move_shell(std::pair<const shell_id_type, Tshell> const& shell)
    {
        typedef Tshell shell_type;
        boost::fusion::at_key<shell_type>(smatm_).update(shell);
    }

    void move_single_particle(single_type& domain,
                              position_type const& new_pos)
    {
        particle_type const& old(domain.particle().second);
        particle_id_pair new_pp(
            domain.particle().first,
            particle_type(old.sid(),
                particle_shape_type(new_pos, old.radius()), old.D()));
        domain.particle().second = new_pp.second;
        (*base_type::world_).update_particle(new_pp);
    }

    // move_domain {{{
    template<typename T>
    void move_domain(AnalyticalSingle<traits_type, T>& domain,
                     position_type const& new_pos,
                     length_type const& new_shell_size)
    {
        typedef typename AnalyticalSingle<traits_type, T>::shell_type shell_type;
        typename shell_type::shape_type new_shape(shape(domain.shell().second));
        shape_position(new_shape) = new_pos;
        shape_size(new_shape) = new_shell_size;
        shell_type const new_shell(domain.id(), new_shape);
        move_shell(std::pair<shell_id_type const, shell_type>(domain.shell().first, new_shell));
    }

    template<typename T>
    void move_domain(AnalyticalSingle<traits_type, T>& domain,
                     position_type const& new_pos)
    {
        move_domain(domain, new_pos, shape_size(shape(domain.shell().second)));
    }
    // }}}

    // remove_domain {{{
    template<typename T>
    void remove_domain(AnalyticalSingle<traits_type, T>& domain)
    {
        typedef T shell_type;
        LOG_DEBUG(("remove domain: %s", boost::lexical_cast<std::string>(domain.id()).c_str()));
        domains_.erase(domain.id());
        boost::fusion::at_key<shell_type>(smatm_).erase(domain.shell().first);
        remove_event(domain);
    }

    template<typename T>
    void remove_domain(AnalyticalPair<traits_type, T>& domain)
    {
        typedef T shell_type;
        LOG_DEBUG(("remove domain: %s", boost::lexical_cast<std::string>(domain.id()).c_str()));
        domains_.erase(domain.id());
        boost::fusion::at_key<shell_type>(smatm_).erase(domain.shell().first);
        remove_event(domain);
    }

    void remove_domain(multi_type& domain)
    {
        LOG_DEBUG(("remove domain: %s", boost::lexical_cast<std::string>(domain.id()).c_str()));
        domains_.erase(domain.id());
        typename multi_type::spherical_shell_id_pair_range r(domain.get_shells());
        spherical_shell_matrix_type& mat(
            boost::fusion::at_key<spherical_shell_type>(smatm_));
        std::for_each(boost::begin(r), boost::end(r),
            compose_unary(
                boost::bind(&spherical_shell_matrix_type::erase, mat, _1),
                select_first<spherical_shell_id_pair>()));
        remove_event(domain);
    }

    void remove_domain(domain_type& domain)
    {
        {
            spherical_single_type* _domain(dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        {
            cylindrical_single_type* _domain(dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        {
            spherical_pair_type* _domain(dynamic_cast<spherical_pair_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        {
            cylindrical_pair_type* _domain(dynamic_cast<cylindrical_pair_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        {
            multi_type* _domain(dynamic_cast<multi_type*>(&domain));
            if (_domain)
            {
                remove_domain(*_domain);
                return;
            }
        }
        throw not_implemented(std::string("unsupported domain type"));
    }
    // }}}

    void add_event(single_type& domain, single_event_kind const& kind)
    {
        boost::shared_ptr<event_type> new_event(
            new single_event(base_type::t_ + domain.dt(), domain, kind));
        domain.event() = std::make_pair(scheduler_.add(new_event), new_event);
        LOG_DEBUG(("add_event: #%d - %s", domain.event().first, boost::lexical_cast<std::string>(domain).c_str()));
    }

    void add_event(pair_type& domain, pair_event_kind const& kind)
    {
        boost::shared_ptr<event_type> new_event(
            new pair_event(base_type::t_ + domain.dt(), domain, kind));
        domain.event() = std::make_pair(scheduler_.add(new_event), new_event);
        LOG_DEBUG(("add_event: #%d - %s", domain.event().first, boost::lexical_cast<std::string>(domain).c_str()));
    }

    void add_event(multi_type& domain)
    {
        boost::shared_ptr<event_type> new_event(
            new multi_event(base_type::t_ + domain.dt(), domain));
        domain.event() = std::make_pair(scheduler_.add(new_event), new_event);
        LOG_DEBUG(("add_event: #%d - %s", domain.event().first, boost::lexical_cast<std::string>(domain).c_str()));
    }

    void update_event(single_type& domain, single_event_kind const& kind)
    {
        LOG_DEBUG(("update_event: #%d", domain.event().first));
        boost::shared_ptr<event_type> new_event(
            new single_event(base_type::t_ + domain.dt(), domain, kind));
        domain.event() = std::make_pair(domain.event().first, new_event);
        scheduler_.update(domain.event());
    }

    void update_event(pair_type& domain, pair_event_kind const& kind)
    {
        LOG_DEBUG(("update_event: #%d", domain.event().first));
        boost::shared_ptr<event_type> new_event(
            new pair_event(base_type::t_ + domain.dt(), domain, kind));
        domain.event() = std::make_pair(domain.event().first, new_event);
        scheduler_.update(domain.event());
    }

    void update_event(multi_type& domain)
    {
        LOG_DEBUG(("update_event: #%d", domain.event().first));
        boost::shared_ptr<event_type> new_event(
            new pair_event(base_type::t_ + domain.dt(), domain));
        domain.event() = std::make_pair(domain.event().first, new_event);
        scheduler_.update(domain.event());
    }

    void remove_event(domain_type& domain)
    {
        LOG_DEBUG(("remove_event: #%d", domain.event().first));
        scheduler_.remove(domain.event().first);
    }

    // create_single {{{
    boost::shared_ptr<single_type> create_single(particle_id_pair const& p)
    {
        single_type* new_single(0);
        domain_id_type did(didgen_());
        species_type const& species((*base_type::world_).get_species(p.second.sid()));
        boost::shared_ptr<structure_type> structure((*base_type::world_).get_structure(species.structure_id()));
        {
            planar_surface_type* _structure(
                dynamic_cast<planar_surface_type*>(structure.get()));
            if (_structure)
            {
                cylindrical_shell_id_pair const new_shell(
                    this->new_shell(did, cylinder_type(
                        p.second.position(), p.second.radius(),
                        normalize(cross_product(
                            _structure->shape().unit_x(),
                            _structure->shape().unit_y())),
                        p.second.radius())));
                new_single = new cylindrical_single_type(did, p, new_shell);
            }
        }
        {
            // Heads up. The cylinder's *size*, not radius, is changed when you 
            // make the cylinder bigger, because of the redefinition of set_radius.

            // The radius of a rod is not more than it has to be (namely the radius 
            // of the particle), so if the particle undergoes an unbinding reaction 
            // we still have to clear the target volume and the move may be 
            // rejected (NoSpace error).
            cylindrical_surface_type* _structure(
                dynamic_cast<cylindrical_surface_type*>(structure.get()));
            if (_structure)
            {
                const cylindrical_shell_id_pair new_shell(
                    this->new_shell(
                        did, cylinder_type(
                            p.second.position(), p.second.radius(),
                            _structure->shape().unit_z(),
                            p.second.radius())));
                new_single = new cylindrical_single_type(did, p, new_shell);
            }
        }
        {
            cuboidal_region_type* _structure(
                dynamic_cast<cuboidal_region_type*>(structure.get()));
            if (_structure)
            {
                spherical_shell_id_pair new_shell(
                    this->new_shell(did, ::shape(p.second)));
                new_single = new spherical_single_type(did, p, new_shell);
            }
        }
        if (!new_single)
        {
            throw not_implemented(
                (boost::format("unsupported structure type: %s") %
                    boost::lexical_cast<std::string>(*structure)).str());;
        }

        boost::shared_ptr<domain_type> const retval(new_single);
        domains_.insert(std::make_pair(did, retval));
        ++domain_count_per_type_[&typeid(*new_single)];
        return boost::dynamic_pointer_cast<single_type>(retval);
    }
    // }}}

    // create_pair {{{
    boost::shared_ptr<pair_type> create_pair(particle_id_pair const& p0,
                                             particle_id_pair const& p1,
                                             position_type const& com,
                                             position_type const& iv,
                                             length_type const& shell_size)
    {
        pair_type* new_pair(0);
        domain_id_type did(didgen_());
        species_type const& species((*base_type::world_).get_species(p0.second.sid()));
        boost::shared_ptr<structure_type> const structure(((*base_type::world_).get_structure(species.structure_id())));
        typename network_rules_type::reaction_rule_vector const& rules(
            (*base_type::network_rules_).query_reaction_rule(
                p0.second.sid(), p1.second.sid()));
        {
            planar_surface_type* _structure(
                dynamic_cast<planar_surface_type*>(structure.get()));
            if (_structure)
            {
                cylindrical_shell_id_pair const new_shell(
                    this->new_shell(did, cylinder_type(
                        com,
                        shell_size,
                        normalize(cross_product(
                            _structure->shape().unit_x(),
                            _structure->shape().unit_y())),
                        std::max(p0.second.radius(), p1.second.radius()))));
                new_pair = new cylindrical_pair_type(did, p0, p1, new_shell,
                                                       iv, rules);
            }
        }
        {
            // The radius of a rod is not more than it has to be (namely the 
            // radius of the biggest particle), so if the particle undergoes 
            // an unbinding reaction we still have to clear the target volume 
            // and the move may be rejected (NoSpace error).

            cylindrical_surface_type* _structure(
                dynamic_cast<cylindrical_surface_type*>(structure.get()));
            if (_structure)
            {
                cylindrical_shell_id_pair const new_shell(
                    this->new_shell(did, cylinder_type(
                        com,
                        shell_size,
                        shape(*_structure).unit_z(),
                        std::max(p0.second.radius(), p1.second.radius()))));
                new_pair = new cylindrical_pair_type(did, p0, p1, new_shell,
                                                     iv, rules);
            }
        }
        {
            cuboidal_region_type* _structure(
                dynamic_cast<cuboidal_region_type*>(structure.get()));
            if (_structure)
            {
                spherical_shell_id_pair new_shell(
                    this->new_shell(did,
                        sphere_type(com, shell_size)));
                new_pair = new spherical_pair_type(did, p0, p1, new_shell,
                                                   iv, rules);
            }
        }

        if (!new_pair)
        {
            throw not_implemented(
                (boost::format("unsupported structure type: %s") %
                    boost::lexical_cast<std::string>(*structure)).str());;
        }

        boost::shared_ptr<domain_type> const retval(new_pair);
        domains_.insert(std::make_pair(did, retval));
        ++domain_count_per_type_[&typeid(*new_pair)];
        return boost::dynamic_pointer_cast<pair_type>(retval);
    }
    // }}}

    // create_multi {{{
    boost::shared_ptr<multi_type> create_multi()
    {
        domain_id_type did(didgen_());
        multi_type* new_multi(new multi_type(did, *this, traits_type::DEFAULT_DT_FACTOR));
        boost::shared_ptr<domain_type> const retval(new_multi);
        domains_.insert(std::make_pair(did, retval));
        ++domain_count_per_type_[&typeid(*new_multi)];
        return boost::dynamic_pointer_cast<multi_type>(retval);
    }
    // }}}

    // draw_r {{{
    template<typename Tgf>
    static length_type draw_r(rng_type& rng,
                              Tgf const& gf,
                              time_type const& dt,
                              length_type const& a,
                              length_type const& sigma = -1.)
    {
        length_type r(0.);
        double rnd(0.);
        try
        {
            do
            {
                rnd = rng.uniform(0., 1.);
                r = gf.drawR(rnd, dt);
            } while (r > a || r <= sigma);
        }
        catch (std::exception const& e)
        {
            throw propagation_error(
                (boost::format(
                    "gf.drawR() failed: %s, rnd=%g, dt=%g, a=%g, sigma=%g") % e.what() % rnd % dt % a % sigma).str());
        }

        return r;
    }
    // }}}

    // draw_theta {{{
    template<typename Tgf>
    static length_type draw_theta(rng_type& rng,
                              Tgf const& gf,
                              time_type const& dt,
                              length_type const& r)
    {
        length_type theta(0.);
        double rnd(0.);
        try
        {
            rnd = rng.uniform(0., 1.);
            theta = gf.drawTheta(rnd, r, dt);
        }
        catch (std::exception const& e)
        {
            throw propagation_error(
                (boost::format(
                    "gf.drawTheta() failed: %s, rnd=%g, dt=%g, r=%g") % e.what() % rnd % dt % r).str());
        }

        return theta;
    }
    // }}}

    // draw_displacement {{{
    position_type draw_displacement(
        AnalyticalSingle<traits_type, spherical_shell_type> const& domain,
        length_type const& r)
    {
        return normalize(
            create_vector<position_type>(
                base_type::rng_.uniform(-1., 1.),
                base_type::rng_.uniform(-1., 1.),
                base_type::rng_.uniform(-1., 1.)), r);
    }

    position_type draw_displacement(
        AnalyticalSingle<traits_type, cylindrical_shell_type> const& domain,
        length_type const& r)
    {
        return multiply(shape(domain.shell().second).unit_z(), r);
    }
    // }}}

    // draw_new_position {{{
    template<typename Tshell>
    position_type draw_new_position(
            AnalyticalSingle<traits_type, Tshell> const& domain,
            time_type const& dt)
    {
        typedef Tshell shell_type;
        typedef typename shell_type::shape_type shape_type;
        typedef typename detail::get_greens_function<shape_type>::type greens_function;
        length_type const r(
            draw_r(
                base_type::rng_,
                greens_function(
                    domain.particle().second.D(),
                    domain.mobility_radius()),
                base_type::dt_,
                domain.mobility_radius()));
        position_type const displacement(draw_displacement(domain, r));
#ifdef DEBUG
        length_type const scale(domain.particle().second.radius());
        BOOST_ASSERT(!feq(length(displacement), std::abs(r), scale));
#endif
        return add(domain.particle().second.position(), displacement);
    }

    position_type draw_new_position(single_type& domain,
                                    time_type const& dt)
    {
        {
            spherical_single_type* _domain(dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
            {
                return draw_new_position(*_domain, dt);
            }
        }
        {
            cylindrical_single_type* _domain(dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
            {
                return draw_new_position(*_domain, dt);
            }
        }
        throw not_implemented(std::string("unsupported domain type"));
    }
    // }}}

    // draw_new_positions {{{
    template<typename Tdraw, typename T>
    boost::array<position_type, 2> draw_new_positions(
        AnalyticalPair<traits_type, T> const& domain,
        time_type const& dt)
    {
        Tdraw d(base_type::rng_, *base_type::world_);
        position_type const new_com(d.draw_com(domain, dt));
        position_type const new_iv(d.draw_iv(domain, dt, domain.iv()));
        D_type const D0(domain.particles()[0].second.D());
        D_type const D1(domain.particles()[1].second.D());
        return array_gen(
            subtract(new_com, multiply(new_iv, D0 / (D0 + D1))),
            add(new_com, multiply(new_iv, D1 / (D0 + D1))));
    }
    // }}}

    // propagate {{{
    /**
     * The difference between a burst and a propagate is that a burst 
     * always takes place before the actual scheduled event for the single, 
     * while propagate_single can be called for an escape event.
     *
     * Another subtle difference is that burst_single always reschedules 
     * (update_event) the single, while just calling propagate does not. 
     * So whoever calls propagate_single directly should reschedule the single 
     * afterwards.
     */
    //template<typename T>
    //void propagate(AnalyticalSingle<traits_type, T>& domain, position_type const& new_pos)
    void propagate(single_type& domain, position_type const& new_pos)
    {
        LOG_DEBUG(("single.dt=%g, single.last_time=%g, self.t=%g",
                   domain.dt(), domain.last_time(), base_type::t_));

        position_type const _new_pos(
            (*base_type::world_).apply_boundary(new_pos));
#ifdef DEBUG
        LOG_DEBUG(("propagate %s: %s => %s",
            boost::lexical_cast<std::string>(domain).c_str(),
            boost::lexical_cast<std::string>(domain.particle().second.position()).c_str(),
            boost::lexical_cast<std::string>(_new_pos).c_str()));

        particle_shape_type const new_particle(new_pos, domain.particle().second.radius());
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                (*base_type::world_).check_overlap(
                    new_particle, domain.particle().first));
        if (overlapped && ::size(*overlapped))
        {
            throw propagation_error("propagate_single: check_overlap failed.");
        }
#endif

        move_single_particle(domain, new_pos);
    }

    template<typename T>
    boost::array<boost::shared_ptr<single_type>, 2>
    propagate(AnalyticalPair<traits_type, T>& domain,
              boost::array<position_type, 2> const& _new_pos)
    {
        boost::array<particle_id_pair, 2> const& particles(domain.particles());
        boost::array<particle_id_pair, 2> new_particles(particles);
        new_particles[0].second.position() = (*base_type::world_).apply_boundary(_new_pos[0]);
        new_particles[1].second.position() = (*base_type::world_).apply_boundary(_new_pos[1]);

        BOOST_ASSERT(
            !(*base_type::world_).check_overlap(
                shape(new_particles[0].second),
                new_particles[0].first, new_particles[1].first));
        BOOST_ASSERT(
            !(*base_type::world_).check_overlap(
                shape(new_particles[1].second),
                new_particles[0].first, new_particles[1].first));
        BOOST_ASSERT(check_pair_pos(domain, new_particles));

        LOG_DEBUG(("propagate: %s => %s, %s => %s",
            boost::lexical_cast<std::string>(particles[0].second.position()).c_str(),
            boost::lexical_cast<std::string>(new_particles[0].second.position()).c_str(),
            boost::lexical_cast<std::string>(particles[1].second.position()).c_str(),
            boost::lexical_cast<std::string>(new_particles[1].second.position()).c_str()));

        remove_domain(domain);

        return array_gen(
            create_single(new_particles[0]),
            create_single(new_particles[1]));
    }
    // }}}


    template<typename Trange>
    void burst_domains(Trange const& domain_ids, boost::optional<std::vector<boost::shared_ptr<domain_type> >&> const& result = boost::optional<std::vector<boost::shared_ptr<domain_type> >&>())
    {
        BOOST_FOREACH(domain_id_type id, domain_ids)
        {
            boost::shared_ptr<domain_type> domain(get_domain(id));
            burst(domain, result);
        }
    }

    // burst {{{
    template<typename T>
    void burst(AnalyticalSingle<traits_type, T>& domain)
    {
        position_type const old_pos(shape_position(shape(domain.shell().second)));
        length_type const old_shell_size(shape_size(shape(domain.shell().second))); 
        length_type const particle_radius(domain.particle().second.radius());

        // Override dt, burst happens before single's scheduled event.
        domain.dt() = base_type::t_ - domain.last_time();

        position_type const new_pos(draw_new_position(domain, domain.dt()));

        propagate(domain, new_pos);

        move_domain(domain, new_pos);

        // Todo. if isinstance(single, InteractionSingle):
        domain.dt() = 0.;
        domain.last_time() = base_type::t_;
        update_event(domain, SINGLE_EVENT_ESCAPE);

        BOOST_ASSERT(
            (*base_type::world_).distance(new_pos, old_pos)
                <= old_shell_size - particle_radius);
        // Displacement check is in NonInteractionSingle.draw_new_position.

        BOOST_ASSERT(shape_size(shape(domain.shell().second)) == particle_radius);
    }

    template<typename T>
    boost::array<boost::shared_ptr<single_type>, 2> burst(AnalyticalPair<traits_type, T>& domain)
    {
        length_type const dt(base_type::t_ - domain.last_time());

        return propagate(domain,
            draw_new_positions<draw_on_burst>(domain, dt));
    }

    void burst(multi_type& domain, boost::optional<std::vector<boost::shared_ptr<domain_type> >&> const& result = boost::optional<std::vector<boost::shared_ptr<domain_type> >&>())
    {
        remove_domain(domain);

        BOOST_FOREACH(particle_id_pair p, domain.get_particles_range())
        {
            boost::shared_ptr<single_type> s(create_single(p));
            add_event(*s, SINGLE_EVENT_ESCAPE);
            if (result)
            {
                result.get().push_back(boost::dynamic_pointer_cast<domain_type>(s));
            }
        }
    }

    void burst(boost::shared_ptr<domain_type> domain, boost::optional<std::vector<boost::shared_ptr<domain_type> >&> const& result = boost::optional<std::vector<boost::shared_ptr<domain_type> >&>())
    {
        LOG_DEBUG(("burst: bursting %s", boost::lexical_cast<std::string>(*domain).c_str()));
        BOOST_ASSERT(base_type::t_ >= domain->last_time());
        BOOST_ASSERT(base_type::t_ <= domain->last_time() + domain->dt());

        {
            spherical_single_type* _domain(dynamic_cast<spherical_single_type*>(domain.get()));
            if (_domain)
            {
                burst(*_domain);
                if (result)
                    result.get().push_back(domain);
                return;
            }
        }
        {
            cylindrical_single_type* _domain(dynamic_cast<cylindrical_single_type*>(domain.get()));
            if (_domain)
            {
                burst(*_domain);
                if (result)
                    result.get().push_back(domain);
                return;
            }
        }
        {
            spherical_pair_type* _domain(dynamic_cast<spherical_pair_type*>(domain.get()));
            if (_domain)
            {
                boost::array<boost::shared_ptr<single_type>, 2> bursted(burst(*_domain));
                if (result)
                {
                    result.get().push_back(boost::dynamic_pointer_cast<domain_type>(bursted[0]));
                    result.get().push_back(boost::dynamic_pointer_cast<domain_type>(bursted[1]));
                }
                return;
            }
        }
        {
            cylindrical_pair_type* _domain(dynamic_cast<cylindrical_pair_type*>(domain.get()));
            if (_domain)
            {
                boost::array<boost::shared_ptr<single_type>, 2> bursted(burst(*_domain));
                if (result)
                {
                    result.get().push_back(boost::dynamic_pointer_cast<domain_type>(bursted[0]));
                    result.get().push_back(boost::dynamic_pointer_cast<domain_type>(bursted[1]));
                }
                return;
            }
        }
        {
            multi_type* _domain(dynamic_cast<multi_type*>(domain.get()));
            if (_domain)
            {
                burst(*_domain, result);
                return;
            }
        }
        throw not_implemented("?");
    }
    // }}}

    // attempt_single_reaction {{{
    bool attempt_single_reaction(single_type& domain)
    {
        const particle_id_pair reactant(domain.particle());
        const species_type reactant_species((*base_type::world_).get_species(reactant.second.sid()));
        reaction_rules const& rules((*base_type::network_rules_).query_reaction_rule(reactant.second.sid()));
        if (::size(rules) == 0)
        {
            return false;
        }

        reaction_rule_type const& r(draw_reaction_rule(rules));

        switch (::size(r.get_products()))
        {
        case 0:
            remove_domain(domain);
            (*base_type::world_).remove_particle(reactant.first);
            (base_type::rrec_)(reaction_record_type(r.id(), array_gen<particle_id_type>(), reactant.first));
            break;
        case 1: 
            {
                species_type const& product_species(
                    (*base_type::world_).get_species(r.get_products()[0]));

                if (reactant_species.radius() < product_species.radius())
                {
                    clear_volume(::shape(reactant.second));
                }

                if ((*base_type::world_).check_overlap(::shape(reactant.second), reactant.first))
                {
                    log_.info("no space for product particle.");
                    throw no_space();
                }

                remove_domain(domain);
                (*base_type::world_).remove_particle(reactant.first);
                particle_id_pair product(
                    (*base_type::world_).new_particle(
                        product_species.id(), reactant.second.position()));
                boost::shared_ptr<single_type> new_domain(create_single(product));
                add_event(*new_domain, SINGLE_EVENT_ESCAPE);
                (base_type::rrec_)(reaction_record_type(r.id(), array_gen(product.first), reactant.first));
            }
            break;
        case 2:
            {
                species_type const& product_species0(
                    (*base_type::world_).get_species(r.get_products()[0]));
                species_type const& product_species1(
                    (*base_type::world_).get_species(r.get_products()[1]));

                D_type const D01(product_species0.D() + product_species1.D());
                length_type r01(product_species0.radius() + product_species1.radius());
                Real const rad(std::max(
                        r01 * (product_species0.D() / D01) + product_species0.radius(),
                        r01 * (product_species1.D() / D01) + product_species1.radius()));
                clear_volume(particle_shape_type(reactant.second.position(), rad));

                particle_shape_type new0, new1;

                int i = num_retries_;
                while (--i >= 0)
                {
                    boost::shared_ptr<structure_type> structure(
                        (*base_type::world_).get_structure(
                            reactant_species.structure_id()));
                    position_type vector(
                        structure->random_vector(
                            r01 * traits_type::world_type::traits_type::MINIMAL_SEPARATION_FACTOR,
                            base_type::rng_));
                    // place particles according to the ratio D1:D2
                    // this way, species with D=0 doesn't move.
                    // FIXME: what if D1 == D2 == 0?
                    for (;;) {
                        new0 = particle_shape_type(
                            (*base_type::world_).apply_boundary(
                                add(reactant.second.position(),
                                    multiply(vector, product_species0.D() / D01))),
                            product_species0.radius());
                        new1 = particle_shape_type(
                            (*base_type::world_).apply_boundary(
                                add(reactant.second.position(),
                                    multiply(vector, product_species1.D() / D01))),
                            product_species1.radius());

                        if ((*base_type::world_).distance(new0.position(), new1.position()) >= r01)
                            break;

                        vector *= 1.0 + 1e-7;
                    }

                    // accept the new positions if there is enough space.
                    if ((!(*base_type::world_).check_overlap(
                            new0, reactant.first)) &&
                        (!(*base_type::world_).check_overlap(
                            new1, reactant.first)))
                        break;
                }
                if (i < 0)
                {
                    log_.info("no space for product particles.");
                    throw no_space();
                }

                remove_domain(domain);
                (*base_type::world_).remove_particle(reactant.first);

                const particle_id_pair
                    pp0(
                        (*base_type::world_).new_particle(
                            product_species0.id(), new0.position())),
                    pp1(
                        (*base_type::world_).new_particle(
                            product_species1.id(), new1.position()));
                // create domains for two particles and add them to
                // the event queue
                add_event(*create_single(pp0), SINGLE_EVENT_ESCAPE);
                add_event(*create_single(pp1), SINGLE_EVENT_ESCAPE);

                (base_type::rrec_)(reaction_record_type(
                    r.id(), array_gen(pp0.first, pp1.first), reactant.first));
            }
            break;
        default:
            throw not_implemented("reactions that produces more than two products are not supported.");
        }
        return true;
    }
    // }}}

    time_type draw_single_reaction_time(species_id_type const& sid)
    {
        reaction_rules const& rules(
            (*base_type::network_rules_).query_reaction_rule(sid));
        rate_type const k_tot(calculate_k_tot(rules));
        if (k_tot == 0.)
        {
            return std::numeric_limits<time_type>::infinity();
        }
        else if (k_tot == std::numeric_limits<rate_type>::infinity())
        {
            return 0.;
        }
        else
        {
            return (1. / k_tot) * std::log(1. / base_type::rng_.uniform(0., 1.));
        }
    }

    template<typename Tshell>
    time_type draw_escape_or_interaction_time(AnalyticalSingle<traits_type, Tshell> const& domain)
    {
        if (domain.particle().second.D() == 0.)
        {
            return std::numeric_limits<time_type>::infinity();
        }
        else
        {
            typedef Tshell shell_type;
            typedef typename shell_type::shape_type shape_type;
            typedef typename detail::get_greens_function<shape_type>::type greens_function;
            return greens_function(domain.particle().second.D(),
                            domain.mobility_radius())
                .drawTime(base_type::rng_.uniform(0., 1.));
        }
    }

    template<typename Tshell>
    std::pair<time_type, pair_event_kind>
    draw_com_escape_or_iv_event_time(AnalyticalPair<traits_type, Tshell> const& domain)
    {
        typedef Tshell shell_type;
        typedef typename shell_type::shape_type shape_type;
        typedef typename detail::get_pair_greens_function<shape_type> pair_greens_functions;
        typedef typename pair_greens_functions::iv_type iv_greens_function;
        typedef typename pair_greens_functions::com_type com_greens_function;
        BOOST_ASSERT(::size(domain.reactions()) == 1);
        time_type const dt_com(
            com_greens_function(domain.D_R(), domain.a_R()).drawTime());
        time_type const dt_iv(
            iv_greens_function(domain.D_tot(), domain.reactions()[0].k(),
                           domain.r0(), domain.sigma(), domain.a_r()).drawTime());
        if (dt_com < dt_iv)
        {
            return std::make_pair(dt_com, PAIR_EVENT_COM_ESCAPE);
        }
        else
        {
            return std::make_pair(dt_iv, PAIR_EVENT_IV);
        }
    }

    template<typename Tshell>
    std::pair<time_type, pair_event_kind>
    draw_single_reaction_time(AnalyticalPair<traits_type, Tshell> const& domain)
    {
        time_type const dt[2] = {
            draw_single_reaction_time(domain.particles()[0].second.sid()),
            draw_single_reaction_time(domain.particles()[1].second.sid())
        };
        if (dt[0] < dt[1])
        {
            return std::make_pair(dt[0], PAIR_EVENT_SINGLE_REACTION_0);
        }
        else
        {
            return std::make_pair(dt[1], PAIR_EVENT_SINGLE_REACTION_1);
        }
    }

    // {{{ determine_next_event
    template<typename Tshell>
    void determine_next_event(AnalyticalSingle<traits_type, Tshell>& domain)
    {
        typedef Tshell shell_type;
        typedef typename shell_type::shape_type shape_type;
        typedef typename detail::get_greens_function<shape_type>::type greens_function;
        time_type const dt_reaction(draw_single_reaction_time(domain.particle().second.sid()));
        time_type const dt_escape_or_interaction(draw_escape_or_interaction_time(domain));
        single_event_kind event_kind;
        if (dt_reaction < dt_escape_or_interaction)
        {
            domain.dt() = dt_reaction;
            event_kind = SINGLE_EVENT_REACTION;
        }
        else
        {
            domain.dt() = dt_escape_or_interaction;
            event_kind = SINGLE_EVENT_ESCAPE;
        }

        domain.last_time() = base_type::t_;
        update_event(domain, event_kind);
    }

    void determine_next_event(single_type& domain) const
    {
        {
            spherical_single_type* _domain(
                dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
            {
                determine_next_event(*_domain);
                return;
            }
        }
        {
            cylindrical_single_type* _domain(
                dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
            {
                determine_next_event(*_domain);
                return;
            }
        }
        throw not_implemented("unsupported domain type");
    }

    template<typename Tshell>
    void determine_next_event(AnalyticalPair<traits_type, Tshell>& domain)
    {
        std::pair<time_type, pair_event_kind> const dt_reaction(draw_single_reaction_time(domain));
        std::pair<time_type, pair_event_kind> const dt_com_escape_or_iv_event(draw_com_escape_or_iv_event_time(domain));
        std::pair<time_type, pair_event_kind> const dt_and_event_pair(
            dt_reaction.first < dt_com_escape_or_iv_event.first ?
                dt_reaction: dt_com_escape_or_iv_event);
        domain.dt() = dt_and_event_pair.first;
        domain.last_time() = base_type::t_;
        update_event(domain, dt_and_event_pair.second);
    }

    void determine_next_event(pair_type& domain) const
    {
        {
            spherical_pair_type* _domain(
                dynamic_cast<spherical_pair_type*>(&domain));
            if (_domain)
            {
                determine_next_event(*_domain);
                return;
            }
        }
        {
            cylindrical_pair_type* _domain(
                dynamic_cast<cylindrical_pair_type*>(&domain));
            if (_domain)
            {
                determine_next_event(*_domain);
                return;
            }
        }
        throw not_implemented("unsupported domain type");
    }
    // }}}

    // get_intruders {{{ 
    std::pair<std::vector<domain_id_type>*,
              std::pair<domain_id_type, length_type> >
    get_intruders(particle_shape_type const& p,
                  domain_id_type const& ignore) const
    {
        typedef intruder_collector<distance_calculator> collector_type;

        collector_type col(
            distance_calculator((*base_type::world_), p.position()), ignore);
        boost::fusion::for_each(smatm_, shell_collector_applier<collector_type>(col, p.position()));
        return std::make_pair(col.intruders.container().get(), col.closest);
    }
    // }}}

    template<typename TdidSet>
    std::pair<domain_id_type, length_type>
    get_closest_domain(position_type const& p, TdidSet const& ignore) const
    {
        typedef closest_object_finder<distance_calculator, TdidSet> collector_type;

        collector_type col(
            distance_calculator((*base_type::world_), p), ignore);
        boost::fusion::for_each(smatm_, shell_collector_applier<collector_type>(col, p));
        return col.closest;
    }

    void restore_domain(single_type& domain)
    {
        std::pair<domain_id_type, length_type> const closest(
            get_closest_domain(
                domain.position(), 
                array_gen(domain.id())));
        restore_domain(domain, closest);
    }

    template<typename T>
    void restore_domain(AnalyticalSingle<traits_type, T>& domain,
                        std::pair<domain_id_type, length_type> const& closest)
    {
        domain_type const& closest_domain(*get_domain(closest.first));
        length_type new_shell_size(0.);

        {
            single_type const* const _closest_domain(
                dynamic_cast<single_type const*>(&closest_domain));
            if (_closest_domain)
            {
                length_type const distance_to_closest(
                    (*base_type::world_).distance(
                        domain.position(), _closest_domain->position()));
                new_shell_size = calculate_single_shell_size(
                        domain, *_closest_domain,
                        distance_to_closest,
                        closest.second);
            } else {
                new_shell_size = closest.second / traits_type::SAFETY;
            }
        }
        LOG_DEBUG(("restore shell: %s (shell_size=%g, dt=%g) closest=%s (distance=%g)",
            boost::lexical_cast<std::string>(domain).c_str(),
            new_shell_size,
            domain.dt(),
            boost::lexical_cast<std::string>(closest_domain).c_str(),
            closest.second));
        move_domain(domain, domain.position(), new_shell_size);

        determine_next_event(domain);
    }

    void restore_domain(single_type& domain,
                        std::pair<domain_id_type, length_type> const& closest)
    {
        {
            spherical_single_type *_domain(
                dynamic_cast<spherical_single_type*>(&domain));
            if (_domain)
                return restore_domain(*_domain, closest);
        }
        {
            cylindrical_single_type *_domain(
                dynamic_cast<cylindrical_single_type*>(&domain));
            if (_domain)
                return restore_domain(*_domain, closest);
        }
        throw not_implemented(std::string("unsupported domain type"));
    }

    template<typename Trange>
    void burst_non_multis(Trange const& domain_ids,
                          std::vector<boost::shared_ptr<domain_type> >& bursted)
    {
        BOOST_FOREACH (domain_id_type id, domain_ids)
        {
            boost::shared_ptr<domain_type> domain(get_domain(id));
            if (dynamic_cast<multi_type*>(domain.get()))
            {
                bursted.push_back(domain);
            }
            else
            {
                burst(domain, bursted);
            }
        }
    }

    template<typename T>
    length_type distance(AnalyticalSingle<traits_type, T> const& domain,
                         position_type const& pos)
    {
        return (*base_type::world_).distance(shape(domain.shell().second), pos);
    }

    template<typename T>
    length_type distance(AnalyticalPair<traits_type, T> const& domain,
                         position_type const& pos)
    {
        return (*base_type::world_).distance(shape(domain.shell().second), pos);
    }

    length_type distance(multi_type const& domain, position_type const& pos)
    {
        length_type retval(std::numeric_limits<length_type>::infinity());
        BOOST_FOREACH (spherical_shell_id_pair const& shell,
                       domain.get_shells())
        {
            length_type const dist((*base_type::world_).distance(
                    shape(shell.second), pos));
            if (retval > dist)
            {
                retval = dist;
            }
        }
        return retval;
    }

    length_type distance(domain_type const& domain, position_type const& pos)
    {
        {
            spherical_single_type const* _domain(
                dynamic_cast<spherical_single_type const*>(&domain));
            if (_domain)
            {
                return distance(*_domain, pos);
            }
        }
        {
            cylindrical_single_type const* _domain(
                dynamic_cast<cylindrical_single_type const*>(&domain));
            if (_domain)
            {
                return distance(*_domain, pos);
            }
        }
        {
            spherical_pair_type const* _domain(
                dynamic_cast<spherical_pair_type const*>(&domain));
            if (_domain)
            {
                return distance(*_domain, pos);
            }
        }
        {
            cylindrical_pair_type const* _domain(
                dynamic_cast<cylindrical_pair_type const*>(&domain));
            if (_domain)
            {
                return distance(*_domain, pos);
            }
        }
        {
            multi_type const* _domain(dynamic_cast<multi_type const*>(&domain));
            if (_domain)
            {
                return distance(*_domain, pos);
            }
        }
        throw not_implemented(std::string("unsupported domain type"));
    }

    boost::optional<pair_type&>
    form_pair(single_type& domain, single_type& possible_partner,
              std::vector<boost::shared_ptr<domain_type> > const& neighbors)
    {
        LOG_DEBUG(("trying to form Pair(%s, %s)",
                    boost::lexical_cast<std::string>(domain).c_str(),
                    boost::lexical_cast<std::string>(possible_partner).c_str()));
        // 1. Determine min shell size.
        length_type const r[] = {
            domain.particle().second.radius(),
            possible_partner.particle().second.radius()
        };
        length_type const sigma(r[0] + r[1]);

        D_type const D[] = {
            domain.particle().second.D(),
            possible_partner.particle().second.D()
        };
        D_type const D01(D[0] + D[1]);

        BOOST_ASSERT(domain.particle().second.position() ==
                     domain.position());

        BOOST_ASSERT(possible_partner.particle().second.position() ==
                     possible_partner.position());

        position_type iv(
                subtract(domain.position(),
                    (*base_type::world_).cyclic_transpose(
                        possible_partner.position(),
                        domain.position())));
        length_type const r0(length(iv));
        length_type const distance_from_sigma(r0 - sigma);
        BOOST_ASSERT(distance_from_sigma >= 0);

        length_type const shell_size[] = {
           r0 * D[0] / D01 + r[0], r0 * D[1] / D01 + r[1]
        };
        length_type const shell_size_margin[] = {
            r[0] * 2,
            r[1] * 2
        };
        size_t const larger_shell_index(
            shell_size[0] + shell_size_margin[0] >=
                shell_size[1] + shell_size_margin[1] ? 0: 1);
        length_type const min_shell_size(shell_size[larger_shell_index]);
        length_type const min_shell_size_margin(shell_size_margin[larger_shell_index]);

        // 2. Check if min shell size not larger than max shell size or 
        // sim cell size.
        position_type com((*base_type::world_).apply_boundary(
            (*base_type::world_).calculate_pair_CoM(
                domain.position(), possible_partner.position(),
                D[0], D[1])));
        length_type const min_shell_size_with_margin(
            min_shell_size + min_shell_size_margin);
        length_type const max_shell_size(
            std::min(this->max_shell_size(),
                     distance_from_sigma * 100
                        + sigma + min_shell_size_margin));

        if (min_shell_size_with_margin >= max_shell_size)
        {
            LOG_DEBUG(("Pair(%s, %s) not formed: min_shell_size %g >="
                       "max_shell_size %g",
                       boost::lexical_cast<std::string>(domain).c_str(),
                       boost::lexical_cast<std::string>(possible_partner).c_str(),
                       min_shell_size_with_margin, max_shell_size));
            return boost::optional<pair_type&>();
        }

        // 3. Check if bursted Singles not too close.
        // The simple check for closest below could miss
        // some of them, because sizes of these Singles for this
        // distance check has to include SINGLE_SHELL_FACTOR, while
        // these burst objects have zero mobility radii.  This is not
        // beautiful, a cleaner framework may be possible.

        domain_type* closest_domain (0);
        length_type closest_shell_distance(std::numeric_limits<length_type>::infinity());
        BOOST_FOREACH (boost::shared_ptr<domain_type> _neighbor, neighbors)
        {
            single_type* const neighbor(
                dynamic_cast<single_type*>(_neighbor.get()));
            if (neighbor && neighbor->id() != possible_partner.id())
            {
                length_type const shell_distance(
                    (*base_type::world_).distance(com, neighbor->position()) -
                        neighbor->particle().second.radius() *
                            (1.0 + traits_type::SINGLE_SHELL_FACTOR));
                if (shell_distance < closest_shell_distance)
                {
                    closest_domain = neighbor;
                    closest_shell_distance = shell_distance;
                }
            }
        }

        BOOST_ASSERT(closest_domain);
        if (closest_shell_distance <= min_shell_size_with_margin)
        {
            LOG_DEBUG(("Pair(%s, %s) not formed: squeezed by burst neighbor %s",
                       boost::lexical_cast<std::string>(domain).c_str(),
                       boost::lexical_cast<std::string>(possible_partner).c_str(),
                       boost::lexical_cast<std::string>(*closest_domain).c_str()));
            return boost::optional<pair_type&>();
        }

        BOOST_ASSERT(closest_shell_distance > 0);

        // 4. Determine shell size and check if closest object not too 
        // close (squeezing).
        {
            std::pair<domain_id_type, length_type> possible_closest(
                get_closest_domain(com, array_gen(domain.id(),
                                                  possible_partner.id())));
            if (possible_closest.second < closest_shell_distance)
            {
                domain_type* const _closest_domain(
                        get_domain(possible_closest.first).get());
                closest_domain = _closest_domain;
                closest_shell_distance = possible_closest.second;
            }
        }

        BOOST_ASSERT(closest_domain);
        BOOST_ASSERT(closest_shell_distance > 0);

        LOG_DEBUG(("Pair closest neighbor: %s %g, "
                   "min_shell_with_margin=%g",
                   boost::lexical_cast<std::string>(*closest_domain).c_str(),
                   closest_shell_distance,
                   min_shell_size_with_margin));

        length_type new_shell_size(0.);

        {
            single_type* const _closest_domain(
                    dynamic_cast<single_type*>(closest_domain));
            if (_closest_domain)
            {
                particle_type const& closest_domain_particle(
                        _closest_domain->particle().second);
                D_type const D_tot(closest_domain_particle.D() + D01);
                length_type const closest_particle_distance(
                    (*base_type::world_).distance(
                        com, closest_domain_particle.position()));
                length_type const closest_min_shell(
                        closest_domain_particle.radius() *
                            (traits_type::SINGLE_SHELL_FACTOR + 1.0));

                // options for shell size:
                // a. ideal shell size
                // b. closest shell is from a bursted single
                // c. closest shell is closer than ideal shell size 
                new_shell_size = std::min(
                    std::min(
                        (D01 / D_tot) * (
                            closest_particle_distance - min_shell_size 
                            - closest_domain_particle.radius())
                        + min_shell_size,
                        closest_particle_distance - closest_min_shell),
                    closest_shell_distance);
            }
            else
            {
                new_shell_size = closest_shell_distance;
            }
            new_shell_size /= traits_type::SAFETY;
            BOOST_ASSERT(new_shell_size < closest_shell_distance);

            if (new_shell_size <= min_shell_size_with_margin)
            {
                LOG_DEBUG(("Pair(%s, %s) not formed: squeezed by %s",
                    boost::lexical_cast<std::string>(domain).c_str(),
                    boost::lexical_cast<std::string>(possible_partner).c_str(),
                    boost::lexical_cast<std::string>(closest_domain).c_str()));
                return boost::optional<pair_type&>();
            }
        }

        // 5. Check if singles would not be better.
        {
            length_type const dist[] = {
                (*base_type::world_).distance(com, domain.position()),
                (*base_type::world_).distance(com, domain.position())
            };

            if (new_shell_size < std::max(
                    dist[0] + r[0] *
                        (1.0 + traits_type::SINGLE_SHELL_FACTOR),
                    dist[1] + r[1] *
                        (1.0 + traits_type::SINGLE_SHELL_FACTOR)) * 1.3)
            {
                LOG_DEBUG(("Pair(%s, %s) not formed: leaving singles are better",
                            boost::lexical_cast<std::string>(domain).c_str(),
                            boost::lexical_cast<std::string>(possible_partner).c_str()));
                return boost::optional<pair_type&>();
            }
        }

        // 6. Ok, Pair makes sense. Create one.
        new_shell_size = std::min(new_shell_size, max_shell_size);

        boost::shared_ptr<pair_type> new_pair(
            create_pair(
                domain.particle(),
                possible_partner.particle(),
                com, iv, new_shell_size));

        determine_next_event(*new_pair);
        BOOST_ASSERT(new_pair->dt() >= 0);

        new_pair->last_time() = base_type::t_;

        remove_domain(domain);
        remove_domain(possible_partner);

        BOOST_ASSERT(
                closest_shell_distance ==
                    std::numeric_limits<length_type>::infinity()
                || new_shell_size < closest_shell_distance);
        BOOST_ASSERT(new_shell_size >= min_shell_size_with_margin);
        BOOST_ASSERT(new_shell_size <= max_shell_size);

        log_.info("%s,\ndt=%g, shell_size=%g, "
                  "closest_shell_distance=%s,\nclosest=%s",
                  boost::lexical_cast<std::string>(*new_pair).c_str(),
                  new_pair->dt(), new_shell_size,
                  closest_shell_distance,
                  boost::lexical_cast<std::string>(closest_domain).c_str());

        check_domain(*new_pair);

        return *new_pair;
    }

    boost::optional<multi_type&>
    form_multi(single_type& domain,
               std::vector<boost::shared_ptr<domain_type> > const& neighbors,
               std::pair<domain_type&, length_type> closest)
    {
        length_type const min_shell_size(
                domain.particle().second.radius() *
                    (1.0 + multi_shell_factor_));

        // Multis shells need to be contiguous.
        if (closest.second > min_shell_size)
        {
            return boost::optional<multi_type&>();
        }

        multi_type* retval(0);
        retval = dynamic_cast<multi_type*>(&closest.first);
        if (!retval)
        {
            retval = create_multi().get();
        }

        add_to_multi(*retval, domain);

        BOOST_FOREACH (boost::shared_ptr<domain_type> neighbor, neighbors)
        {
            length_type const dist(distance(*neighbor, domain.position()));
            if (dist < min_shell_size)
                add_to_multi_recursive(*retval, domain); 
        }

        return *retval;
    }

    bool add_to_multi(multi_type& multi, single_type& single)
    {
        LOG_DEBUG(("add to multi: %s => %s",
                boost::lexical_cast<std::string>(single).c_str(),
                boost::lexical_cast<std::string>(multi).c_str()));

        if (!multi.add_particle(single.particle()))
            return false;

        spherical_shell_id_pair sid_shell_pair(
            new_shell(
                multi.id(),
                sphere_type(
                    single.particle().second.position(),
                    single.particle().second.radius() *
                        (1. + multi_shell_factor_))));
        multi.add_shell(sid_shell_pair);
        remove_domain(single);

        return true;
    }

    void add_to_multi(multi_type& multi, multi_type& other_multi)
    {
        LOG_DEBUG(("add to multi: %s => %s",
                boost::lexical_cast<std::string>(other_multi).c_str(),
                boost::lexical_cast<std::string>(multi).c_str()));
        
        // merge other_multi into multi. other_multi will be removed.
        spherical_shell_matrix_type& mat(
            boost::fusion::at_key<spherical_shell_type>(smatm_));
        BOOST_FOREACH (spherical_shell_id_pair const& _shell,
                       other_multi.get_shells())
        {
            typename spherical_shell_matrix_type::iterator const i(
                mat.find(_shell.first));
            BOOST_ASSERT(i != mat.end());
            spherical_shell_type& shell((*i).second);
            shell.did() = multi.id();
            multi.add_shell(spherical_shell_id_pair(_shell.first, shell));
        }

        BOOST_FOREACH (particle_id_pair const& particle,
                       other_multi.get_particles_range())
        {
            multi.add_particle(particle);
        }

        domains_.erase(other_multi.id());
        remove_event(other_multi);
    }

    void add_to_multi_recursive(multi_type& multi, domain_type& domain)
    {
        {
            single_type* single(dynamic_cast<single_type*>(&domain));
            if (single)
            {
                particle_shape_type const new_shell(
                    single->particle().second.position(),
                    single->particle().second.radius() *
                        (1.0 + multi_shell_factor_));

                boost::scoped_ptr<std::vector<domain_id_type> > neighbors(
                    get_neighbor_domains(new_shell, single->id()));

                std::vector<boost::shared_ptr<domain_type> > bursted;
                burst_non_multis(*neighbors, bursted);

                BOOST_FOREACH (domain_id_type neighbor_id, *neighbors)
                {
                    boost::shared_ptr<domain_type> neighbor(get_domain(neighbor_id));
                    length_type const dist(distance(*neighbor, single->position()));
                    if (dist < new_shell.radius())
                        add_to_multi_recursive(multi, domain); 
                }
                return;
            }
        }
        {
            multi_type* other_multi(dynamic_cast<multi_type*>(&domain));
            if (other_multi)
            {
                add_to_multi(multi, *other_multi);
            }
        }
    }

    boost::optional<domain_type&> form_pair_or_multi(
        single_type& domain,
        std::vector<boost::shared_ptr<domain_type> > const& neighbors)
    {
        BOOST_ASSERT(!neighbors.empty());

        domain_type* possible_partner(0);
        length_type length_to_possible_partner(
                std::numeric_limits<length_type>::infinity());
        BOOST_FOREACH (boost::shared_ptr<domain_type> neighbor, neighbors)
        {
            length_type const dist(distance(*neighbor, domain.position()));
            if (dist < length_to_possible_partner)
            {
                possible_partner = neighbor.get();
                length_to_possible_partner = dist;
            }
        }

        // First, try forming a Pair.
        {
            single_type* const _possible_partner(
                    dynamic_cast<single_type*>(possible_partner));
            if (_possible_partner)
            {
                boost::optional<pair_type&> new_pair(
                    form_pair(domain, *_possible_partner, neighbors));
                if (new_pair)
                {
                    return new_pair.get();
                }
            }
        }

        // If a Pair is not formed, then try forming a Multi.
        {
            boost::optional<multi_type&> new_multi(
                    form_multi(domain, neighbors,
                               std::pair<domain_type&, length_type>(
                                    *possible_partner,
                                    length_to_possible_partner)));
            if (new_multi)
            {
                return new_multi.get();
            }
        }
        return boost::optional<domain_type&>();
    }

    void fire_event(single_event const& event)
    {
        single_type& domain(event.domain());
        BOOST_ASSERT(
            (domain.dt() + domain.last_time() - base_type::t_)
                <= 1e-18 * base_type::t_);
        ++single_step_count_[event.kind()];
        switch (event.kind())
        {
        case SINGLE_EVENT_REACTION:
            LOG_DEBUG(("fire_single: single reaction (%s)", boost::lexical_cast<std::string>(domain.id()).c_str()));
            propagate(domain, draw_new_position(domain, domain.dt()));
            try
            {
                attempt_single_reaction(domain);
            }
            catch (no_space const&)
            {
                LOG_DEBUG(("single reaction rejected"));
            }
            break;

        case SINGLE_EVENT_ESCAPE:
            LOG_DEBUG(("fire_single: single escape (%s)", boost::lexical_cast<std::string>(domain.id()).c_str()));

            // handle immobile case
            if (domain.D() == 0.)
            {
                determine_next_event(domain);
                domain.last_time() = base_type::t_;
                return;
            }

            if (domain.dt() != 0.)
                propagate(domain, draw_new_position(domain, domain.dt()));
            length_type const min_shell_radius(domain.particle().second.radius() * (1. + single_shell_factor_));
            {
                std::vector<domain_id_type>* intruders;
                std::pair<domain_id_type, length_type> closest;
                boost::tie(intruders, closest) = get_intruders(
                    particle_shape_type(domain.position(), min_shell_radius),
                    domain.id());
                boost::scoped_ptr<std::vector<domain_id_type> > _(intruders);

                LOG_DEBUG(("intruders: %s, closest: %s (dist=%g)",
                    intruders ?
                        boost::algorithm::join(
                            make_transform_iterator_range(
                                *intruders,
                                boost::bind(
                                    &boost::lexical_cast<std::string,
                                                         domain_id_type>, _1)),
                            std::string(", ")).c_str():
                        "(none)",
                    boost::lexical_cast<std::string>(closest.first).c_str(),
                    closest.second));
                if (intruders)
                {
                    std::vector<boost::shared_ptr<domain_type> > bursted;
                    burst_non_multis(*intruders, bursted);
                    if (form_pair_or_multi(domain, bursted))
                        return;
                    // if nothing was formed, recheck closest and restore shells.
                    restore_domain(domain);
                    BOOST_FOREACH (boost::shared_ptr<domain_type> _single, bursted)
                    {
                        boost::shared_ptr<single_type> single(
                            boost::dynamic_pointer_cast<single_type>(_single));
                        if (!single)
                            continue;
                        restore_domain(*single);
                    }
                } else {
                    restore_domain(domain, closest);
                }
                LOG_DEBUG(("%s (dt=%g)",
                    boost::lexical_cast<std::string>(domain).c_str(),
                    domain.dt()));
                add_event(domain, SINGLE_EVENT_ESCAPE);
            }
        }
    }

    template<typename Tshell>
    GreensFunction3DRadAbs::EventKind
    draw_iv_event_type(AnalyticalPair<traits_type, Tshell> const& domain)
    {
        typedef Tshell shell_type;
        typedef typename shell_type::shape_type shape_type;
        typedef typename detail::get_pair_greens_function<shape_type>::iv_type iv_greens_function;
        // Draw actual pair event for iv at very last minute.
        BOOST_ASSERT(::size(domain.reactions()) == 1);
        reaction_rule_type const& r(domain.reactions()[0]);
        iv_greens_function const gf(domain.D_tot(), r.k(), domain.r0(), domain.sigma(), domain.a_r());

        double const rnd(base_type::rng_.uniform(0, 1.));
        return gf.drawEventType(rnd, domain.dt());
    }

    void fire_event(pair_event const& event)
    {
        {
            spherical_pair_type* _domain(dynamic_cast<spherical_pair_type*>(&event.domain()));
            if (_domain)
            {
                fire_event(*_domain, event.kind());
                return;
            }
        }
        {
            cylindrical_pair_type* _domain(dynamic_cast<cylindrical_pair_type*>(&event.domain()));
            if (_domain)
            {
                fire_event(*_domain, event.kind());
                return;
            }
        }
    }

    template<typename T>
    void fire_event(AnalyticalPair<traits_type, T>& domain, pair_event_kind kind)
    {
        typedef AnalyticalSingle<traits_type, T> corresponding_single_type;
        check_domain(domain);

        ++pair_step_count_[kind];
        LOG_DEBUG(("fire_pair: event_kind %d", kind));

        //  1. Single reaction
        //  2. Pair reaction
        //  3a. IV escape
        //  3b. com escape

        switch (kind)
        {
        case PAIR_EVENT_SINGLE_REACTION_0: 
        case PAIR_EVENT_SINGLE_REACTION_1:
            {
                int const index(kind == PAIR_EVENT_SINGLE_REACTION_0 ? 0 : 1);
                int const theother_index(1 - index);
                particle_id_pair const& reacting_particle(
                        domain.particles()[index]);
                position_type const old_CoM(domain.position());
                LOG_DEBUG(("pair: single reaction %s", boost::lexical_cast<std::string>(reacting_particle.first).c_str()));

                boost::array<boost::shared_ptr<single_type>, 2> const new_single(burst(domain));

                add_event(*new_single[theother_index], SINGLE_EVENT_ESCAPE);
                try
                {
                    attempt_single_reaction(*new_single[index]);
                }
                catch (no_space const&)
                {
                    LOG_DEBUG(("single reaction rejected"));
                }
            }
            break;

        case PAIR_EVENT_COM_ESCAPE:
            {
                time_type const dt(domain.dt());
                boost::array<position_type, 2> const new_pos(
                    draw_new_positions<draw_on_com_escape>(
                        domain, dt));
                boost::array<boost::shared_ptr<single_type>, 2> const new_single(
                    propagate(domain, new_pos));

                add_event(*new_single[0], SINGLE_EVENT_ESCAPE);
                add_event(*new_single[1], SINGLE_EVENT_ESCAPE);
            }
            break;
        
        case PAIR_EVENT_IV:
            {
                // Draw actual pair event for iv at very last minute.
                switch (draw_iv_event_type(domain))
                {
                case GreensFunction3DRadAbs::IV_REACTION:
                    {
                        LOG_DEBUG(("iv_reaction"));
                        BOOST_ASSERT(::size(domain.reactions()) == 1);
                        reaction_rule_type const& r(domain.reactions()[0]);

                        switch (::size(r.get_products()))
                        {
                        case 1:
                            {
                                species_type const& new_species(
                                    (*base_type::world_).get_species(
                                        r.get_products()[0]));

                                // calculate new R
                                position_type const new_com(
                                    (*base_type::world_).apply_boundary(
                                        draw_on_iv_reaction(
                                            base_type::rng_,
                                            *base_type::world_).draw_com(
                                                domain, domain.dt())));
                           
                                BOOST_ASSERT(
                                    (*base_type::world_).distance(
                                        domain.shell().second.position(),
                                        new_com) + new_species.radius()
                                    < shape(domain.shell().second).radius());

                                (*base_type::world_).remove_particle(domain.particles()[0].first);
                                (*base_type::world_).remove_particle(domain.particles()[1].first);

                                particle_id_pair const new_particle(
                                    (*base_type::world_).new_particle(
                                        new_species.id(), new_com));
                                boost::shared_ptr<single_type> new_single(
                                    create_single(new_particle));
                                add_event(*new_single, SINGLE_EVENT_ESCAPE);

                                (base_type::rrec_)(reaction_record_type(
                                    r.id(),
                                    array_gen(new_particle.first),
                                    domain.particles()[0].first,
                                    domain.particles()[1].first));

                                LOG_DEBUG(("product: %s", boost::lexical_cast<std::string>(new_single).c_str()));
                            }
                            break;
                        default:
                            throw not_implemented("num products >= 2 not supported.");
                        }
                        remove_domain(domain);
                    }
                    break;
                case GreensFunction3DRadAbs::IV_ESCAPE:
                    {
                        time_type const dt(domain.dt());
                        boost::array<position_type, 2> const new_pos(
                            draw_new_positions<draw_on_iv_escape>(
                                domain, dt));
                        boost::array<boost::shared_ptr<single_type>, 2> const new_single(
                            propagate(domain, new_pos));

                        add_event(*new_single[0], SINGLE_EVENT_ESCAPE);
                        add_event(*new_single[1], SINGLE_EVENT_ESCAPE);
                    }
                    break;
                }
            }
            break;
        }
    }

    void fire_event(multi_event& event)
    {
        multi_type& domain(event.domain());
        domain.step();
        LOG_DEBUG(("fire_multi: %s", domain.last_event()));
        multi_step_count_[domain.last_event()]++; 
        switch (domain.last_event())
        {
        case multi_type::REACTION:
            (base_type::rrec_)(domain.last_reaction());
            burst(domain);
            break;
        case multi_type::ESCAPE:
            burst(domain);
            break;
        case multi_type::NONE:
            add_event(domain);
            break;
        }
    }

    void fire_event(event_type& event)
    {
        {
            single_event* _event(dynamic_cast<single_event*>(&event));
            if (_event)
            {
                fire_event(*_event);
                return;
            }
        }
        {
            pair_event* _event(dynamic_cast<pair_event*>(&event));
            if (_event)
            {
                fire_event(*_event);
                return;
            }
        }
        {
            multi_event* _event(dynamic_cast<multi_event*>(&event));
            if (_event)
            {
                fire_event(*_event);
                return;
            }
        }
        throw not_implemented(std::string("unsupported domain type"));
    }

    void _step(time_type const& dt)
    {
#ifdef DEBUG
        check();
#endif
        ++base_type::num_steps_;

        event_id_pair_type ev(scheduler_.pop());
        base_type::t_ = ev.second->time();

        log_.info("%d: t=%g dt=%g event=%s rejectedmoves=%d",
                  base_type::num_steps_, base_type::t_, base_type::dt_,
                  boost::lexical_cast<std::string>(ev.second).c_str(),
                  rejected_moves_);

        fire_event(*ev.second);

        time_type const next_time(scheduler_.top().second->time());
        base_type::dt_ = next_time - base_type::t_; 

        if (base_type::dt_ == 0.)
        {
            ++zero_step_count_;
            if (zero_step_count_ >= std::max(scheduler_.size(), static_cast<std::size_t>(10u)))
            {
                throw illegal_state("too many dt=zero steps. simulator halted?");
            }
        }
        else
        {
            zero_step_count_ = 0;
        }
    }

    void check_shell_matrix() const
    {
        typedef std::map<domain_id_type, std::vector<shell_id_type> >
            domain_shell_association;

        domain_shell_association did_map;

        boost::fusion::for_each(smatm_,
                domain_shell_map_builder<domain_shell_association>(
                    (*base_type::world_), did_map));

        typename domain_type::size_type shell_population(0);

        BOOST_FOREACH (typename event_scheduler_type::value_type const& value,
                       scheduler_.events())
        {
            domain_type const& domain(dynamic_cast<domain_event_base&>(*value.second).domain());
            typename domain_type::size_type const num_shells(domain.num_shells());

            shell_population += num_shells;
            typename std::vector<shell_id_type> const& shell_ids(
                did_map[domain.id()]);
            BOOST_ASSERT(static_cast<typename domain_type::size_type>(
                    ::size(shell_ids)) == num_shells);
        }

        std::size_t matrix_population(0);
        boost::fusion::for_each(smatm_,
                shell_population_summer(matrix_population));
        BOOST_ASSERT(shell_population == matrix_population);
    }

    void check_domains() const
    {
        std::set<event_id_type> event_ids;
        BOOST_FOREACH (typename domain_map::value_type const& domain_id_pair,
                       domains_)
        {
            event_ids.insert(domain_id_pair.second->event().first);
        }

        BOOST_FOREACH (event_id_pair_type const& event_id_pair,
                       scheduler_.events())
        {
            BOOST_ASSERT(event_ids.find(event_id_pair.first) != event_ids.end());
            event_ids.erase(event_id_pair.first);
        }

        // self.domains always include a None  --> this can change in future
        BOOST_ASSERT(event_ids.empty());
    }

    void check_every_domain() const
    {
        BOOST_FOREACH (event_id_pair_type const& event, scheduler_.events())
        {
            domain_event_base& _event(
                dynamic_cast<domain_event_base&>(*event.second));
            check_domain(_event.domain());
        }
    }

    void check_event_stoichiometry() const
    {
        std::size_t event_population(0);
        BOOST_FOREACH (event_id_pair_type const& event, scheduler_.events())
        {
            domain_event_base& _event(
                dynamic_cast<domain_event_base&>(*event.second));
            event_population += _event.domain().multiplicity();
        }

        BOOST_ASSERT((*base_type::world_).num_particles() == event_population);
    }

    void check() const
    {
        BOOST_ASSERT(base_type::t_ >= 0.0);
        BOOST_ASSERT(base_type::dt_ >= 0.0);

        check_shell_matrix();
        check_domains();
        check_event_stoichiometry();
        check_every_domain(); 
    }

    template<typename T>
    void check_domain(AnalyticalSingle<traits_type, T> const& domain) const
    {
        std::pair<domain_id_type, length_type> closest(
            get_closest_domain(domain.position(), array_gen(domain.id())));
        BOOST_ASSERT(shape_size(shape(domain.shell().second)) <= user_max_shell_size_);
        BOOST_ASSERT(shape_size(shape(domain.shell().second)) <= max_shell_size());
        BOOST_ASSERT(closest.second > shape_size(shape(domain.shell().second)));
    }

    template<typename T>
    void check_domain(AnalyticalPair<traits_type, T> const& domain) const
    {
        std::pair<domain_id_type, length_type> closest(
            get_closest_domain(domain.position(), array_gen(domain.id())));
        BOOST_ASSERT(shape_size(shape(domain.shell().second)) <= user_max_shell_size_);
        BOOST_ASSERT(shape_size(shape(domain.shell().second)) <= max_shell_size());
        BOOST_ASSERT(closest.second > shape_size(shape(domain.shell().second)));
    }

    void check_domain(multi_type const& domain) const
    {
        BOOST_FOREACH (typename multi_type::shell_id_pair const& shell,
                       domain.get_shells())
        {
            std::pair<domain_id_type, length_type> closest(
                get_closest_domain(shape_position(shape(shell.second)),
                                   array_gen(domain.id())));
            BOOST_ASSERT(shape_size(shape(shell.second)) <= user_max_shell_size_);
            BOOST_ASSERT(shape_size(shape(shell.second)) <= max_shell_size());
            BOOST_ASSERT(closest.second > shape_size(shape(shell.second)));
        }
    }

    void check_domain(domain_type const& domain) const
    {
        {
            spherical_single_type const* const _domain(
                dynamic_cast<spherical_single_type const*>(&domain));
            if (_domain)
            {
                check_domain(*_domain);
            }
        }
        {
            cylindrical_single_type const* const _domain(
                dynamic_cast<cylindrical_single_type const*>(&domain));
            if (_domain)
            {
                check_domain(*_domain);
            }
        }
        {
            spherical_pair_type const* const _domain(
                dynamic_cast<spherical_pair_type const*>(&domain));
            if (_domain)
            {
                check_domain(*_domain);
            }
        }
        {
            cylindrical_pair_type const* const _domain(
                dynamic_cast<cylindrical_pair_type const*>(&domain));
            if (_domain)
            {
                check_domain(*_domain);
            }
        }
        {
            multi_type const* const _domain(
                dynamic_cast<multi_type const*>(&domain));
            if (_domain)
            {
                check_domain(*_domain);
            }
        }
    }

    static rate_type calculate_k_tot(reaction_rules const& rules)
    {
        using namespace boost::lambda;
        using boost::lambda::_1;
        using boost::lambda::bind;
        rate_type k_tot(0.);
        std::for_each(boost::begin(rules), boost::end(rules),
            var(k_tot) += bind(&reaction_rule_type::k, _1));
        return k_tot;
    }

    reaction_rule_type const& draw_reaction_rule(reaction_rules const& rules)
    {
        const rate_type k_tot(calculate_k_tot(rules));
        const rate_type t(base_type::rng_.uniform(0., 1.) * k_tot);
        rate_type a(0.);
        BOOST_FOREACH(reaction_rule_type const& r, rules)
        {
            a += r.k();
            if (a > t)
                return r;
        }
        throw std::exception(); // should never happen
    }

    template<typename T1, typename T2>
    static position_type
    adjust_iv_with_old_iv(T1 const& new_iv, T2 const& old_iv)
    {
        length_type const angle(std::acos(old_iv[2] / length(old_iv)));
        if (std::fmod(angle, M_PI) != 0.0)
        {
            position_type const rotation_axis(
                normalize(position_type(-old_iv[1], old_iv[0], 0.)));
            return rotate_vector(new_iv, rotation_axis, angle);
        }
        else if (angle == 0.)
        {
            return new_iv;
        }
        else
        {
            return position_type(new_iv[0], new_iv[1], -new_iv[1]); 
        }
    }

    template<typename T>
    bool check_pair_pos(AnalyticalPair<traits_type, T> const& domain,
                        boost::array<particle_id_pair, 2> const& new_particles)
    {
        length_type const new_distance(
            (*base_type::world_).distance(new_particles[0].second.position(),
                                          new_particles[1].second.position()));
        length_type const r01(new_particles[0].second.radius() +
                              new_particles[1].second.radius());

        if (new_distance <= r01)
        {
            log_.warn(
                "rejected move: pair=%s, radii=%g, particle_distance=%g",
                boost::lexical_cast<std::string>(domain).c_str(),
                r01, new_distance);
            return false;
        }

        // particles within mobility radius.
        position_type const& com(shape(domain.shell().second).position());
        length_type const radius(shape(domain.shell().second).radius());
        length_type const d[2] = {
            (*base_type::world_).distance(com, new_particles[0].second.position()) + new_particles[0].second.radius(),
            (*base_type::world_).distance(com, new_particles[1].second.position()) + new_particles[1].second.radius()
        };
        if (d[0] > radius || d[1] > radius)
        {
            log_.warn(
                "rejected move: new particle(s) out of protective sphere: pair=%s, radii=%g, d0=%g, d1=%g",
                boost::lexical_cast<std::string>(domain).c_str(),
                d[0], d[1]);
            return false;
        }
        return true;
    }

    template<typename T>
    static PairGreensFunction* choose_pair_greens_function(
            AnalyticalPair<traits_type, T> const& domain,
            time_type const& t)
    {
        length_type const r0(domain.r0());
        length_type const distance_from_sigma(r0 - domain.sigma());
        length_type const distance_from_shell(domain.a_r() - r0);
        length_type const threshold_distance(
            traits_type::CUTOFF_FACTOR * std::sqrt(6. * domain.D_tot() * t));

        BOOST_ASSERT(::size(domain.reactions()) == 1);
        if (distance_from_sigma < threshold_distance)
        {
            if (distance_from_shell < threshold_distance)
            {
                // near both a and sigma;
                // use GreensFunction3DRadAbs
                LOG_DEBUG(("GF: normal"));
                return new GreensFunction3DRadAbs(
                    domain.D_tot(), domain.reactions()[0].k(),
                    r0, domain.sigma(), domain.a_r());
            }
            else
            {
                // near sigma; use GreensFunction3DRadInf
                LOG_DEBUG(("GF: only sigma"));
                return new GreensFunction3DRadInf(
                    domain.D_tot(), domain.reactions()[0].k(),
                    r0, domain.sigma());
            }
        }
        else
        {
            if (distance_from_shell < threshold_distance)
            {
                // near a;
                LOG_DEBUG(("GF: only a"));
                return new GreensFunction3DAbs(
                    domain.D_tot(), r0, domain.a_r());
            }
            else
            {
                // distant from both a and sigma; 
                LOG_DEBUG(("GF: free"));
                return new GreensFunction3D(domain.D_tot(), r0);
            }
        }
    }

    static length_type calculate_single_shell_size(
            single_type const& single,
            single_type const& closest,
            length_type const& distance,
            length_type const& shell_distance)
    {
        length_type const min_radius0(single.particle().second.radius());
        D_type const D0(single.particle().second.D());
        if (D0 == 0)
            return min_radius0;
        length_type const min_radius1(closest.particle().second.radius());
        D_type const D1(closest.particle().second.D());
        length_type const min_radius01(min_radius0 + min_radius1);
        length_type const sqrtD0(std::sqrt(D0));

        return std::max(std::min(sqrtD0 / (sqrtD0 + std::sqrt(D1))
                            * (distance - min_radius01) + min_radius0,
                            shell_distance / traits_type::SAFETY),
                        min_radius0);
    }

protected:
    int const num_retries_;
    length_type const user_max_shell_size_;

    domain_map domains_;
    spherical_shell_matrix_type ssmat_;
    cylindrical_shell_matrix_type csmat_;
    shell_matrix_map_type smatm_;
    shell_id_generator shidgen_;
    domain_id_generator didgen_;
    event_scheduler_type scheduler_;
    std::map<single_event_kind, int> single_step_count_;
    std::map<pair_event_kind, int> pair_step_count_;
    std::map<typename multi_type::event_kind, int> multi_step_count_;
    std::map<std::type_info const*, int> domain_count_per_type_;
    length_type single_shell_factor_;
    length_type multi_shell_factor_;
    unsigned int rejected_moves_;
    unsigned int zero_step_count_;
    static Logger& log_;
};

template<typename Ttraits>
inline char const* retrieve_domain_type_name(
    AnalyticalSingle<Ttraits,
        Shell<Sphere<typename Ttraits::world_traits_type::length_type>,
              typename Ttraits::domain_id_type> > const&)
{
    return "SphericalSingle";
}

template<typename Ttraits>
inline char const* retrieve_domain_type_name(
    AnalyticalSingle<Ttraits,
        Shell<Cylinder<typename Ttraits::world_traits_type::length_type>,
              typename Ttraits::domain_id_type> > const&)
{
    return "CylindricalSingle";
}

template<typename Ttraits>
inline char const* retrieve_domain_type_name(
    AnalyticalPair<Ttraits,
        Shell<Sphere<typename Ttraits::world_traits_type::length_type>,
              typename Ttraits::domain_id_type> > const&)
{
    return "SphericalPair";
}

template<typename Ttraits>
inline char const* retrieve_domain_type_name(
    AnalyticalPair<Ttraits,
        Shell<Cylinder<typename Ttraits::world_traits_type::length_type>,
              typename Ttraits::domain_id_type> > const&)
{
    return "CylindricalPair";
}



template<typename Ttraits_>
Logger& EGFRDSimulator<Ttraits_>::log_(Logger::get_logger("EGFRDSimulator"));

#endif /* EGFRDSIMULATOR_HPP */
